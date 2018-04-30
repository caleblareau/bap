import click
import os
import os.path
import sys
import shutil
import yaml
import random
import string
import itertools
import time
import pysam

from pkg_resources import get_distribution
from subprocess import call, check_call
from .bapHelp import *
from .bapProjectClass import *
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs

@click.command()
@click.version_option()

@click.argument('mode', type=click.Choice(['bam', 'c1fastq', 'check', 'support']))

@click.option('--input', '-i', help='Input for bap; varies by which mode is specified; see documentation')
@click.option('--output', '-o', default="bap_out", help='Output directory for analysis; see documentation.')

@click.option('--ncores', '-c', default = "detect", help='Number of cores to run the main job in parallel.')
@click.option('--reference-genome', '-r', default = "hg19", help='Support for built-in genome; choices are hg19, mm9, hg38, mm10, hg19_mm10_c and hg19-mm10 (species mix)')

@click.option('--cluster', default = "",  help='Message to send to Snakemake to execute jobs on cluster interface; see documentation.')
@click.option('--jobs', default = "0",  help='Max number of jobs to be running concurrently on the cluster interface.')

@click.option('--minimum-barcode-fragments', '-bf', default = 500, help='Minimum number of fragments to be thresholded for doublet merging.')
@click.option('--minimum-cell-fragments', '-cf', default = 500, help='Minimum number of unique to be thresholded for final output.')
@click.option('--minimum-jaccard-fragments', '-ji', default = 0.025, help='Minimum jaccard index for collapsing bead barcodes to cell barcodes')
@click.option('--one-to-one', '-oo', is_flag=True, help='Enforce that each bead barcode maps to one unique drop barcode (making this merging useless)')

@click.option('--extract-mito', '-em', is_flag=True, help='Extract mitochondrial DNA too?.')
@click.option('--keep-temp-files', '-z', is_flag=True, help='Keep all intermediate files.')

@click.option('--bedtools-genome', '-bg', default = "", help='Path to bedtools genome file; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--blacklist-file', '-bl', default = "", help='Path to bed file of blacklist; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--tss-file', '-ts', default = "", help='Path bed file of transcription start sites; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--mito-chromosome', '-mc', default = "default", help='Name of the mitochondrial chromosome; only necessary to specify if the reference genome is unknown but will overwrite default if needbe')
@click.option('--r-path', default = "", help='Path to R; by default, assumes that R is in PATH')

@click.option('--drop-tag', '-dt', default = "DB", help='New tag in the .bam file(s) that will be the name of the drop barcode.')
@click.option('--barcode-tag', '-bt', default = "XB", help='Tag in the .bam file(s) that point to the bead barcode; valid for bam mode.')
@click.option('--bam-name', '-bn', default="default", help='Name for the sam')

@click.option('--bwa-path', default = "", help='Path to bwa; by default, assumes that bwa is in PATH; only needed for "c1fastq" mode.')
@click.option('--bwa-index', '-bi', default = "", help='Path to the bwa (.fa) index; should be specified as if you were calling bwa (with file index prefix); only needed for "c1fastq" mode.')


def main(mode, input, output, ncores, reference_genome,
	cluster, jobs, minimum_barcode_fragments, minimum_cell_fragments, minimum_jaccard_fragments, one_to_one,
	extract_mito, keep_temp_files,
	bedtools_genome, blacklist_file, tss_file, mito_chromosome, r_path, 
	drop_tag, barcode_tag, bam_name,
	bwa_path, bwa_index):
	
	"""
	bap: Bead-based scATAC-seq data Processing \n
	Caleb Lareau, clareau <at> broadinstitute <dot> org \n
	
	mode = ['bam', 'c1fastq', 'check', 'support']\n
	"""
	
	__version__ = get_distribution('bap').version
	script_dir = os.path.dirname(os.path.realpath(__file__))

	click.echo(gettime() + "Starting bap pipeline v%s" % __version__)
	
	# Determine which genomes are available
	rawsg = os.popen('ls ' + script_dir + "/anno/bedtools/*.sizes").read().strip().split("\n")
	supported_genomes = [x.replace(script_dir + "/anno/bedtools/chrom_", "").replace(".sizes", "") for x in rawsg]  
	
	# Determine number of cores in main job
	if(ncores == "detect"):
		ncores = str(available_cpu_count())
	else:
		ncores = str(ncores)
	
	# Parameterize optional snakemake configuration
	snakeclust = ""
	njobs = int(jobs)
	if(njobs > 0 and cluster != ""):
		snakeclust = " --jobs " + str(jobs) + " --cluster '" + cluster + "' "
		click.echo(gettime() + "Recognized flags to process jobs on a cluster.")
	
	#------------
	# C1 Analysis
	#------------
	
	if(mode == "c1fastq"):
		click.echo(gettime() + "Preprocessing data as if it were from a C1 experiment...")
		
		# Figure out input
		samples, fastq1, fastq2 = inferSampleVectors(input)
		
		# Make output folders
		of = output; logs = of + "/logs"; fin = of + "/final"; trim = of + "/01_trimmed"; 
		aligned = of + "/02_aligned_reads"; processed = of + "/03_processed_reads"
		
		folders = [of, logs, fin, trim, aligned, processed, 
			of + "/.internal/parseltongue", of + "/.internal/samples",
			logs + "/bwa", logs + "/trim"]
	
		mkfolderout = [make_folder(x) for x in folders]

		# Create internal README files 
		if not os.path.exists(of + "/.internal/README"):
			with open(of + "/.internal/README" , 'w') as outfile:
				outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")
		if not os.path.exists(of + "/.internal/parseltongue/README"):	
			with open(of + "/.internal/parseltongue/README" , 'w') as outfile:
				outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")
		if not os.path.exists(of + "/.internal/samples/README"):
			with open(of + "/.internal" + "/samples" + "/README" , 'w') as outfile:
				outfile.write("This folder creates samples to be interpreted by Snakemake; don't modify it.\n\n")
			
		# Set up sample bam plain text file
		for i in range(len(samples)):
			with open(of + "/.internal/samples/" + samples[i] + ".fastqs.txt" , 'w') as outfile:
				outfile.write(fastq1[i] + "\t" + fastq2[i])
		
		# Set up dictionary with all the goodies
		d = {}
		d["bwa"] = get_software_path('bwa', bwa_path)
		d["bwa_index"] = bwa_index
		d["script_dir"] = script_dir
		d["barcode_tag"] = barcode_tag
		d["output"] = output
		
		y_s = of + "/.internal/parseltongue/c1fastq.object.yaml"
		with open(y_s, 'w') as yaml_file:
			yaml.dump(dict(d), yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
		
		# Trim, align, annotate, and merge via snakemake
		snakecmd_c1fastq = 'snakemake'+snakeclust+' --snakefile '+script_dir+'/bin/snake/Snakefile.bapc1fastq.c1fastq --cores '+ncores+' --config cfp="' + y_s + '" -T'
		os.system(snakecmd_c1fastq)
		
		sys.exit(gettime() + 'All done')
	
	if(mode == "support"):
		'''
		Show supported genomes and then bow out
		'''
		click.echo(gettime() + "List of built-in genomes supported in bap:")
		click.echo(gettime() + str(supported_genomes))
		sys.exit(gettime() + 'Specify one of these genomes or provide your own files (see documentation).')
		
	# Verify dependencies and set up an object to do all the dirty work
	p = bapProject(script_dir, supported_genomes, mode, input, output, ncores, reference_genome,
		cluster, jobs, minimum_barcode_fragments, minimum_cell_fragments, minimum_jaccard_fragments, one_to_one,
		extract_mito, keep_temp_files,
		bedtools_genome, blacklist_file, tss_file, mito_chromosome, r_path, 
		drop_tag, barcode_tag, bam_name,
		bwa_path, bwa_index)
	
	if(reference_genome in ["hg19-mm10", "hg19_mm10_c"]):
		speciesMix = True
	else:
		speciesMix = False
	
	# Make a counts table from user-supplied peaks and bam files
	if(mode == 'bam'):
	
		click.echo(gettime() + "Attempting to parse supplied .bam files for analysis.")
	
		# Make sure that the .bam file is valid
		if(not os.path.exists(input)):
			sys.exit("Cannot parse supplied .bam file in --input")
		if(not os.path.exists(input + ".bai")):
			sys.exit("Index supplied .bam before proceeding")
		
		# Make output folders
		of = output; logs = of + "/logs"; fin = of + "/final"; mito = of + "/mito"; temp = of + "/temp"
		temp_filt_split = temp + "/filt_split"; temp_frag_overlap = temp + "/frag_overlap"
		temp_drop_barcode = temp + "/drop_barcode"
		
		folders = [of, logs, fin, mito, temp,
			temp_filt_split, temp_frag_overlap,temp_drop_barcode, 
			of + "/.internal/parseltongue", of + "/.internal/samples"]
	
		mkfolderout = [make_folder(x) for x in folders]
		
		# Create internal README files 
		if not os.path.exists(of + "/.internal/README"):
			with open(of + "/.internal/README" , 'w') as outfile:
				outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")
		if not os.path.exists(of + "/.internal/parseltongue/README"):	
			with open(of + "/.internal/parseltongue/README" , 'w') as outfile:
				outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")
		if not os.path.exists(of + "/.internal/samples/README"):
			with open(of + "/.internal" + "/samples" + "/README" , 'w') as outfile:
				outfile.write("This folder creates samples to be interpreted by Snakemake; don't modify it.\n\n")
		
		#-------------------------------------------------------
		# Step 1- Filter and split input .bam file by chromosome
		#-------------------------------------------------------
		line1 = 'python ' +script_dir+'/bin/python/11_quantBarcode_Filt.py --input '+p.bamfile
		line2 = ' --name ' + p.bam_name + ' --output ' + temp_filt_split + ' --barcode-tag ' 
		line3 = p.barcode_tag + ' --min-fragments ' + str(p.minimum_barcode_fragments)
		line4 = " --bedtools-genome " +p.bedtoolsGenomeFile + " --ncores " + str(ncores)
			
		filt_split_cmd = line1 + line2 + line3 + line4
		os.system(filt_split_cmd)
		
		#-------------------------------------------
		# Step 2 / 3- Process fragments by Snakemake
		#-------------------------------------------
		
		# Round trip the .yaml of user configuration
		y_s = of + "/.internal/parseltongue/bap.object.bam.yaml"
		with open(y_s, 'w') as yaml_file:
			yaml.dump(dict(p), yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
			
		snakecmd_chr = 'snakemake'+snakeclust+' --snakefile '+script_dir+'/bin/snake/Snakefile.bap.chr --cores '+ncores+' --config cfp="' + y_s + '" -T'
		os.system(snakecmd_chr)

		#-----------------------------------
		# Step 4 - QC stats / outside Snakemake since not-essential
		#-----------------------------------
		barcodeTranslateFile = p.output + "/final/" + p.bam_name + ".barcodeTranslate.tsv"
		barcodeQuantFile =  p.output + "/final/" + p.bam_name + ".barcodequants.csv"
		finalBamFile = p.output + "/final/" + p.bam_name + ".bap.bam"
		
		if(speciesMix):
			qc_R = script_dir + "/bin/R/16a_qualityControlReport_SM.R"
		else:
			qc_R = script_dir + "/bin/R/16_qualityControlReport.R"

		r_callQC = " ".join([p.R+"script", qc_R, finalBamFile, barcodeTranslateFile, barcodeQuantFile, p.tssFile, p.drop_tag, p.blacklistFile])
		print(r_callQC)
		os.system(r_callQC)

		#-----------------------------------
		# Step 5 - Process mitochondria
		#-----------------------------------
		if(extract_mito):
			click.echo(gettime() + "Creating a new .bam file for mitochondria.")
			
			dict_file = fin+"/"+p.bam_name+".barcodeTranslate.tsv"
			
			line1 = 'python ' +script_dir+'/bin/python/15_processMito.py --input '+p.bamfile
			line2 = ' --output ' + mito + "/" + p.bam_name + ".mito.bam" + " --mitochr " + p.mitochr 
			line3 = ' --bead-barcode ' + p.barcode_tag + ' --drop-barcode ' + p.drop_tag + " --dict-file " + dict_file
			mito_cmd = line1 + line2 + line3
			
			os.system(mito_cmd)
		
		#-------------------------------------------------------
		# Final-- remove intermediate files if necessary
		#-------------------------------------------------------
		if keep_temp_files:
			click.echo(gettime() + "Temporary files not deleted since --keep-temp-files was specified.")
		else:
			byefolder = of
			shutil.rmtree(byefolder + "/.internal")
			shutil.rmtree(byefolder + "/temp")

			if(not extract_mito):
				shutil.rmtree(byefolder + "/mito")
			
			click.echo(gettime() + "Intermediate files successfully removed.")
		
			
	# If we get to this stage in the processing, then everything worked
	if (mode == "check"):
		click.echo(gettime() + "Dependencies and user-reported file paths OK")
	
	click.echo(gettime() + "Complete.")
	

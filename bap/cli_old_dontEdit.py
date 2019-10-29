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
from .bap1ProjectClass import *
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs

@click.command()
@click.version_option()

@click.argument('mode', type=click.Choice(['bam', 'check', 'support']))

@click.option('--input', '-i', help='Input for bap; varies by which mode is specified; most certainly a .bam file with an index')
@click.option('--output', '-o', default="bap_out", help='Output directory for analysis; this is where everything is housed.')
@click.option('--name', '-n', default="default", help='Name for all of the output files (default: uses the .bam prefix)')

@click.option('--ncores', '-c', default = "detect", help='Number of cores to run the main job in parallel.')
@click.option('--reference-genome', '-r', default = "hg19", help='Support for built-in genome; choices are hg19, mm9, hg38, mm10, hg19_mm10_c and hg19-mm10 (species mix)')

@click.option('--cluster', default = "",  help='Message to send to Snakemake to execute jobs on cluster interface; see documentation.')
@click.option('--jobs', default = "0",  help='Max number of jobs to be running concurrently on the cluster interface.')
@click.option('--peak-file', '-pf', default = "", help='If supplied, compute FRIP (in QC stats) and generate Summarized Experiment')

@click.option('--minimum-barcode-fragments', '-bf', default = 0, help='Minimum number of fragments to be thresholded for doublet merging.')
@click.option('--minimum-cell-fragments', '-cf', default = 500, help='Minimum number of unique to be thresholded for final output.')
@click.option('--barcode-whitelist', '-wl', default = "", help='File path of a whitelist of bead barcodes (one per line) to be used in lieu of a fixed threshold.')

@click.option('--minimum-jaccard-index', '-ji', default = 0.0, help='Minimum jaccard index for collapsing bead barcodes to cell barcodes')
@click.option('--nc-threshold', '-nc', default = 6, help='Number of barcodes that a paired-end read must be observed for the read to be filtered.')
@click.option('--one-to-one', '-oo', is_flag=True, help='Enforce that each bead barcode maps to one unique drop barcode (making this merging useless)')
@click.option('--barcoded-tn5',  is_flag=True, help='Process data knowing that the barcodes were generated with a barcoded Tn5')

@click.option('--extract-mito', '-em', is_flag=True, help='Extract mitochondrial DNA too?.')
@click.option('--keep-temp-files', '-z', is_flag=True, help='Keep all intermediate files.')
@click.option('--mapq', '-mq', default = 30, help='Minimum mapping quality to keep read for downstream analyses')

@click.option('--bedtools-genome', '-bg', default = "", help='Path to bedtools genome file; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--blacklist-file', '-bl', default = "", help='Path to bed file of blacklist; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--tss-file', '-ts', default = "", help='Path bed file of transcription start sites; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--mito-chromosome', '-mc', default = "default", help='Name of the mitochondrial chromosome; only necessary to specify if the reference genome is unknown but will overwrite default if needbe')
@click.option('--r-path', default = "", help='Path to R; by default, assumes that R is in PATH')

@click.option('--drop-tag', '-dt', default = "DB", help='New tag in the .bam file(s) that will be the name of the drop barcode.')
@click.option('--bead-tag', '-bt', default = "XB", help='Tag in the .bam file(s) that point to the bead barcode; valid for bam mode.')


def main(mode, input, output, name, ncores, reference_genome,
	cluster, jobs, peak_file,
	minimum_barcode_fragments, minimum_cell_fragments, barcode_whitelist,
	minimum_jaccard_index, nc_threshold, one_to_one, barcoded_tn5,
	extract_mito, keep_temp_files, mapq, 
	bedtools_genome, blacklist_file, tss_file, mito_chromosome, r_path, 
	drop_tag, bead_tag):
	
	"""
	bap: Bead-based scATAC-seq data Processing \n
	Caleb Lareau, clareau <at> broadinstitute <dot> org \n
	
	mode = ['bam', 'check', 'support']\n
	"""
	
	__version__ = get_distribution('bap-atac').version
	script_dir = os.path.dirname(os.path.realpath(__file__))

	click.echo(gettime() + "Starting bap pipeline v%s" % __version__)
	click.echo(gettime() + "WARNING: this is bap1-- you should really use bap2 unless you are trying to reproduce old results!")
	click.echo(gettime() + "NOTE: The same arguments apply-- see `bap2 --help`")

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
	

	if(mode == "support"):
		'''
		Show supported genomes and then bow out
		'''
		click.echo(gettime() + "List of built-in genomes supported in bap:")
		click.echo(gettime() + str(supported_genomes))
		sys.exit(gettime() + 'Specify one of these genomes or provide your own files (see documentation).')
		
	# Verify dependencies and set up an object to do all the dirty work
	p = bap1Project(script_dir, supported_genomes, mode, input, output, name, ncores, reference_genome,
		cluster, jobs, peak_file,
		minimum_barcode_fragments, minimum_cell_fragments, barcode_whitelist,
		minimum_jaccard_index, nc_threshold, one_to_one, barcoded_tn5, 
		extract_mito, keep_temp_files, mapq, 
		bedtools_genome, blacklist_file, tss_file, mito_chromosome, r_path, 
		drop_tag, bead_tag)
	
	if(reference_genome in ["hg19-mm10", "hg19_mm10_c", "hg19-mm10_nochr"]):
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
		
		# Verify that the supplied reference genome and the bam have overlapping chromosomes
		idxstats = pysam.idxstats(input)
		chrs_bam = [x.split("\t")[0] for x in idxstats.split("\n")]
		chrs_ref = []
		with open(p.bedtoolsGenomeFile) as f:
			for line in f:
				chrs_ref.append(line.split("\t")[0])
		
		# Quant the overlaps and exit if no overlap
		chrs_both = intersection(chrs_bam, chrs_ref)
		if(len(chrs_both) > 0):
			click.echo(gettime() + "Found " + str(len(chrs_both)) + " chromosomes for analysis (including mitochondria).")
		else:
			sys.exit("Found no overlapping chromosomes between bam and reference. Check reference genome specification with the --reference-genome flag.")
					
		# Make output folders
		of = output; logs = of + "/logs"; fin = of + "/final"; mito = of + "/mito"; temp = of + "/temp"
		temp_filt_split = temp + "/filt_split"; temp_frag_overlap = temp + "/frag_overlap"; knee = of + "/knee"
		temp_drop_barcode = temp + "/drop_barcode"
		
		folders = [of, logs, fin, mito, temp,
			temp_filt_split, temp_frag_overlap,temp_drop_barcode, knee,
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
		line1 = 'python ' +script_dir+'/bin/python/10_quantBarcode_Filt.py --input '+p.bamfile
		line2 = ' --name ' + p.name + ' --output ' + temp_filt_split + ' --barcode-tag ' 
		line3 = p.bead_tag + ' --min-fragments ' + str(p.minimum_barcode_fragments)
		line4 = " --bedtools-genome " +p.bedtoolsGenomeFile + " --ncores " + str(ncores) + " --mapq " + str(mapq)
		line5 = " --barcode-whitelist " + p.barcode_whitelist + " --knee-call " + script_dir+'/bin/R/10b_knee_execute.R'
			
		filt_split_cmd = line1 + line2 + line3 + line4 + line5
		print(filt_split_cmd)
		os.system(filt_split_cmd)
		
		# Verify that we got output and fail if not
		val = file_len(p.output + "/final/" + p.name + ".barcodequants.csv") -1
		click.echo(gettime() + "Found " + str(val) + " abundant barcodes passing filter.")
		
		if(val < 2):
			sys.exit(gettime() + "Found insufficient # of barcodes for analyis; check reference genome or lower --minimum-barcode-fragments")
		
		#----------------------------------------
		# Step 2 - Process fragments by Snakemake
		#----------------------------------------
		
		# Round trip the .yaml of user configuration
		y_s = of + "/.internal/parseltongue/bap.object.bam.yaml"
		with open(y_s, 'w') as yaml_file:
			yaml.dump(dict(p), yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
			
		snakecmd_chr = 'snakemake'+snakeclust+' --snakefile '+script_dir+'/bin/snake/Snakefile.bap.chr --cores '+ncores+' --config cfp="' + y_s + '"'
		os.system(snakecmd_chr)
		
		#-------------------------------------------------
		# Step 3 - Handle aggregation of NC filtered reads
		#-------------------------------------------------
		nc_R = script_dir + "/bin/R/15_ncReadCleanup.R"
		barcodeQuantFile =  p.output + "/final/" + p.name + ".barcodequants.csv"
		
		# Irrespective, modify the barcode quantification sumstats
		nc_R_call = " ".join([p.R+"script", nc_R, p.output, str(keep_temp_files), barcodeQuantFile])
		os.system(nc_R_call)

		#-----------------------------------
		# Step 4 - QC stats / outside Snakemake since not-essential
		#-----------------------------------
		barcodeTranslateFile = p.output + "/final/" + p.name + ".barcodeTranslate.tsv"
		qcStats16File =  p.output + "/final/" + p.name + ".QCstats.csv"
		finalBamFile = p.output + "/final/" + p.name + ".bap.bam"
		qc_R = script_dir + "/bin/R/16_qualityControlReport_SE.R"
		
		# Determine user flags for proper QC stuff		
		if(p.peakFile == ""):
			peakFileGo = "none"
		else:
			peakFileGo = p.peakFile
		
		if(speciesMix):
			speciesMixGo = "yes"
		else:
			speciesMixGo = "no"
			
		if(one_to_one):
			oneToOneGo = "yes"
		else:
			oneToOneGo = "no"
			
		click.echo(gettime() + "generating QC report + summarized experiment file...")
		# Full system call to R script
		r_callQC = " ".join([p.R+"script", qc_R, finalBamFile, barcodeTranslateFile, barcodeQuantFile, p.tssFile, p.drop_tag, p.blacklistFile, peakFileGo, speciesMixGo, oneToOneGo])
		os.system(r_callQC)

		#-----------------------------------
		# Step 5 - Process mitochondria
		#-----------------------------------
		if(extract_mito):
			click.echo(gettime() + "Creating a new .bam file for mitochondria.")
			
			dict_file = fin+"/"+p.name+".barcodeTranslate.tsv"
			
			line1 = 'python ' +script_dir+'/bin/python/17_processMito.py --input '+p.bamfile
			line2 = ' --output ' + mito + "/" + p.name + ".mito.bam" + " --mitochr " + p.mitochr 
			line3 = ' --bead-barcode ' + p.bead_tag + ' --drop-barcode ' + p.drop_tag + " --dict-file " + dict_file
			mito_cmd = line1 + line2 + line3
			
			os.system(mito_cmd)
		
		#-----------------------------------
		# Step 6 - Make knee plots
		#-----------------------------------
		bapParamsFile =  p.output + "/knee/" + p.name + ".bapParams.csv"
		beadBarcodesFile = p.output + "/knee/" + p.name + ".barcodeQuantSimple.csv"
		implicatedBarcodeFile = p.output + "/final/" + p.name + ".implicatedBarcodes.csv.gz"
		
		click.echo(gettime() + "generating knee plots associated with beads.")
		kneePlot_R = script_dir + "/bin/R/19_makeKneePlots.R"
		r_callKneePlot = " ".join([p.R+"script", kneePlot_R, bapParamsFile, beadBarcodesFile, implicatedBarcodeFile])
		os.system(r_callKneePlot)
		
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
	

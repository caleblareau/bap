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

@click.argument('mode', type=click.Choice(['bam', 'support']))

@click.option('--input', '-i', default = ".", required=True, help='Input for bap; varies by which mode is specified; see documentation')
@click.option('--output', '-o', default="bap_out", help='Output directory for analysis; see documentation.')
@click.option('--name', '-n', default="bap",  help='Prefix for project name')
@click.option('--ncores', '-c', default = "detect", help='Number of cores to run the main job in parallel.')

@click.option('--bowtie2-index', '-bi', default = "", required=True, help='Path to the bowtie2 index; should be specified as if you were calling bowtie2 (with file index prefix)')
@click.option('--cluster', default = "",  help='Message to send to Snakemake to execute jobs on cluster interface; see documentation.')
@click.option('--jobs', default = "0",  help='Max number of jobs to be running concurrently on the cluster interface.')

@click.option('--peaks-file', '-pf', default = "", help='Path to a pre-defined peaks file; option is only useful in `counts` mode')
@click.option('--by-rgid', is_flag = True, help='Option to makes counts table RGID aware; option is only useful in `counts` mode')

@click.option('--peak-width', '-pw', default = "250", help='Fixed width value of resulting peaks from summit calling / padding. 250, 500 recommended.')
@click.option('--keep-duplicates', '-kd', is_flag=True, help='Keep optical/PCR duplicates')
@click.option('--extract-mito', '-em', is_flag=True, help='Extract mitochondrial DNA too?.')
@click.option('--reference-genome', '-rg', default = "", help='Support for built-in genome; choices are hg19, mm9, hg38, mm10, hg19_mm10_c (species mix)')

@click.option('--keep-temp-files', '-z', is_flag=True, help='Keep all intermediate files.')
@click.option('--skip-fastqc', '-sf', is_flag=True, help='Throw this flag to skip fastqc on the trimmed fastqs; will only run if software is discovered in the path.')
@click.option('--overwrite', '-ov', is_flag=True, help='Throw this flag to always retrim, realign, reprocess, etc. each sample even if it already exists in output (see documentation).')

@click.option('--bedtools-genome', '-bg', default = "", help='Path to bedtools genome file; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--blacklist-file', '-bl', default = "", help='Path to bed file of blacklist; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--tss-file', '-ts', default = "", help='Path bed file of transcription start sites; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--macs2_genome_size', '-mg', default = "", help='String passed to macs2 for ; overrides default if --reference-genome flag is set and is necessary for non-supported genomes.')
@click.option('--bs-genome', '-bs', default = "", help='String corresponding to the R/Bioconductor package for BS genome of build; overrides default if --reference-genome flag is set and is necessary for non-supported genomes..')

@click.option('--bedtools-path', default = "", help='Path to bedtools; by default, assumes that bedtools is in PATH')
@click.option('--bowtie2-path', default = "", help='Path to bowtie2; by default, assumes that bowtie2 is in PATH')
@click.option('--macs2-path', default = "", help='Path to macs2; by default, assumes that macs2 is in PATH')
@click.option('--samtools-path', default = "", help='Path to samtools; by default, assumes that samtools is in PATH')
@click.option('--R-path', default = "", help='Path to R; by default, assumes that R is in PATH')

def main(mode, input, output, name, ncores, bowtie2_index,
	cluster, jobs, peaks_file, by_rgid,
	peak_width, keep_duplicates, max_javamem, trash_mito, reference_genome,
	clipl, clipr, py_trim, keep_temp_files, skip_fastqc, overwrite,
	bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome, 
	bedtools_path, bowtie2_path, java_path, macs2_path, samtools_path, r_path):
	
	"""
	bap: a toolkit for PROcessing ATAC-seq data. \n
	Caleb Lareau, clareau <at> broadinstitute <dot> org
	
	modes = ['bam', 'support']\n
	"""
	
	__version__ = get_distribution('bap').version
	script_dir = os.path.dirname(os.path.realpath(__file__))

	click.echo(gettime() + "Starting bap pipeline v%s" % __version__)
	
	# Determine which genomes are available
	rawsg = os.popen('ls ' + script_dir + "/anno/bedtools/*.sizes").read().strip().split("\n")
	supported_genomes = [x.replace(script_dir + "/anno/bedtools/chrom_", "").replace(".sizes", "") for x in rawsg]  

	if(mode == "support"):
		click.echo(gettime() + "List of built-in genomes supported in bap:")
		click.echo(gettime() + str(supported_genomes))
		sys.exit(gettime() + 'Specify one of these genomes or provide your own files (see documentation).')
		
	
	# Take a collection of summits files and return a consensus set of peaks	
	if(mode == 'summitsToPeaks'):
		click.echo(gettime() + "Starting inference of peaks from summits.")
		
		# Need chromosome sizes and blacklist
		bedtoolsGenomeFile, blacklistFile = getBfiles(bedtools_genome, blacklist_file, reference_genome, script_dir, supported_genomes)
		
		# Figure out which samples to process
		bedFiles = os.popen("ls " + input.rstrip("/") + "/*summits.bed*").read().strip().split("\n")
		if(len(bedFiles) < 1):
			sys.exist("No summit *summits.bed* files found; QUITTING")
		else:
			click.echo(gettime() + "Calling peaks from these samples:")
			click.echo(gettime() + str(bedFiles))
		
		# Verify dependencies
		R = get_software_path('R', r_path)
		check_R_packages(['data.table', 'GenomicRanges', 'tools'], R)
		
		# Execute software
		make_folder(output)
		summitRcall = " ".join([R +"script", script_dir + "/bin/R/summitsToCleanPeaks.R", ",".join(bedFiles), peak_width, blacklistFile, bedtoolsGenomeFile, str(999999999), str(0.01), output, name])
		os.system(summitRcall)	
		click.echo(gettime() + "Completed peak inference from summit files.")
		sys.exit()
	
	# TO DO: 
	# Make a mode to handle split-pool data
	if(mode == 'indexSplit'):
		sys.exit("Mode does not actually work yet")
	
	# Make a counts table from user-supplied peaks and bam files
	if(mode == 'counts'):
		click.echo(gettime() + "Attempting to assemble counts table from user-specified input.")
		
		# Verify dependencies
		R = get_software_path('R', r_path)
		check_R_packages(['chromVAR', 'SummarizedExperiment', 'tools'], R)
		
		# Make sure that there are samples to process / there is a peak file
		bamfiles = os.popen("ls " + input.rstrip("/") + "/*.bam").read().strip().split("\n")
		if(len(bamfiles) < 1):
			sys.exist("No sample *.bam files found in user-specified input; QUITTING")
		else:
			click.echo(gettime() + "Making a counts table from these samples:")
			click.echo(gettime() + str(bamfiles))
		if(os.path.isfile(peaks_file)):
			click.echo(gettime() + "Found peaks file: " + peaks_file)
		
		# Execute software
		make_folder(output)
		countsRcall = " ".join([R +"script", script_dir + "/bin/R/makeCountsTable.R", input, peaks_file, str(by_rgid), output, name])
		os.system(countsRcall)	
		click.echo(gettime() + "Completed peak inference from summit files.")
		sys.exit()
	
	p = bapProject(script_dir, supported_genomes, mode, input, output, name, ncores, bowtie2_index,
		cluster, jobs, peak_width, keep_duplicates, max_javamem, trash_mito, reference_genome,
		clipl, clipr, py_trim, keep_temp_files, skip_fastqc, overwrite,
		bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome, 
		bedtools_path, bowtie2_path, java_path, macs2_path, samtools_path, r_path)
	
	if (mode == "check"):
		click.echo(gettime() + "Dependencies and user-reported file paths OK")
		click.echo("\nbap will process the following samples / files with bulk / single specified: \n")
		print("Sample", "Fastq1", "Fastq2")
		for x in range(len(p.samples)):
			print(p.samples[x], p.fastq1[x], p.fastq2[x])
		click.echo("\nIf this table doesn't look right, consider specifying a manually created sample input table (see documentation).\n")
		sys.exit(gettime() + "Successful check complete; QUITTING.")
	
	# Single or bulk processing
	if(mode == "single" or mode == "bulk"):

		# Potentially submit jobs to cluster		
		if(ncores == "detect"):
			ncores = str(available_cpu_count())
		else:
			ncores = str(ncores)
		
		snakeclust = ""
		njobs = int(jobs)
		if(njobs > 0 and cluster != ""):
			snakeclust = " --jobs " + jobs + " --cluster '" + cluster + "' "
			click.echo(gettime() + "Recognized flags to process jobs on a computing cluster.")		
		
		# Make output folders
		of = output; logs = of + "/logs"; fin = of + "/final"; trim = of + "/01_trimmed"; 
		aligned = of + "/02_aligned_reads"; processed = of + "/03_processed_reads";
		qc = of + "/04_qc"
		
		folders = [of, logs, fin, trim, aligned, processed, qc,
			of + "/.internal/parseltongue", of + "/.internal/samples",
			logs + "/bowtie2", logs + "/trim", logs + "/macs2",
			of + "/03_processed_reads/temp", fin + "/plots"]
	
		mkfolderout = [make_folder(x) for x in folders]
		
		make_folder(logs + "/picard")
		make_folder(logs + "/picard/inserts")
		make_folder(logs + "/tss")
		make_folder(logs + "/samples")
		make_folder(of + "/mito")
		
		if not keep_duplicates:
			make_folder(logs + "/picard/markdups")
		if not skip_fastqc:
			make_folder(logs + "/fastqc")
			
		if(mode == "bulk"):
			make_folder(of + "/final/bams")
			make_folder(of + "/final/summits")
			make_folder(of + "/04_qc/macs2_each")
		if(mode == "single"):
			make_folder(of + "/03_processed_reads/bams")

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
		
		# Create promoter file:
		ptss = of + "/.internal/promoter.tss.bed"
		if not os.path.exists(ptss):
			os.system('''awk '{print $1"\t"$2-2000"\t"$3+2000"\t"$4}' '''+ p.tssFile + " > " + ptss)
			
		# Set up sample bam plain text file
		for i in range(len(p.samples)):
			with open(of + "/.internal/samples/" + p.samples[i] + ".fastqs.txt" , 'w') as outfile:
				outfile.write(p.fastq1[i] + "\t" + p.fastq2[i])
		
		y_s = of + "/.internal/parseltongue/bap.object.yaml"
		with open(y_s, 'w') as yaml_file:
			yaml.dump(dict(p), yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
		
		snakecmd_scatter = 'snakemake'+snakeclust+' --snakefile '+script_dir+'/bin/snake/Snakefile.bap.scatter --cores '+ncores+' --config cfp="' + y_s + '" -T'
		os.system(snakecmd_scatter)
		
		if(mode == 'single'):
			
			# Merge into one .bam file:
			finalmergedbam = fin + "/"+p.name+".merged.bam"
			if not os.path.isfile(finalmergedbam):
				os.system(p.samtools + " merge " +finalmergedbam+" "+ of + "/03_processed_reads/bams/*.bam")
				pysam.index(finalmergedbam)
			
		snakecmd_gather = 'snakemake --snakefile '+script_dir+'/bin/snake/Snakefile.bap.gather --cores '+ncores+' --config cfp="' + y_s + '" -T'
		os.system(snakecmd_gather)
		
		if keep_temp_files:
			click.echo(gettime() + "Temporary files not deleted since --keep-temp-files was specified.")
		else:
			if(mode == "bulk" or mode == "single"):
				byefolder = of
			
			shutil.rmtree(byefolder + "/.internal")
			shutil.rmtree(byefolder + "/01_trimmed")
			shutil.rmtree(byefolder + "/02_aligned_reads")
			shutil.rmtree(byefolder + "/03_processed_reads")
			shutil.rmtree(byefolder + "/04_qc")
			if(trash_mito):
				shutil.rmtree(byefolder + "/mito")
			
			click.echo(gettime() + "Intermediate files successfully removed.")
		
	click.echo(gettime() + "Complete.")
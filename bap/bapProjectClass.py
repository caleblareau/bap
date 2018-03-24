import click
import os
import os.path
import sys
import shutil
import random
import string
import itertools
import time
import platform
from ruamel import yaml
from .bapHelp import *

def getBfiles(bedtools_genome, blacklist_file, reference_genome, script_dir, supported_genomes):

	'''
	Function that isn't actually a bapProject specific function.
	Used to collate the built-in genomes with the possibility that the
	user specified another genome.
	'''
	
	# Handle bedtools
	if(bedtools_genome == "" and reference_genome == ""):
		sys.exit("ERROR: bap needs either the bedtools genome or a correctly specified reference genome to get peaks from summit files; QUITTING")
	elif any(reference_genome == s for s in supported_genomes):
		bedtoolsGenomeFile = script_dir + "/anno/bedtools/chrom_" + reference_genome + ".sizes"
	else:
		if(os.path.isfile(bedtools_genome)):
			bedtoolsGenomeFile = bedtools_genome
		else: 
			sys.exit("Could not find the bedtools genome file: %s" % bedtools_genome)
	
	# Handle blacklist	
	if(blacklist_file == "" and reference_genome == ""):
		sys.exit("ERROR: bap needs either a blacklist bed file or a correctly specified reference genome to get peaks from summit files; QUITTING")
	elif any(reference_genome == s for s in supported_genomes):
		blacklistFile = script_dir + "/anno/blacklist/" + reference_genome + ".full.blacklist.bed"
	else:
		if(os.path.isfile(bedtools_genome)):
			blacklistFile = bedtools_genome
		else: 
			sys.exit("Could not find the blacklist file: %s" % bedtools_genome)
			
	return(bedtoolsGenomeFile, blacklistFile)

class bapProject():
	def __init__(self, script_dir, supported_genomes, mode, input, output, name, ncores, bowtie2_index,
		cluster, jobs, peak_width, keep_duplicates, max_javamem, trash_mito, reference_genome,
		clipl, clipr, py_trim, keep_temp_files, skip_fastqc, overwrite,
		bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome, 
		bedtools_path, bowtie2_path, java_path, macs2_path, samtools_path, r_path):

		# Figure out operating system just for funzies; not presently used
		self.os = "linux"
		if(platform.platform()[0:5]=="Darwi"):
			self.os = "mac"
	
		#----------------------------------
		# fastq processing specific options
		#----------------------------------
		if(mode == "fastq"):	
			
			# Need to align with bowtie2
			self.bowtie2 = get_software_path('bowtie2', bowtie2_path)
			
			# verify bowtie2 index
			bwt2idxfiles = os.popen("ls " + bowtie2_index + "*.bt2*").read().strip().split("\n")
			if(len(bwt2idxfiles) < 6):
				sys.exit("ERROR: cannot find bowtie2 index; specify with --bowtie2-index and make sure to add the prefix along with the folder path")
			else:
				self.bowtie2_index = bowtie2_index
				
		#-----------------------------------------
		# Verify R and all of its goodies are here
		#-----------------------------------------
		R = get_software_path('R', r_path)
		check_R_packages(['Rsamtools', 'GenomicAlignments', 'GenomicRanges', 'BiocParallel', 'dplyr', 'SummarizedExperiment'], R)
		self.R = R
		
		#----------------------------------
		# Assign straightforward attributes
		#----------------------------------
		self.extract_mito = extract_mito
		self.cluster = cluster
		self.jobs = jobs
		self.output = output
		self.mode = mode
		self.script_dir = script_dir
		
		# Collect samples / fastq lists
		self.samples, self.fastq1, self.fastq2 = inferSampleVectors(input)
		
		# remove samples that have been previously processed
		if not overwrite:
			self.samples, self.fastq1, self.fastq2 = filterExistingSamples(self.samples, self.fastq1, self.fastq2, output)
		
		#------------------------
		# Handle reference genome
		#------------------------
		self.reference_genome = reference_genome
		if any(self.reference_genome == s for s in supported_genomes):
			click.echo(gettime() + "Found designated reference genome: %s" % self.reference_genome)
			
			self.tssFile = script_dir + "/anno/TSS/" + self.reference_genome + ".refGene.TSS.bed"
			self.blacklistFile = script_dir + "/anno/blacklist/" + self.reference_genome + ".full.blacklist.bed"
			self.bedtoolsGenomeFile = script_dir + "/anno/bedtools/chrom_" + self.reference_genome + ".sizes"

		else: 
			click.echo(gettime() + "Could not identify this reference genome: %s" % self.reference_genome)
			click.echo(gettime() + "Attempting to infer necessary input files from user specification.")
			necessary = [bedtools_genome, blacklist_file, tss_file, macs2_genome_size, bs_genome]
			if '' in necessary:
				if reference_genome == '':
					sys.exit("ERROR: specify valid reference genome with --reference-genome flag; QUITTING")
				else:
					sys.exit("ERROR: non-supported reference genome specified so these five must be validly specified: --bedtools-genome, --blacklist-file, --tss-file; QUITTING")
		
		#------------------------------		
		# Make sure all files are valid
		#------------------------------	
		if(bedtools_genome != ""):
			if(os.path.isfile(bedtools_genome)):
				self.bedtoolsGenomeFile = bedtools_genome
			else: 
				sys.exit("Could not find the bedtools genome file: %s" % bedtools_genome)
				
		if(blacklist_file != ""):
			if(os.path.isfile(blacklist_file)):
				self.blacklistFile = blacklist_file
			else: 
				sys.exit("Could not find the blacklist bed file: %s" % blacklist_file)
		
		if(tss_file != ""):	
			if(os.path.isfile(tss_file)):
				self.tssFile = tss_file
			else: 
				sys.exit("Could not find the transcription start sites file: %s" % tss_file)		
		
	#--------------------------------------------------------------------------------
	# Define a method to dump the object as a .yaml/dictionary for use in other files
	#--------------------------------------------------------------------------------
	def __iter__(self):
		
		# Purposefully skip samples, fastq1, fastq2 -- will put individual samples there in call
		yield 'script_dir', self.script_dir
		yield 'tssFile', self.tssFile
		yield 'blacklistFile', self.blacklistFile
		yield 'bedtoolsGenomeFile', self.bedtoolsGenomeFile
		yield 'BSgenome', self.BSgenome
		yield 'macs2_genome_size', self.macs2_genome_size
		yield 'clipl', self.clipl
		yield 'clipr', self.clipr
		yield 'peak_width', self.peak_width
		yield 'max_javamem', self.max_javamem
		yield 'trash_mito', self.trash_mito
		yield 'keep_duplicates', self.keep_duplicates
		yield 'cluster', self.cluster
		yield 'jobs', self.jobs
		yield 'name', self.name
		yield 'output', self.output
		yield 'mode', self.mode
		yield 'skip_fastqc', self.skip_fastqc
		yield 'bedtools', self.bedtools
		yield 'bowtie2', self.bowtie2
		yield 'bowtie2_index', self.bowtie2_index
		yield 'java', self.java
		yield 'macs2', self.macs2
		yield 'samtools', self.samtools
		yield 'R', self.R
		
	
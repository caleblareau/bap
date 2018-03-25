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
	def __init__(self, script_dir, supported_genomes, mode, input, output, ncores, reference_genome,
		cluster, jobs, minimum_barcode_fragments, minimum_cell_fragments,
		extract_mito, keep_temp_files,
		bedtools_genome, blacklist_file, tss_file, r_path, 
		barcode_tag, bam_name,
		bowtie2_path, bowtie2_index):
		
				
		#----------------------------------
		# Assign straightforward attributes
		#----------------------------------
		self.script_dir = script_dir
		self.mode = mode
		self.output = output
		self.cluster = cluster
		self.jobs = jobs
		self.minimum_barcode_fragments = minimum_barcode_fragments
		self.minimum_cell_fragments = minimum_cell_fragments
		self.extract_mito = extract_mito
		
		self.barcode_tag = barcode_tag
		
		# Figure out operating system just for funzies; not presently used
		self.os = "linux"
		if(platform.platform()[0:5]=="Darwi"):
			self.os = "mac"
		
		if(mode == "bam"):
			self.bamfile = input
			self.bowtie2 = "NA"
			self.bowtie2_index = "bowtie2_index"
			self.bam_name = bam_name
			
			if(bam_name == "default"):
				filename, file_extension = os.path.splitext(self.bamfile)
				self.bam_name = os.path.basename(filename)
		
		#----------------------------------
		# fastq processing specific options -- needs improvement
		#----------------------------------
		if(mode == "fastq"):	
			
			self.bamfile = "NA"
			self.bam_name = "NA"
			
			# Need to align with bowtie2
			self.bowtie2 = get_software_path('bowtie2', bowtie2_path)
			
			# verify bowtie2 index
			bwt2idxfiles = os.popen("ls " + bowtie2_index + "*.bt2*").read().strip().split("\n")
			if(len(bwt2idxfiles) < 6):
				sys.exit("ERROR: cannot find bowtie2 index; specify with --bowtie2-index and make sure to add the prefix along with the folder path")
			else:
				self.bowtie2_index = bowtie2_index
			
			# Collect samples / fastq lists
			# self.samples, self.fastq1, self.fastq2 = inferSampleVectors(input)

				
		#------------------------------------------
		# Verify R and all of its packages are here
		#------------------------------------------
		R = get_software_path('R', r_path)
		check_R_packages(['Rsamtools', 'GenomicAlignments', 'GenomicRanges', 'BiocParallel', 'dplyr', 'SummarizedExperiment'], R)
		self.R = R

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
		
		yield 'script_dir', self.script_dir
		yield 'mode', self.mode
		yield 'output', self.output
		yield 'bamfile', self.bamfile

		yield 'cluster', self.cluster
		yield 'jobs', self.jobs
		yield 'minimum_barcode_fragments', self.minimum_barcode_fragments
		yield 'minimum_cell_fragments', self.minimum_cell_fragments
		
		yield 'extract_mito', self.extract_mito
		yield 'tssFile', self.tssFile
		yield 'blacklistFile', self.blacklistFile
		yield 'bedtoolsGenomeFile', self.bedtoolsGenomeFile
		yield 'R', self.R
		
		yield 'barcode_tag', self.barcode_tag
		yield 'bam_name', self.bam_name
		
		yield 'bowtie2', self.bowtie2
		yield 'bowtie2_index', self.bowtie2_index
		
	
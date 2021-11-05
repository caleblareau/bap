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
from pkg_resources import get_distribution
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

def mitoChr(reference_genome, mito_chromosome):

	if(mito_chromosome != "default"):
		return(mito_chromosome)
	else:
		if(reference_genome in ["hg19", "mm10", "hg38", "mm10"]):
			return("chrM")
		elif(reference_genome in ["GRCh37", "GRCh38", "GRCm37", "GRCm38"]):
			return("MT")
		elif(reference_genome in ["hg38_mm10_s"]):
			return("Homo_sapiens.MT")
		elif(reference_genome in ["hg38.mm10"]):
			return("hg38.chrM")
		elif(reference_genome == "hg19_mm10_c"):
			return("humanM")
		else:
			return("hg19_chrM")

class bap2Project():
	def __init__(self, script_dir, supported_genomes, mode, input, output, name, ncores, reference_genome,
		cluster, jobs, peak_file, barcode_prior,
		minimum_barcode_fragments, barcode_whitelist,
		minimum_jaccard_index, nc_threshold, regularize_threshold, one_to_one, barcoded_tn5, keep_read_names,
		extract_mito, s_temp_files, snakemake_stdout, mapq, max_insert, all_pairs,
		bedtools_genome, blacklist_file, tss_file, mito_chromosome,
		r_path, bedtools_path, samtools_path, bgzip_path, tabix_path, snakemake_path,
		drop_tag, bead_tag, speciesMix):
		
		#----------------------------------
		# Assign straightforward attributes
		#----------------------------------
		self.bap_version = get_distribution('bap-atac').version
		self.script_dir = script_dir
		self.mode = mode
		self.output = output
		self.cluster = cluster
		self.jobs = jobs
		self.mapq = mapq
		self.ncores = ncores

		
		if(minimum_jaccard_index > 1):
			sys.exit("Cannot specify jaccard index > 1; user specified : %s" % str(minimum_jaccard_index))
		
		if(minimum_jaccard_index < 0):
			sys.exit("Cannot specify a nevative jaccard index; user specified : %s" % str(minimum_jaccard_index))
		
		
		# Handle potential whitelist
		if(barcode_whitelist == ""):
			barcode_whitelist = "none"
		else:
			if(os.path.isfile(barcode_whitelist)):
				pass
			else: 
				sys.exit("Could not find the bead whitelist file: %s" % barcode_whitelist)
				
		if(barcode_prior == ""):
			barcode_prior = "none"
		else:
			if(os.path.isfile(barcode_prior)):
				pass
			else: 
				sys.exit("Could not find the barcode prior file: %s" % barcode_prior)
		self.barcode_prior = barcode_prior
		self.barcode_whitelist = barcode_whitelist
		self.nc_threshold = nc_threshold
		self.regularize_threshold = regularize_threshold

		self.minimum_barcode_fragments = minimum_barcode_fragments
		self.minimum_jaccard_index = minimum_jaccard_index
		self.one_to_one = one_to_one
		self.barcoded_tn5 = barcoded_tn5
		self.keep_read_names = keep_read_names
		self.extract_mito = extract_mito
		self.drop_tag = drop_tag
		self.bead_tag = bead_tag
		self.max_insert = max_insert
		
		# Figure out operating system just for funzies; not presently used
		self.os = "linux"
		if(platform.platform()[0:5]=="Darwi"):
			self.os = "mac"
		
		if(mode == "bam"):
			self.bamfile = input
			self.bwa = "NA"
			self.bwa_index = "bwa_index"
			self.name = name
			
			if(name == "default"):
				filename, file_extension = os.path.splitext(self.bamfile)
				self.name = os.path.basename(filename)
		
		if(all_pairs):
			self.proper_pair = "false"
		else:
			self.proper_pair = "true"
		
		#------------------------------------------
		# Verify R and all of its packages are here
		#------------------------------------------
		R = get_software_path('R', r_path)
		#check_R_packages(['Rsamtools', 'GenomicAlignments', 'GenomicRanges', 'dplyr'], R)
		self.R = R
		
		bedtools = get_software_path('bedtools', bedtools_path)
		self.bedtools = bedtools
		samtools = get_software_path('samtools', samtools_path)
		self.samtools = samtools
		bgzip = get_software_path('bgzip', bgzip_path)
		self.bgzip = bgzip
		tabix = get_software_path('tabix', tabix_path)
		self.tabix = tabix
		snakemake = get_software_path('snakemake', snakemake_path)
		self.snakemake = snakemake
		#------------------------
		# Handle reference genome
		#------------------------
		self.reference_genome = reference_genome
		if any(self.reference_genome == s for s in supported_genomes):
			
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
		
		if(reference_genome in ["hg19-mm10", "hg19_mm10_c"]):
			self.speciesMix = "yes"
		else:
			self.speciesMix = "none"
		
		#------------------------------		
		# Make sure all files are valid
		#------------------------------	
		if(bedtools_genome != ""):
			if(os.path.isfile(bedtools_genome)):
				click.echo(gettime() + "Importing custom reference genome files.")
				self.bedtoolsGenomeFile = bedtools_genome
			else: 
				sys.exit("Could not find the bedtools genome file: %s" % bedtools_genome)
		else:
			# Cleaned up case
			click.echo(gettime() + "Found designated reference genome: %s" % self.reference_genome)

				
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
		
		if(peak_file != ""):	
			if(os.path.isfile(peak_file)):
				self.peakFile = peak_file
			else: 
				sys.exit("Could not find the peaks file: %s" % peak_file)
		else:
			self.peakFile = ""
		
		self.mitochr = mitoChr(reference_genome, mito_chromosome)
		
	#--------------------------------------------------------------------------------
	# Define a method to dump the object as a .yaml/dictionary for use in other files
	#--------------------------------------------------------------------------------
	def __iter__(self):
		
		yield 'bap_version', self.bap_version
		yield 'script_dir', self.script_dir
		yield 'mode', self.mode
		yield 'output', self.output
		yield 'bamfile', self.bamfile
		yield 'name', self.name
		yield 'ncores', self.ncores
		
		yield 'cluster', self.cluster
		yield 'jobs', self.jobs
		yield 'peakFile', self.peakFile
		yield 'barcode_prior', self.barcode_prior
		
		yield 'minimum_barcode_fragments', self.minimum_barcode_fragments
		yield 'barcode_whitelist', self.barcode_whitelist

		yield 'minimum_jaccard_index', self.minimum_jaccard_index
		yield 'nc_threshold', self.nc_threshold
		yield 'regularize_threshold', self.regularize_threshold
		yield 'one_to_one', self.one_to_one
		yield 'barcoded_tn5', self.barcoded_tn5
		yield 'keep_read_names', self.keep_read_names
		yield 'proper_pair', self.proper_pair
		
		yield 'extract_mito', self.extract_mito
		yield 'mapq', self.mapq
		
		yield 'tssFile', self.tssFile
		yield 'blacklistFile', self.blacklistFile
		yield 'bedtoolsGenomeFile', self.bedtoolsGenomeFile
		yield 'mitochr', self.mitochr
		
		yield 'R', self.R
		yield 'samtools', self.samtools
		yield 'bedtools', self.bedtools
		yield 'bgzip', self.bgzip
		yield 'tabix', self.tabix
		yield 'snakemake', self.snakemake
		
		yield 'drop_tag', self.drop_tag
		yield 'bead_tag', self.bead_tag
		yield 'max_insert', self.max_insert
		yield 'speciesMix', self.speciesMix
	
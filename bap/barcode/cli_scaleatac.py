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
import glob

import csv
import re
from itertools import groupby
from ..bapHelp import *
from .scaleHelp import *

from pkg_resources import get_distribution
from subprocess import call, check_call

from multiprocessing import Pool, freeze_support
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

@click.command()
@click.version_option()

@click.option('--fastqs', '-f', help='Path of folder created by mkfastq or bcl2fastq; can be comma separated that will be collapsed into one output.')
@click.option('--sample', '-s', help='Prefix of the filenames of FASTQs to select; can be comma separated that will be collapsed into one output.')
@click.option('--output', '-o', default = "Scale_debarcoded", help='A unique run id, used to name output.')

@click.option('--ncores', '-c', default=4, help='Number of cores to be used in parallel de-barcoding. Default = 4')
@click.option('--nreads', '-n', default=10000000, help='Number of reads to be processed in a chunk (both for an output file unit and in parallel processing). Default = 10M')

				
def main(fastqs, sample, output, ncores, nreads):
	
	"""
	scale-barcode: De-barcode samples from scale atac product \n
	Caleb Lareau, clareau <at> stanford <dot> edu \n
	Trims, processes Tn5 barcode, and corrects bead barcode in one shot \n
	"""
	
	__version__ = get_distribution('bap-atac').version
	script_dir = os.path.dirname(os.path.realpath(__file__))
	click.echo(gettime() + "Starting de-barcoding of scale data from bap pipeline v%s" % __version__)
	
	# Define bead barcodes
	barcodesfilepath = script_dir + '/whitelist/737K-cratac-v1.txt.gz'
	with gzip.open(barcodesfilepath, "rt") as my_file:
		barcodesR = my_file.readlines()
	barcodes = [barcode.rstrip() for barcode in barcodesR]
	click.echo(gettime() + "Found and imported " + str(len(barcodes)) + " bead barcodes")	
	global barcodes_set 
	barcodes_set = set(barcodes)

	#----------
	# Deal with the sample parsing
	#----------

	# Function to verify R1/R2/I2 are present for nominated samples
	def verify_sample_from_R1(list_of_R1s):
		verified_R1s = []
		for R1file in list_of_R1s:
			R2file = R1file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
			I2file = R1file.replace("_R1_001.fastq.gz", "_I2_001.fastq.gz")
			if(os.path.exists(R2file) and os.path.exists(I2file)):
				verified_R1s.append(R1file)
		return(verified_R1s)

	# identify all sequencing data that should be parsed for conversion
	def parse_directories(fastqs):
		list_folders = fastqs.split(",")
		list_samples = sample.split(",")

		all_R1s = []
	
		# Look into all supplied folders for specific files:
		for path_to_folder in list_folders:
		
			# Look at all of the possible sample names
			for sample_name_one in list_samples:
				matching_R1s = glob.glob(path_to_folder+"/*" + sample_name_one + "*" + "_R1_001.fastq.gz")
				for file in matching_R1s:
					all_R1s.append(file)
		verified_R1s = verify_sample_from_R1(all_R1s)
		return(verified_R1s)

	# Import files
	R1s_for_analysis = parse_directories(fastqs)
	click.echo(gettime() + "Processing these fastq samples: ")
	for r in R1s_for_analysis:
		click.echo(" - " + r.replace("_R1_001.fastq.gz", ""))
	click.echo(gettime() + "Processing... "+"\n\n")
	
	outfq1file = output + "_R1.fastq.gz"
	outfq2file = output + "_R2.fastq.gz"
	missed_bead = 0
	missed_tn5 = 0
	missed_both = 0
	npass = 0
	
	with gzip.open(outfq1file, "wt") as out_f1:
		with gzip.open(outfq2file, "wt") as out_f2:
			for R1file in R1s_for_analysis:
				R2file = R1file.replace("_R1_001.fastq.gz", "_R2_001.fastq.gz")
				I2file = R1file.replace("_R1_001.fastq.gz", "_I2_001.fastq.gz")
		
				# Read in fastq in chunks the size of the maximum user tolerated number
				it1 = batch_iterator(FastqGeneralIterator(gzip.open(R1file, "rt")), nreads)
				it2 = batch_iterator(FastqGeneralIterator(gzip.open(R2file, "rt")), nreads)
				it3 = batch_iterator(FastqGeneralIterator(gzip.open(I2file, "rt")), nreads)
		
				for i, batch_R1 in enumerate(it1):
					batch_R2 = it2.__next__()
					batch_R3 = it3.__next__()
			
					pool = Pool(processes=ncores)
					pm = pool.map(debarcode_correct_trio, zip(batch_R1, batch_R2, batch_R3))
					pool.close()
				
					# process and write out
					fq_data = list(map(''.join, zip(*[item.pop(0) for item in pm])))
					mismapped_counts = list(map(sum, zip(*[item.pop(0) for item in pm])))
					out_f1.writelines(fq_data[0])
					out_f2.writelines(fq_data[1])
					missed_bead += mismapped_counts[0]
					missed_tn5 += mismapped_counts[1]
					missed_both += mismapped_counts[2]
					npass += mismapped_counts[3]
					
	nfail = missed_tn5 + missed_bead - missed_both
	total = nfail + npass
	print("Total reads that passed: "+str(npass) + "; Total reads that failed: "+str(nfail))
	print("Overall success: "+str(round(npass/(npass+nfail)*100, 2))+"% ")

	print("Of the " +str(nfail) +" reads that failed: "+str(missed_tn5)+" had an incompatible Tn5 barcode")
	print("Of the " +str(nfail) +" reads that failed: "+str(missed_bead)+" had an incompatible bead barcode")
	
if __name__ == "__main__":
    main()

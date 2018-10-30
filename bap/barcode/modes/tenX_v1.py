#!/usr/bin/env python

# Caleb Lareau, Broad Institute
# Implemented: 13 June 2018
# This program will demultiplex BioRad
# scATAC from v2.1 scheme

##### IMPORT MODULES #####
import os
import re
import regex
import sys
import gzip
from barcodeHelp import * # local python script

from optparse import OptionParser
from multiprocessing import Pool, freeze_support
from itertools import repeat

from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator

from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from fuzzysearch import find_near_matches

#### OPTIONS ####
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process raw .fastq reads and make data suitable for downstream processes"

opts.add_option("-a", "--fastq1", help="<Read1> Accepts fastq or fastq.gz")
opts.add_option("-b", "--fastq2", help="<Read2> Accepts fastq or fastq.gz")
opts.add_option("-i", "--fastqI", help="<Index read> Accepts fastq or fastq.gz")
opts.add_option("-f", "--barcodesFile", help="<gzip of the universe of valid barcodes to check")

opts.add_option("-n", "--nreads", default = 5000000, help="Number of reads in each split output file")
opts.add_option("-c", "--ncores", default = 4, help="Number of cores for parallel processing.")

opts.add_option("-x", "--nmismatches", default=1, help="Number of mismatches")
opts.add_option("-o", "--output", help="Output sample convention")

options, arguments = opts.parse_args()

print(options)

# return usage information if no argvs given
if len(sys.argv)==1:
	os.system(sys.argv[0]+" --help")
	sys.exit()

##### INPUTS #####
a = options.fastq1
b = options.fastq2
c = options.fastqI
outname = options.output
o = options.output
bf = options.barcodesFile
cpu = int(options.ncores)
n = int(options.nreads)
n_mismatch = int(options.nmismatches)

# Parse input files
extension = a.split('.')[-1]
if extension == "fastq" or extension == "fq":
	sys.exist("Quitting... GZIP your .fastq files!")
elif extension == "gz":
	print("Found supplied .fastq.gz files")
else:
	sys.exit("ERROR! The input files (-a , -b) a *.fastq.gz")


# Define global variables
dumb = "N"*16

# Define barcodes
with gzip.open(bf, "rt") as my_file:
	barcodesR = my_file.readlines()
barcodes = [barcode.rstrip() for barcode in barcodesR]
print("Found and imported " + str(len(barcodes)) + " barcodes")	

bcPartFile = bf.replace(".txt.gz", ".parts.txt")

part1 = []; part2 = []
with open(bcPartFile, "r") as bcpf:
	line = bcpf.readline()
	while line:
		line = line.strip().split(",")
		bcpart = line[0]
		bcidx = line[1]	
		if(bcidx == "1"):
			part1.append(bcpart)
		if(bcidx == "2"):
			part2.append(bcpart)
		line = bcpf.readline()
		
#------------------------------

def debarcode_10X(trio):
	"""
	Function that is called in parallel
	"""
	# Parse out inputs
	listRead1 = trio[0]; listRead2 = trio[1]; indexRead=trio[2]
	
	# parameters to return
	fq1 = ""
	fq2 = ""
	mm_quant = ""

	npass = 0
	nfail = 0
	
	# Grab attributes
	title1 = listRead1[0]; sequence1 = listRead1[1]; quality1 = listRead1[2]
	title2 = listRead2[0]; sequence2 = listRead2[1]; quality2 = listRead2[2]
	index_sequence = indexRead[1]
	try:
		bc0 = prove_barcode_simple(index_sequence, barcodes)
		if(bc0 == "NA"):
			bc1, mm1 = prove_barcode(index_sequence[0:8], part1, n_mismatch)
			bc2, mm2 = prove_barcode(index_sequence[9:15], part2, n_mismatch)
			bc0 = prove_barcode_simple(bc1+bc2, barcodes)
			# If we still can't find it after slicing it up, fail
			if(bc0 == "NA"):
				return([["", ""], [0, 1], [bc0+",0,0"]]) # pass, fail in the middle
			else:
				# We found a satisfactory, imperfect match
				fq1 = formatRead(bc0 + "_" + title1, sequence1, quality1)
				fq2 = formatRead(bc0 + "_" + title2, sequence2, quality2)
				return([[fq1, fq2], [1, 0], [bc0+","+mm1+","+mm2+","+mm3]]) # pass, fail in the middle
				
		else:
			# Perfect match!
			fq1 = formatRead(bc0 + "_" + title1, sequence1, quality1)
			fq2 = formatRead(bc0 + "_" + title2, sequence2, quality2)
			return([[fq1, fq2], [1, 0], [bc0+",0,0"]]) # pass, fail in the middle

	except:
		return([["", ""], [0, 1], [bc0+",0,0"]]) # pass, fail in the middle


# Define variables to keep track of things that fail
npass = 0
nfail = 0

with gzip.open(a, "rt") as f1:
	with gzip.open(b, "rt") as f2:
		with gzip.open(c, "rt") as f3:
		
			# Establish iterators
			it1 = batch_iterator(FastqGeneralIterator(f1), n)
			it2 = batch_iterator(FastqGeneralIterator(f2), n)
			it3 = batch_iterator(FastqGeneralIterator(f3), n)
		
			# iterate over batches of length n
			for i, batch1 in enumerate(it1):
				batch2 = it2.__next__()
				index = it3.__next__()
				output = o +  "-c" + str(i+1).zfill(3)
			
				# parallel process the barcode processing and accounting of failures.
				pool = Pool(processes=cpu)
				pm = pool.map(debarcode_10X, zip(batch1, batch2, index))
				pool.close()
			
				# Aggregate output
				fqs = list(map(''.join, zip(*[item.pop(0) for item in pm])))
				counts = list(map(sum, zip(*[item.pop(0) for item in pm])))
				mm_values = list(map(''.join, zip(*[item.pop(0) for item in pm])))
			
				# Increment for QC
				npass = npass + counts[0]
				nfail = nfail + counts[1]
			
				# Export one chunk in parallel
				filename1 = output +'_1.fastq.gz'
				filename2 = output +'_2.fastq.gz'
				filenameMM = output +'_mismatches.csv.gz'
			
				pool = Pool(processes=3)
				toke = pool.starmap(chunk_writer_gzip, [(filename1, fqs[0]), (filename2, fqs[1]), (filenameMM, mm_values)])
				pool.close()
			
with open(o + "-parse" + '.sumstats.log', 'w') as logfile:
	# give summary statistics
	logfile.write("\nParsing read pairs:\n" + a + "\n" + b + "\n")
	logfile.write("\n"+str(npass)+" reads parsed with barcodes ("+str(round(npass/(npass+nfail)*100, 2))+"% success)\n")
	logfile.write("Total reads that failed: "+str(nfail)+"\n\n")

	
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

opts.add_option("-n", "--nreads", default = 5000000, help="Number of reads in each split output file")
opts.add_option("-c", "--ncores", default = 4, help="Number of cores for parallel processing.")

opts.add_option("-j", "--constant1", default = "TATGCATGAC", help="Barcode Constant 1")
opts.add_option("-k", "--constant2", default="AGTCACTGAG", help="Barcode Constant 2")
opts.add_option("-l", "--nextera", default="TCGTCGGCAGCGTC", help="Nextera Adaptor Sequence")
opts.add_option("-m", "--me", default="AGATGTGTATAAGAGACAG", help="ME Sequence")

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
outname = options.output
o = options.output

cpu = int(options.ncores)
n = int(options.nreads)

c1 = options.constant1
c2 = options.constant2
nxt = options.nextera
me = options.me
n_mismatch = int(options.nmismatches)

# Infer the length from the adaptors
c1_len = len(c1)
c2_len = len(c2)
nxt_len = len(nxt)
me_len = len(me)

# Parse input files
extension = a.split('.')[-1]
if extension == "fastq" or extension == "fq":
	sys.exist("Quitting... GZIP your .fastq files!")
elif extension == "gz":
	print("Found supplied .fastq.gz files")
else:
	sys.exit("ERROR! The input files (-a , -b) a *.fastq.gz")


# Define global variables
dumb = "N"*7 + "_" + "N"*7 + "_" + "N"*7 
dumb2 = "N"*21

# Define barcodes
barcodes = ["GGACGAC","GCAGTGT","GAGAGGT","GAACCGT","GGTTAGT","GCCTTTG","GATAGAC","GTGGTAG","GTAATAC","CGAGGTC","CATCAGT","CCAAGCT","CCTTAGG","CACGGAC","CAGGCGG","CCGAACC","CACTTCT","CTGGCAT","CGATTAC","TCGTTCT","TGCTACT","TTCCTCT","TACTTTC","TGAATCC","TAGTACC","TTATCAT","TGATTGT","TGGCAAC","TGTTTAG","AGTTTCT","ATGGTGT","ATTGCCT","ACTCAAT","AGACCAT","AGCGAAT","ACCTACC","AGATAGG","AAGGTTC","AGGCATG","GTGGCGC","GGTCGTA","GTGTCCA","GAGGACA","GTCCTTC","GAGCGTG","GATCACC","GTTGATG","CATACGC","CTGCGCC","CGTAGCC","CGCGGCG","CATCTTA","CCAGTCA","CGTTTGA","CCACTTG","CTAACTC","CGAGTGG","TCCTGGC","TGACCGC","TAAGGTA","TCGCGCA","TCATACA","TAAGAGG","TGGAAGG","TCCGCTC","TAACGCC","TGCGTTG","TCGGATG","AGCCGCC","ACACGCG","ACTACGA","AATGGCC","ATGTTCC","ACGTTGG","AGACTTC","ATATAAC","ATAGTTG","GCACAGC","GACAATA","GAATCAA","GCTCCAA","GCGTAGA","GGAAGTT","GGAGCCT","GAATATG","GGTTCAC","CTAGAGC","CGTGATA","CGCCTAA","CGATGCA","CTTGCGA","CCATAAT","CCTATGT","CGCGCTT","CCGCGAT","CGGCCAG","TTGAGGC","TTTCCTA","TCAGCAA","TCCTTAA","TGGACCA","TAGTGTT","TATACTT","TGTCGCT","TACGCAT","TTGTAAG","TGTAGTG","AGTAAGC","ATGAATA","AACGTAA","AATTCCA","AATGATT","AAGTTAT","ACAGCTT","AGCTGAG","ACAGTAC","GGCAGGC","GCGCACG","GAGCTAA","GGTAACA","GCTAATT","GTCGGTT","GGTGTTT","GCGACTC","CTTACCG","CTATTCG","CTAAGAA","CACGCCA","CGGAGGA","CTTGTCC","CTCATTT","CGGATCT","CAGAATT","CGCAATC","TGCGAGC","TTAAGCG","TCTTGTA","TACCGAA","TTCTGCA","TCCAGTT","TGGCCTT","TCGGCGT","TCTGAAC","TCGACAG","AAGCAGC","ATTCACG","AAGTGCG","ATAGGCA","ATTCGTT","ACGTATT","ACCGGCT","AATTGGT","ATTATTC","AACGGTG","GAGTTGC","GGCGGAA","GTTAGGA","GTGCATT","GCCTCGT","GCTTTAT","GTGTGTC","GGCGTCC","CTCTTGC","CGGCTGC","CGGTACG","CGTACAA","CACATGA","CCGGTTT","CGACACT","CCTCCTT","CATGTAT","CTTCATC","CAGAGAG","TATGTGC","TCAAGAC","TTGGTTA","TGGTGAA","TTACAGA","TGAGATT","TTTGGTC","TTGGACT","TTCGTAC","TGAGGAG","ACCATGC","AGAGACC","AGCAACG","ACGAGAA","AACCACA","AACTCTT","ATGAGCT","AGGACGT","AGGATAC"]

#------------------------------

def extract_barcode_v2(sequence1):
	'''
	Function to extract barcodes
	'''

	# Parse out barcodes if we can ID the constants
	try:
		
		# use some approximate, yet generous, indices to facilitate faster matching
		c1_hit = find_near_matches(c1, sequence1[7:25],  max_l_dist=2) 
		c2_hit = find_near_matches(c2, sequence1[23:42], max_l_dist=2)
		nxt_hit = find_near_matches(nxt, sequence1[33:65],  max_l_dist=2)
		me_hit = find_near_matches(me, sequence1[55:], max_l_dist=2)
		
		# Now grab the barcodes
		bc1, mm1 = prove_barcode(sequence1[0:7], barcodes, n_mismatch)
		bc2, mm2 = prove_barcode(sequence1[c1_hit[0][1]+7:c2_hit[0][0]+23], barcodes, n_mismatch)
		bc3, mm3 = prove_barcode(sequence1[c2_hit[0][1]+23:nxt_hit[0][0]+33], barcodes, n_mismatch)
		seq = sequence1[me_hit[0][1]+55:]
		
		return(bc1 + "_" + bc2 + "_" + bc3, seq, str(mm1)+","+str(mm2)+","+str(mm3))
	except:
		return(dumb, sequence1, "0,0,0")
	

def debarcode_v2(duo):
	"""
	Function that is called in parallel
	"""
	# Parse out inputs
	listRead1 = duo[0]; listRead2 = duo[1]
	
	# parameters to return
	fq1 = ""
	fq2 = ""
	mm_quant = ""

	nbc1 = 0
	nbc2 = 0
	nbc3 = 0

	npass = 0
	nfail = 0
	
	# Grab attributes
	title1 = listRead1[0]; sequence1 = listRead1[1]; quality1 = listRead1[2]
	title2 = listRead2[0]; sequence2 = listRead2[1]; quality2 = listRead2[2]
	
	# Return the barcode with underscores + the biological sequence learned	
	barcode, sequence1, mm = extract_barcode_v2(sequence1)
	quality1 = quality1[-1*len(sequence1):]
	
	three = barcode.split("_")
	barcode = "".join(three)

	if("NNNNNNN" in three or len(sequence1) < 10):

		# Something went wrong
		nfail = nfail + 1
		
		if(barcode != dumb2):
			if("NNNNNNN" == three[0]):
				nbc1 += 1
			if("NNNNNNN" == three[1]):
				nbc2 += 1
			if("NNNNNNN" == three[2]):
				nbc3 += 1

	else:
		npass = 1
		fq1 = formatRead("".join(three) + "_" + title1, sequence1, quality1)
		fq2 = formatRead("".join(three) + "_" + title2, sequence2, quality2)
		mm_quant = mm_quant + "".join(three) + "," + mm +"\n"
	return([[fq1, fq2], [nbc1, nbc2, nbc3, npass, nfail], [mm_quant]])


# Define variables to keep track of things that fail
nbc1 = 0
nbc2 = 0
nbc3 = 0

npass = 0
nfail = 0

with gzip.open(a, "rt") as f1:
	with gzip.open(b, "rt") as f2:
		
		# Establish iterators
		it1 = batch_iterator(FastqGeneralIterator(f1), n)
		it2 = batch_iterator(FastqGeneralIterator(f2), n)
		
		# iterate over batches of length n
		for i, batch1 in enumerate(it1):
			batch2 = it2.__next__()
			output = o +  "-c" + str(i+1).zfill(3)
			
			# parallel process the barcode processing and accounting of failures.
			pool = Pool(processes=cpu)
			pm = pool.map(debarcode_v2, zip(batch1, batch2))
			pool.close()
			
			# Aggregate output
			fqs = list(map(''.join, zip(*[item.pop(0) for item in pm])))
			counts = list(map(sum, zip(*[item.pop(0) for item in pm])))
			mm_values = list(map(''.join, zip(*[item.pop(0) for item in pm])))
			
			# Increment for QC
			nbc1 = nbc1 + counts[0]
			nbc2 = nbc2 + counts[1]
			nbc3 = nbc3 + counts[2]
			npass = npass + counts[3]
			nfail = nfail + counts[4]
			
			
			# Export one chunk in parallel
			filename1 = output +'_1.fastq.gz'
			filename2 = output +'_2.fastq.gz'
			#filenameMM = output +'_mismatches.csv.gz'
			
			pool = Pool(processes=2)
			toke = pool.starmap(chunk_writer_gzip, [(filename1, fqs[0]), (filename2, fqs[1])])
			pool.close()
			
with open(o + "-parse" + '.sumstats.log', 'w') as logfile:
	# give summary statistics
	logfile.write("\nParsing read pairs:\n" + a + "\n" + b + "\n")
	logfile.write("\n"+str(npass)+" reads parsed with barcodes ("+str(round(npass/(npass+nfail)*100, 2))+"% success)\n")
	logfile.write("Total reads that failed: "+str(nfail)+"\n\n")
	logfile.write("\nOf reads that could not be parsed that had valid, detectable constants:\n")
	logfile.write(str(nbc1) + " had a bad BC1 barcode sequence\n")
	logfile.write(str(nbc2) + " had a bad BC2 barcode sequence\n")
	logfile.write(str(nbc3) + " had a bad BC3 barcode sequence\n")
	
	
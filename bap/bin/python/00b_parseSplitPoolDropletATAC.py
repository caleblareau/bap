#!/usr/bin/env python

# Caleb Lareau, Broad Institute

# The following program will trim barcoded reads from the R1 of a .fastq file
# 
# With additional (soft/hard)-clipping at the left and right

# Built-in barcodes
bc_v1 = ["AAAGAA","AACAGC","AACGTG","AAGCCA","AAGTAT","AATTGG","ACAAGG","ACCCAA","ACCTTC",
	"ACGGAC","ACTGCA","AGACCC","AGATGT","AGCACG","AGGTTA","AGTAAA","AGTCTG","ATACTT","ATAGCG",
	"ATATAC","ATCCGG","ATGAAG","ATTAGT","CAACCG","CAAGTC","CACCAC","CACTGT","CAGACT","CAGGAG","CATAGA",
	"CCACGC","CCGATG","CCGTAA","CCTCTA","CGAAAG","CGAGCA","CGCATA","CGGCGT","CGGTCC","CGTTAT","CTAGGT",
	"CTATTA","CTCAAT","CTGTGG","CTTACG","CTTGAA","GAAATA","GAAGGG","GACTCG","GAGCTT","GAGGCC","GAGTGA",
	"GATCAA","GCCAGA","GCCGTT","GCGAAT","GCGCGG","GCTCCC","GCTGAG","GCTTGT","GGACGA","GGATTG","GGCCAT",
	"GGGATC","GGTAGG","GGTGCT","GTACAG","GTCCTA","GTCGGC","GTGGTG","GTTAAC","GTTTCA","TAAGCT","TAATAG",
	"TACCGA","TAGAGG","TATTTC","TCAGTG","TCATCA","TCCAAG","TCGCCT","TCGGGA","TCTAGC","TGAATT","TGAGAC",
	"TGCGGT","TGCTAA","TGGCAG","TGTGTA","TGTTCG","TTAAGA","TTCGCA","TTCTTG","TTGCTC","TTGGAT","TTTGGG"
	]
	
spbc_v1 = ["ATT","AGC","CTG","CGA","GTC","GCT","TAG","TGT"]

##### IMPORT MODULES #####
import os
import re
import sys
import gzip
import string
import Levenshtein
import itertools
from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from optparse import OptionParser
from fuzzysearch import find_near_matches

#### OPTIONS ####
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process raw .fastq reads "
opts = OptionParser(usage=usage)
opts.add_option("-a", "--fastq1", help="<Read1> Accepts fastq or fastq.gz")
opts.add_option("-b", "--fastq2", help="<Read2> Accepts fastq or fastq.gz")
opts.add_option("-c", "--barcodes", default = "bc_v1", help="Style of barcodes used; only bc_v1 available currently.")
opts.add_option("-p", "--split-pool-barcodes", default = "spbc_v1", help="Style of barcodes used for split pool; only spbc_v1 available currently.")
opts.add_option("-j", "--constant1", default = "TAGCCATCGCATTGC", help="Barcode Constant 1")
opts.add_option("-k", "--constant2", default="TACCTCTGAGCTGAA", help="Barcode Constant 2")
opts.add_option("-l", "--nextera", default="TCGTCGGCAGCGTC", help="Nextera Adaptor Sequence")
opts.add_option("-m", "--me", default="AGATGTGTATAAGAGACAG", help="ME Sequence")
opts.add_option("-o", "--out", help="Output sample convention")
options, arguments = opts.parse_args()
print(options)

# return usage information if no argvs given
if len(sys.argv)==1:
	os.system(sys.argv[0]+" --help")
	sys.exit()

##### INPUTS #####
f1_in = options.fastq1
f2_in = options.fastq2
outname = options.out

barcodes = options.barcodes
sp_barcodes = options.split_pool_barcodes

c1 = options.constant1
c2 = options.constant2
nxt = options.nextera
me = options.me

# Infer the length from the adaptors
c1_len = len(c1)
c2_len = len(c2)
nxt_len = len(nxt)
me_len = len(me)

#Process barcode input
if(barcodes == "bc_v1"):
	barcodes = bc_v1
if(sp_barcodes == "spbc_v1"):
	sp_barcodes = spbc_v1

barcode1_len = 6
barcode2_len = 6
barcode3_len = 6
barcode4_len = 3

##### OUTPUTS #####
extension = f1_in.split('.')[-1]
if extension == "fastq" or extension == "fq":
	left = open(f1_in)
	right = open(f2_in)
elif extension == "gz":
	left = gzip.open(f1_in, 'rt')
	right = gzip.open(f2_in, 'rt')
else:
	sys.exit("ERROR! The input files (-a , -b) must be a .fq, .fastq, or .fastq.gz")

# Define output file handles
r1_write = gzip.open(outname + '_1.dc.fastq.gz', 'wt')
r2_write = gzip.open(outname + '_2.dc.fastq.gz', 'wt')
fb_write = gzip.open(outname + '_failed_seqs.txt.gz', 'wt')

##### DEFINE FUNCTIONS #####

def print_to_stderr(msg, newline=True):
    '''
    Wrapper to eventually write to stderr
    '''
    sys.stderr.write(str(msg))
    if newline:
        sys.stderr.write('\n')

def proveBarcode(bc):
	'''
	Function that takes a putative barcode and returns the nearest valid one
	'''
	if(bc in barcodes):
		return(bc)
	else:
		eo = process.extractOne(bc, barcodes)
		if(eo[1] >= 83): # 83 comes from 5/6... the score is the score homology
			return(eo[0])
		else:
			return("N"*barcode1_len)

def proveBarcode4(bc):
	'''
	Function that takes a putative barcode from split-pool and returns the nearest valid one
	'''
	if(bc in sp_barcodes):
		return(bc)
	else:
		eo = process.extractOne(bc, sp_barcodes)
		if(eo[1] >= 66): # 66 comes from 2/3... the score is the score homology
			return(eo[0])
		else:
			return("N"*barcode4_len)


def provev11(loi):
	'''
	Function that takes indices of c1, c2, nxt, and me starts
	and returns whether or not it's a valid configuration given the barcode
	length and the attributes of lengths based on version 1 of the barcoding scheme
	'''
	tone = abs((loi[0] + c1_len + barcode2_len + 1) - loi[1]) <= 1
	ttwo =  abs((loi[1] + c2_len + barcode3_len + 1) - loi[2]) <= 1
	tthree = abs((loi[2] + nxt_len + barcode4_len + 1) - loi[3]) <= 1
	tfour = True # need to fix this
	if(tone and ttwo and tthree and tfour):
		return(True, [loi[0] - barcode1_len, loi[1] - barcode2_len, loi[2] - barcode3_len, loi[3] - barcode4_len, loi[3]])
	else:
		return(False, [tone, ttwo, tthree, tfour])

# Need to adjust line in here depending on the version
def resolveIndex(c1_hit, c2_hit, nxt_hit, me_hit):
	'''
	Function that takes hits and then attempts to reconcile the structure for the barcodes.
	If we can resolve the matches, returns true and a list of indices
	If we can't resolve the matches, return a False and a list of match True/False and
	whether or not it was just a digestion error
	'''
	if(len(c1_hit) < 1 or len(c2_hit) < 1 or len(nxt_hit) < 1 or len(me_hit) < 1):
		return(False, [len(c1_hit) > 0, len(c2_hit) > 0, len(nxt_hit) > 0, len(me_hit) > 0, 1])  # the last number is for a hit on the barcode
	else:
		# Define indicies of all that matched that were possible
		ll = [c1_hit, c2_hit, nxt_hit, me_hit]
		indices = [[item[0] for item in match] for match in ll ]
		
		# Make all pairwise combinations and evaluate
		possibilities = list(itertools.product(*indices))
		fnout = list([provev11(x) for x in possibilities])
		verified = [item[0] for item in fnout]
		if any(verified):
			idx = [i for i, x in enumerate(verified) if x]
			idxbest = fnout[idx[0]][1]
			return(True,idxbest)
		else:
			return(False, [len(c1_hit) > 0, len(c2_hit) > 0, len(nxt_hit) > 0, len(me_hit) > 0, 0]) # the last number is for a hit on the barcode

# Define empty variables for counting hits
c1_hitn = 0
c2_hitn = 0
nxt_hitn = 0
me_hitn = 0
barcode_failn = 0
nfail = 0
npass = 0

##### Loop
while 1:

	# process the first file
	seqhead1 = left.readline().rstrip()
	if not seqhead1: break
	seq1 = left.readline().rstrip()
	qualhead1 = left.readline().rstrip()
	qual1 = left.readline().rstrip()

	# process the second file
	seqhead2 = right.readline().rstrip()
	seq2 = right.readline().rstrip()
	qualhead2 = right.readline().rstrip()
	qual2 = right.readline().rstrip()

	# Return all hits up to 2 matches
	c1_hit = find_near_matches(c1, seq1, max_l_dist=2) 
	c2_hit = find_near_matches(c2, seq1, max_l_dist=2)
	nxt_hit = find_near_matches(nxt, seq1, max_l_dist=2)
	me_hit = find_near_matches(me, seq1, max_l_dist=2)
	
	idx = resolveIndex(c1_hit, c2_hit, nxt_hit, me_hit)
	
	# Adjust sequence and write to fastq if successful
	if(idx[0]):
	
		# Determine barcode from coordinate matches
		bc1 = proveBarcode(seq1[max(idx[1][0],0):(idx[1][0]+barcode1_len)])
		bc2 = proveBarcode(seq1[idx[1][1]:(idx[1][1]+barcode2_len)])
		bc3 = proveBarcode(seq1[idx[1][2]:(idx[1][2]+barcode3_len)])
		bc4 = proveBarcode4(seq1[idx[1][3]:(idx[1][3]+barcode4_len)])
		barcode =  bc1 + bc2 + bc3 + "-" + bc4

		# Check to see if any barcodes failed
		if(bc1 != "N"*barcode1_len and bc2 != "N"*barcode2_len and bc3 != "N"*barcode3_len and bc4 != "N"*barcode4_len):
			# Add barcode to sequence header
			seqhead1 = seqhead1.replace(" ", "") + "_CB:" + barcode
			seqhead2 = seqhead2.replace(" ", "") + "_CB:" + barcode
		
			# Cut out the barcode 
			seq1 = seq1[(idx[1][4] + me_len):len(seq1)]
			qual1 = qual1[(idx[1][4] + me_len):len(qual1)]
		
			# Write out to file of sequences that pass
			r1_write.write(seqhead1+"\n");r1_write.write(seq1+"\n")
			r1_write.write(qualhead1+"\n");r1_write.write(qual1+"\n")
			r2_write.write(seqhead2+"\n");r2_write.write(seq2+"\n")
			r2_write.write(qualhead2+"\n");r2_write.write(qual2+"\n")
			npass += 1
		else:
			barcode_failn += 1
			nfail += 1
			fb_write.write(seq1 + " - barcode not found:" + barcode + "\n" )
		
	else:
		# Index attributes of what was successful
		c1_hitn += idx[1][0]*1
		c2_hitn += idx[1][1]*1
		nxt_hitn += idx[1][2]*1
		me_hitn += idx[1][3]*1
		nfail += 1	
		fb_write.write(seq1 + " - sequence structure not found \n")

r1_write.close();r2_write.close(); fb_write.close()
left.close();right.close()

with open(outname+ '.log', 'w') as logfile:
	# give summary statistics
	logfile.write("\nParsing paired end reads:\n" + f1_in + "\n" + f2_in + "\n")
	logfile.write("\n"+str(npass)+" reads parsed with barcodes ("+str(round(npass/(npass+nfail)*100, 2))+"% success)\n")
	logfile.write(str(nfail)+" reads that could not be parsed\n")
	logfile.write("\nOf reads that could not be parsed:\n")
	logfile.write(str(round((1.0-c1_hitn/nfail) * 100)) + "% had a bad C1 sequence\n")
	logfile.write(str(round((1.0-c2_hitn/nfail) * 100)) + "% had a bad C2 sequence\n")
	logfile.write(str(round((1.0-nxt_hitn/nfail) * 100)) + "% had a bad Nextera sequence\n")
	logfile.write(str(round((1.0-me_hitn/nfail) * 100)) + "% had a bad ME sequence\n")
	logfile.write(str(round((barcode_failn/nfail) * 100)) + "% had a bad barcode sequence\n")
	logfile.write("Total reads that failed: "+str(nfail)+"\n\n")



#!/usr/bin/env python

# Original uthor: Jason Buenrostro
# Modified: Caleb Lareau, Broad Institute

# The following program will trim PE reads
# With additional (soft/hard)-clipping at the left and right

##### IMPORT MODULES #####
import os
import re
import sys
import gzip
import string
import Levenshtein
from optparse import OptionParser

#### OPTIONS ####
opts = OptionParser()
usage = "usage: %prog [options] [inputs] This will trim adapters + additional clipping"
opts = OptionParser(usage=usage)
opts.add_option("-a", help="<Read1> Accepts fastq or fastq.gz")
opts.add_option("-b", help="<Read2> Accepts fastq or fastq.gz")
opts.add_option("-o", default = ".", help="Output directory for fastq files")
opts.add_option("-q", default = ".", help="Output directory for log file")
opts.add_option("-s", default = "sample1", help="Sample name")
options, arguments = opts.parse_args()

# return usage information if no argvs given
if len(sys.argv)==1:
	os.system(sys.argv[0]+" --help")
	sys.exit()

##### INPUTS AND OUTPUTS #####
p1_in = options.a
p2_in = options.b
outdir = options.o
logoutdir = options.q
sample = options.s

##### DEFINE FUNCTIONS #####
complement = str.maketrans('ATCGN', 'TAGCN')
def reverse_complement(sequence):
	return sequence.upper().translate(complement)[::-1]

# Align with mismatch, find first and move on, assumes only one
def fuzz_align(s_seq,l_seq,mismatch):
	for i, base in enumerate(l_seq):  # loop through equal size windows
		l_subset = l_seq[i:i+len(s_seq)]
		dist = Levenshtein.distance(l_subset, s_seq)
		if dist <= mismatch:  # find first then break
			return i, dist

#check for file type and open input file
extension = p1_in.split('.')[-1]
if extension == "fastq" or extension == "fq":
	left = open(p1_in)
	right = open(p2_in)
elif extension == "gz":
	left = gzip.open(p1_in, 'rt')
	right = gzip.open(p2_in, 'rt')
else:
	sys.exit("ERROR! The input file2 must be a .fq, .fastq, or .fastq.gz")

##### SCRIPT #####
# initialize variables
i=0;j=0;k=0;tot_b=0
n=20  # match seq
mismatch=1  # only allow 0-1 mismatches for now

r1_write = gzip.open(outdir + "/" + sample + '_1.trim.fastq.gz', 'wt')
r2_write = gzip.open(outdir + "/" + sample + '_2.trim.fastq.gz', 'wt')

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

	# align reads to themselves
	i = i+1  # total reads
	rc_seq2 = reverse_complement(seq2[0:n])
	idx = seq1.rfind(rc_seq2) # look for perfect match
	if idx > 0:
		j = j+1  # 0 mismatchs
	elif mismatch>0:
		hold = fuzz_align(rc_seq2,seq1,mismatch)  # else allow for mismatch
		if hold:
			idx,mis=hold
			if mis == 1:
				k=k+1  # 1 mismatch

	# trim reads if idx exist
	if idx > 0:
		# keep track on how much trimming
		tot_b = tot_b+len(seq2[idx+n:-1]) #track total bases trimmed 
			
		# trim data
		seq1 = seq1[0:idx+n-1] # modified to sub1 because some aligners (bowtie) dont like perfectly overlapping reads
		seq2 = seq2[0:idx+n-1]
		qual1 = qual1[0:idx+n-1]
		qual2 = qual2[0:idx+n-1]
		
	r1_write.write(seqhead1+"\n");r1_write.write(seq1+"\n")
	r1_write.write(qualhead1+"\n");r1_write.write(qual1+"\n")
	r2_write.write(seqhead2+"\n");r2_write.write(seq2+"\n")
	r2_write.write(qualhead2+"\n");r2_write.write(qual2+"\n")

r1_write.close();r2_write.close()
left.close();right.close()

with open(logoutdir + "/" + sample+ '.trim.log', 'w') as logfile:

	# give summary statistics
	logfile.write(str(i)+" sequences\n")
	logfile.write(str(j)+" mismatches0\n")
	logfile.write(str(k)+" mismatches1\n")
	logfile.write(str(round(tot_b/(j+k),2))+" averageTrimLength\n")


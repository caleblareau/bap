#!/usr/bin/env python

import sys
import re
import os
import pysam
from optparse import OptionParser
from contextlib import contextmanager

# Parse out data
opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files and generate the CB (Cell Barcode) ID .bam files"
opts = OptionParser(usage=usage)
opts.add_option("-i", "--input", help="Filename of new .bam file to be generated")
options, arguments = opts.parse_args()

inbam = options.input

bam = pysam.AlignmentFile(inbam, "rb")
basename = re.sub(".bam$", "", os.path.basename(inbam))

# Function to open lots of files
@contextmanager
def multi_file_manager(files, mode='rt'):
    """ Open multiple files and make sure they all get closed. """
    #print(files)
    files = [pysam.AlignmentFile(file, "wb", template = bam) for file in files]
    yield files
    for file in files:
        file.close()

def make_folder(folder):
	"""
	Function to only make a given folder if it does not already exist
	"""
	if not os.path.exists(folder):
		os.makedirs(folder)

# Function to get tag value
def getTag(intags):
    '''
    Checks for aligner-specific read tags and filters
	'''
    for tg in intags:
    	if('XB' == tg[0]):
    		return(tg[1])
    return("NNNNNNNN")

tn5 = ["AAAGAA","AACAGC","AACGTG","AAGCCA","AAGTAT","AATTGG","ACAAGG","ACCCAA","ACCTTC","ACGGAC","ACTGCA","AGACCC","AGATGT","AGCACG","AGGTTA","AGTAAA","AGTCTG","ATACTT","ATAGCG","ATATAC","ATCCGG","ATGAAG","ATTAGT","CAACCG","CAAGTC","CACCAC","CACTGT","CAGACT","CAGGAG","CATAGA","CCACGC","CCGATG","CCGTAA","CCTCTA","CGAAAG","CGAGCA","CGCATA","CGGCGT","CGGTCC","CGTTAT","CTAGGT","CTATTA","CTCAAT","CTGTGG","CTTACG","CTTGAA","GAAATA","GAAGGG","GACTCG","GAGCTT","GAGGCC","GAGTGA","GATCAA","GCCAGA","GCCGTT","GCGAAT","GCGCGG","GCTCCC","GCTGAG","GCTTGT","GGACGA","GGATTG","GGCCAT","GGGATC","GGTAGG","GGTGCT","GTACAG","GTCCTA","GTCGGC","GTGGTG","GTTAAC","GTTTCA","TAAGCT","TAATAG","TACCGA","TAGAGG","TATTTC","TCAGTG","TCATCA","TCCAAG","TCGCCT","TCGGGA","TCTAGC","TGAATT","TGAGAC","TGCGGT","TGCTAA","TGGCAG","TGTGTA","TGTTCG","TTAAGA","TTCGCA","TTCTTG","TTGCTC","TTGGAT","TTTGGG"]

# Make folders with two letters
junk = [make_folder(basename + "_" + i) for i in list(set(bc[0:2] for bc in tn5))]
files = [basename + "_" + bc[0:2] + "/" + basename + "." + bc + ".bam" for bc in tn5]

# Open all the output files and spit out the filtered data
# Based on where the tn5 barcode matches
with multi_file_manager(files) as fopen:
	for read in bam.fetch():
		tag = getTag(read.tags)
		tn5_bc = tag[-6:]
		if(tn5_bc in tn5):
			idx = tn5.index(tn5_bc)
			fopen[idx].write(read)
			

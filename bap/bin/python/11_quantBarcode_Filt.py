#!/usr/bin/env python

import sys
import re
import pysam
from optparse import OptionParser
from collections import Counter
from contextlib import contextmanager

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files filter based on prevalent cell barcodes"
opts = OptionParser(usage=usage)
opts.add_option("--input", "-i", help="Name of the first .bam file")
opts.add_option("--barcodeTag", default = 'CB', help="Name of the first .bam file")
opts.add_option("--min-fragments", "-m", default = 100, help="Minimum number of fragments for consideration")

options, arguments = opts.parse_args()

bamname = options.input
barcodeTag = options.barcodeTag
minFrag = options.min_fragments

bam = pysam.AlignmentFile(bamname, "rb")

def getBarcode(intags):
    '''
    Parse out the barcode per-read
	'''
    for tg in intags:
    	if(barcodeTag == tg[0]):
    		return(tg[1])
    return("NA")
    

def filterReadBarcodes(intags, bc):
    '''
    Checks for aligner-specific read tags and filters
	'''
    for tg in intags:
    	if(barcodeTag == tg[0] and tg[1] in bc):
    		return(tg[1])
    return("NA")

# Loop over bam and extract the barcode of interest
barcodes = ['NA'] 
for read in bam:
	barcodes.append(getBarcode(read.tags))
bam.close()

barcodes = Counter(barcodes)
barcodes = {x : barcodes[x] for x in barcodes if barcodes[x] >= minFrag }

#-------
# Loop back through, filter for positive barcodes, split by chr
#-------

# Function to open lots of files
@contextmanager
def multi_file_manager(files, mode='rt'):
    """ Open multiple files and make sure they all get closed. """
    #print(files)
    files = [pysam.AlignmentFile(file, "wb", template = bam) for file in files]
    yield files
    for file in files:
        file.close()


# Import Barcodes and strip whitespace
bc = list(barcodes.keys())
chrs = ["chr" + x for x in ["1","2","3","4"]]
files = ["out." + chr + ".bam" for chr in chrs]

# Open all the output files and spit out the filtered data
bam = pysam.AlignmentFile(bamname, "rb")
with multi_file_manager(files) as fopen:
	for read in bam:
		if(filterReadBarcodes(read.tags, bc) != "NA" and read.reference_name in chrs):
				idx = chrs.index(read.reference_name)
				fopen[idx].write(read)
bam.close()

# Write out barcode file
bcfile = open("out.barcodequants.csv", "w") 
for k, v in barcodes.items():
    bcfile.write(k +","+ str(v)+"\n")
bcfile.close() 

for bamfilechr in files:
	pysam.index(bamfilechr)




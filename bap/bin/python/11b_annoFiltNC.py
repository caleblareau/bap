#!/usr/bin/env python

import sys
import re
import pysam
from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files and annotate with the number of occurences of a particular read across the experiment."
opts = OptionParser(usage=usage)
opts.add_option("--input", help="Filename of the existing .bam file")
opts.add_option("--output", help="Filename of new .bam file to be generated")
opts.add_option("--dict-file", help="Filepath of the dictionary to convert annotate with NC count")
opts.add_option("--nc-filt", help="Minimum permissible number of barcodes for the read to be detected without filtering it out")
opts.add_option("--bead-barcode", help="bead barcode tag in bam")

options, arguments = opts.parse_args()

bamname = options.input
outname = options.output
dictfile = options.dict_file
ncfilt = options.nc_filt
bb = options.bead_barcode

def get_barcode(intags):
	'''
	Parse out the bead barcode per-read
	'''
	for tg in intags:
		if(bb == tg[0]):
			return(tg[1])
	return("NA")

# Handle a dictionary of read-name : count pairs
d = {}
with open(dictfile) as f:
  for line in f:
    tok = line.split()
    d[tok[0]] = int(tok[1].strip())

bam = pysam.AlignmentFile(bamname, "rb")
out = pysam.AlignmentFile(outname, "wb", template = bam)

# Loop over bam and selectively write to bam de-duplicating droplet barcodes 
# Retain the same style as before 
bp = 0
bp_count = 0
barcode_blacklist_count = dict()

try:
	for read in bam:
		read_name = read.query_name
		read_count = d.get(read_name, 0)
		bead_bc = get_barcode(read.tags)
		
		# Check to verify that the read NC value is less than that which is permissible
		if(read_count <= int(ncfilt)):
			read.tags = read.tags + [("NC", read_count)]
			out.write(read)
		else:
			# Read is blacklisted; add a count to the barcode
			barcode_blacklist_count[bead_bc] = barcode_blacklist_count.get(bead_bc, 0) + 1
			
			
except OSError: # Truncated bam file from previous iteration handle
	print('Finished parsing bam')
	
bam.close()
out.close()

if(len(barcode_blacklist_count) > 0):
	with open(dictfile.replace("_ncRead.tsv", "_ncBarcode.tsv"), 'w') as file:
		for key, value in barcode_blacklist_count.items():
			file.write(key + "\t" + str(value) + "\n")


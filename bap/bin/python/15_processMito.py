#!/usr/bin/env python

import sys
import re
import pysam
from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files and annotate with a drop barcode with extracted mitochondria."
opts = OptionParser(usage=usage)
opts.add_option("--input", help="Filename of the existing .bam file")
opts.add_option("--output", help="Filename of new .bam file to be generated")
opts.add_option("--bead-barcode", help="Read tag of the bead barcode")
opts.add_option("--drop-barcode", help="Read tag of the drop barcode (to be added)")
opts.add_option("--dict-file", help="Filepath of the dictionary to convert bead-barcode to drop-barcode")
opts.add_option("--mitochr", help="Character string for the mitochondria")


options, arguments = opts.parse_args()

bamname = options.input
outname = options.output
bb = options.bead_barcode
db = options.drop_barcode
dictfile = options.dict_file
mitochr = options.mitochr

def getBarcode(intags):
	'''
	Parse out the bead barcode per-read
	'''
	for tg in intags:
		if(bb == tg[0]):
			return(tg[1])
	return("NA")

# Handle a dictionary of bead-barcode : drop-barcode pairs
d = {}
with open(dictfile) as f:
  for line in f:
    tok = line.split()
    d[tok[0]] = tok[1].strip()

bam = pysam.AlignmentFile(bamname, "rb")
out = pysam.AlignmentFile(outname, "wb", template = bam)

# Loop over bam and selectively write to bam de-duplicating droplet barcodes 
barcode_bp = ['NA']
bp = 0

try:
	for read in bam.fetch(mitochr):
		bead_bc = getBarcode(read.tags)
		drop_bc = d.get(bead_bc, "NA")
	
		# Handle droplet barcodes that we want to consider writing out
		if(drop_bc != "NA"):
	
			# New base pair -- no duplicates
			if(read.reference_start != bp):
				bp = read.reference_start
				barcode_bp = [drop_bc]
				read.tags = read.tags + [(db, drop_bc)]
				out.write(read)
				
			else:
				if( not drop_bc in barcode_bp):
					barcode_bp.append(drop_bc)
					read.tags = read.tags + [(db, drop_bc)]
					out.write(read)

except OSError: # Truncated bam file from previous iteration handle
	print('Finished parsing bam')
	
bam.close()
out.close()
pysam.index(outname)


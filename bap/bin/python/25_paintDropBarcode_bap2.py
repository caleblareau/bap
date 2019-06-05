#!/usr/bin/env python

import sys
import re
import pysam
import gzip
from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files and annotate with a drop barcode"
opts = OptionParser(usage=usage)
opts.add_option("--input", help="Filename of the existing .bam file")
opts.add_option("--output", help="Filename of new .bam file to be generated")
opts.add_option("--bead-barcode", help="Read tag of the bead barcode")
opts.add_option("--drop-barcode", help="Read tag of the drop barcode (to be added)")
opts.add_option("--dict-file", help="Filepath of the dictionary to convert bead-barcode to drop-barcode")
opts.add_option("--hq-frags", help="Filepath of the dictionary to file of deduplicated fragments to retain")

options, arguments = opts.parse_args()

bamname = options.input
outname = options.output
bb = options.bead_barcode
db = options.drop_barcode
dictfile = options.dict_file
hqfrags = options.hq_frags

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

# Import a set of read names that are in the HQ deduplicated set
keep_reads = set( (x.split("\t"))[4].strip() for x in open(hqfrags, "rt") )

# Iterate through bam
bam = pysam.AlignmentFile(bamname, "rb")
out = pysam.AlignmentFile(outname, "wb", template = bam)

try:
	for read in bam:
		bead_bc = getBarcode(read.tags)
		drop_bc = d.get(bead_bc, "NA")
	
		# Handle droplet barcodes that we want to consider writing out
		if(drop_bc != "NA"):
			rn = read.query_name
			
			# Look up corresponding read name
			if(rn in keep_reads):
				read.tags = read.tags + [(db, drop_bc)]
				out.write(read)

except OSError: # Truncated bam 
	print('Finished parsing bam')
	
bam.close()
out.close()

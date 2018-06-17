#!/usr/bin/env python

import sys
import re
import pysam
from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files and annotate with a drop barcode"
opts = OptionParser(usage=usage)
opts.add_option("--input", help="Filename of the existing .bam file")
opts.add_option("--output", help="Filename of new .bam file to be generated")
opts.add_option("--bead-barcode", help="Read tag of the bead barcode")
opts.add_option("--drop-barcode", help="Read tag of the drop barcode (to be added)")
opts.add_option("--dict-file", help="Filepath of the dictionary to convert bead-barcode to drop-barcode")

options, arguments = opts.parse_args()

bamname = options.input
outname = options.output
bb = options.bead_barcode
db = options.drop_barcode
dictfile = options.dict_file

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
# Retain the same style as before 
bp = 0
bp_count = 0
barcode_bp_dict = dict()

try:
	for read in bam:
		bead_bc = getBarcode(read.tags)
		drop_bc = d.get(bead_bc, "NA")
	
		# Handle droplet barcodes that we want to consider writing out
		if(drop_bc != "NA"):
			
			# Append the new tag to the read
			read.tags = read.tags + [(db, drop_bc)]
			
			# New base pair -- no duplicates; write out the dictionary and update
			if(read.reference_start != bp):
				
				# Write out old base pair if we have things to write
				if(len(barcode_bp_dict) > 0):
					for key, value in barcode_bp_dict.items():
						#value.tags = value.tags + [("NC", bp_count)]
						out.write(value)
				
				# Now update to the new base pair... in part by wiping the dictionary
				barcode_bp_dict = dict()
				bp = read.reference_start
				barcode_bp_dict[drop_bc] = read
				bp_count = 1
				
				# Else: same base pair -- do one of two things
				# 1) if we've seen the barcode before, then keep only the first sorted value
				# 2) if we haven't seen the barcode before, then append it
				
			else:
				# Still at the same base pair; verify that we haven't seen this barcode before
				bp_count += 1
				if(drop_bc in barcode_bp_dict.keys()):
					old_read = barcode_bp_dict.get(drop_bc)
					old_name = old_read.query_name
					new_name = read.query_name
					
					# Keep newer read if it's a lesser read than the old one
					if(new_name < old_name):
						barcode_bp_dict[drop_bc] = read
					# else keep the old read-- no code needed
					
				# New barcode-- append the read directly
				else:
					barcode_bp_dict[drop_bc] = read

except OSError: # Truncated bam file from previous iteration handle
	print('Finished parsing bam')
	
bam.close()
out.close()

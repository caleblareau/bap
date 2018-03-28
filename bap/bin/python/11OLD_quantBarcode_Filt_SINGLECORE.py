#!/usr/bin/env python

import sys
import re
import pysam
from optparse import OptionParser
from collections import Counter
from contextlib import contextmanager

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to process aligned .bam files to 1) filter based on prevalent cell barcodes;  2) split based on valid chromosomes; and 3) filter for unique fragments"
opts = OptionParser(usage=usage)

opts.add_option("--input", "-i", help="Name of the .bam file to parse")
opts.add_option("--name", "-n", help="Name of the set of .bam files to collate")
opts.add_option("--output", "-o", help="Path to the output directory for these")

opts.add_option("--barcode-tag", default = 'CB', help="Name of the first .bam file")
opts.add_option("--min-fragments", default = 100, help="Minimum number of fragments for barcode consideration")
opts.add_option("--bedtools-genome", help="Filepath to bedtools genome.")
opts.add_option("--min_chromosome_size", default = 1000000, help="Minimum chromosome size (in bp) to be considered for fragment overlap analysis.")

options, arguments = opts.parse_args()

bamname = options.input
name = options.name
out = options.output

barcodeTag = options.barcode_tag
minFrag = int(options.min_fragments)
bedtoolsGenomeFile = options.bedtools_genome
minChromosomeSize = int(options.min_chromosome_size)

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
	tag = getBarcode(intags)
	if(tag in bc):
		return(tag)
	else:
		return("NA")

# Loop over bam and extract the barcode of interest
barcodes = ['NA'] 
barcode_bp = ['NA']
bp = 0

# Loop through; count and de-duplicate for counts
for read in bam:
	read_barcode = getBarcode(read.tags)
	
	# New base pair -- no duplicates
	if(read.reference_start != bp):
		bp = read.reference_start
		barcode_bp = [read_barcode]
		barcodes.append(read_barcode)
	
	# Same base pair -- verify that it's not existing barcodes
	else:
		if( not read_barcode in barcode_bp):
			barcodes.append(read_barcode)
			barcode_bp.append(read_barcode)
			
bam.close()

# Define barcode numbers and filter
barcodes = Counter(barcodes)
barcodes = {x : barcodes[x] for x in barcodes if barcodes[x] >= minFrag }
bc = list(barcodes.keys())

#-------
# Loop back through, filter for positive barcodes, split by chr
#-------

# Function to open lots of files
@contextmanager
def multi_file_manager(files, mode='rt'):
	"""
	Open multiple files and make sure they all get closed.
	"""
	files = [pysam.AlignmentFile(file, "wb", template = bam) for file in files]
	yield files
	for file in files:
		file.close()

# Determine chromosomes for overlaps with fragments based on bedtools genome

chrlens = {}
with open(bedtoolsGenomeFile) as f:
	for line in f:
		tok = line.split("\t")
		chrlens[tok[0]] = tok[1].strip()

chrlenpass = {x : chrlens[x] for x in chrlens if int(chrlens[x]) >= minChromosomeSize }
chrs = list(chrlenpass.keys())

bamchrfiles = [out + "/" + name + "." + chr + ".bam" for chr in chrs]
bamchrrouter = open(out.replace("temp/filt_split", ".internal/samples") + "/" + name + ".chrbam.txt", "w") 
for v in bamchrfiles:
	bamchrrouter.write(v+"\n")
bamchrrouter.close() 

# Open all the output files and spit out the filtered data
barcode_bp = ['NA']
bp = 0
print(chrs)
bam = pysam.AlignmentFile(bamname, "rb")
with multi_file_manager(bamchrfiles) as fopen:
	for read in bam:
		read_barcode = filterReadBarcodes(read.tags, bc)
		try:
			read_chr = read.reference_name
		except ValueError:
			read_chr = "NA"
		if(read_barcode != "NA" and read_chr in chrs):
			idx = chrs.index(read_chr)
			# New base pair -- no duplicates; write out and update
			if(read.reference_start != bp):
				bp = read.reference_start
				barcode_bp = [read_barcode]
				fopen[idx].write(read)
				
			# Same base pair -- verify that it's not existing barcodes
			else:
				# Still at the same base pair; verify that we haven't seen this barcode before
				if( not read_barcode in barcode_bp):
					barcode_bp.append(read_barcode)
					fopen[idx].write(read)
bam.close()

# Write out barcode file
bcfile = open(out.replace("temp/filt_split", "final") + "/" + name + ".barcodequants.csv", "w") 
for k, v in barcodes.items():
	bcfile.write(k +","+ str(v)+"\n")
bcfile.close() 

# Do a final pass where we index errrything
for bamfilechr in bamchrfiles:
	pysam.index(bamfilechr)

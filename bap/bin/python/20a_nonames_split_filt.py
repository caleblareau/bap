import sys
import re
import pysam
import shutil
import os 
import gzip

from multiprocessing import Pool
from optparse import OptionParser
from collections import Counter
from contextlib import contextmanager

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to process aligned .bam files to 1) split by chromosome and 2) report read ID with bead barcode"
opts = OptionParser(usage=usage)

opts.add_option("--input", "-i", help="Name of the .bam file to parse")
opts.add_option("--name", "-n", help="Name of the set of .bam files to collate")
opts.add_option("--output", "-o", help="Path to the output directory for these")

opts.add_option("--mapq", default = 30, help="Minimum mapq for a read to be kept")
opts.add_option("--mito-chr", default = "chrM", help="Designation of mtDNA chromosome")
opts.add_option("--ncores", default = 4, help="Number of cores for parallel processing")
opts.add_option("--bedtools-reference-genome", default = "", help="Reference genome sizes from bedtools")

options, arguments = opts.parse_args()

bamname = options.input
name = options.name
out = options.output

minmapq = float(options.mapq)
mitochr = options.mito_chr

cpu = int(options.ncores)
bedtoolsGenomeFile = options.bedtools_reference_genome
def intersection(lst1, lst2):
	lst3 = [value for value in lst1 if value in lst2]
	return lst3

def listDictToCounter(lst):
	dct = Counter()
	for d in lst:
		for k, v in d.items():
			dct[k] += v 
	return(dct)

def getBarcode(intags):
	'''
	Parse out the barcode per-read
	'''
	for tg in intags:
		if(barcodeTag == tg[0]):
			return(tg[1])
	return("NA")


#---------------------------------------------------------
# Function for writing the bam split by chromosome
#---------------------------------------------------------
def writeSplitBam(two):
	chrom = two[0]
	new_bam = two[1]
	idx = chrs.index(chrom)

	# Split into per-chromosome bam files
	pysam.view(bamname, chrom, '-b' ,'-o', new_bam, catch_stdout=False)
	pysam.index(new_bam)
	return(chrom)

# Handle the chromosomes
chrlens = {}
with open(bedtoolsGenomeFile) as f:
	for line in f:
		tok = line.split("\t")
		chrlens[tok[0]] = tok[1].strip()

chrlenpass = {x : chrlens[x] for x in chrlens }
chrs = list(chrlenpass.keys())
	
# Final loop to write out passing reads
outbam = [out + "/" + name + "." + chr + ".raw.bam" for chr in chrs]

pool = Pool(processes=cpu)
toy_out = pool.map(writeSplitBam, zip(chrs, outbam))
pool.close()

# Make some routing files
bamchrfiles = [out + "/" + name + "." + chr + ".raw" +".bam" for chr in chrs if chr != mitochr]
bamchrrouter = open(out.replace("temp/filt_split", ".internal/samples") + "/" + name + ".chrbam.txt", "w") 
for v in bamchrfiles:
	bamchrrouter.write(v+"\n")
bamchrrouter.close() 

bamchrrouter2 = open(out.replace("temp/filt_split", ".internal/samples") + "/" + name + ".mitochrbam.txt", "w")
bamchrrouter2.write(out + "/" + name + "." + mitochr + ".raw" +".bam"+"\n")
bamchrrouter2.close() 





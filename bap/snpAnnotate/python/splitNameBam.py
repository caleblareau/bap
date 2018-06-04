import sys
import re
import pysam
import shutil
import os 

from multiprocessing import Pool
from optparse import OptionParser
from collections import Counter
from contextlib import contextmanager

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to process aligned .bam files to 1) filter based on prevalent cell barcodes;  2) split based on valid chromosomes; and 3) filter for unique fragments (at the barcode level)"
opts = OptionParser(usage=usage)

opts.add_option("--input", "-i", help="Name of the .bam file to parse")
opts.add_option("--out", "-o", help="Path to the output directory for these files")
opts.add_option("--ncores", default = 4, help="Number of cores for parallel processing.")
opts.add_option("--chrfile", help="Path to chrs")

options, arguments = opts.parse_args()

bamname = options.input
out = options.out
cpu = int(options.ncores)
chrs = [line.strip() for line in open(options.chrfile, 'r')]

# Handle temporary directory structure
temp = out + "/temp"
temp_split = temp + "/01_split"; temp_namesort = temp + "/02_namesort"


def processBam(chr):
	'''
	Function to extract chromosome 
	Here, file is the output per chromosome
	'''	
	temp1 = temp_split + "/splitBam." + chr + ".bam"
	temp2 = temp_namesort + "/nameSort." + chr + ".bam"
	
	split_cmd = "samtools view -b " + bamname + " > " + temp1
	os.system(split_cmd)
	pysam.index(temp1)
	
	namesort_cmd = "samtools sort -n " + temp1 + " -o " + temp2
	os.system(namesort_cmd)


pool = Pool(processes=cpu)
toy_out = pool.map(processBam, chrs)
pool.close()



#!/usr/bin/env python

import sys
import re
import pysam
from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Script to process bam files aligned to human/mouse genome and extract only human genome"
opts = OptionParser(usage=usage)
opts.add_option("--input", help="Filename of the new .bam file to be generated")
opts.add_option("--output", help="Filename of new .bam file to be generated")
options, arguments = opts.parse_args()

bamname = options.input
outname = options.output

bam = pysam.AlignmentFile(bamname, "rb")
hg19template = pysam.AlignmentFile("Exp47_CD8_merged_all.bam", "rb")
out = pysam.AlignmentFile(outname, "wb", template = hg19template)

# Loop over bam and extract the sequence 
for read in bam:
	name = read.reference_name
	ss = name.split("_")
	if(len(ss) > 0):
		if(ss[0] == "hg19"):
			read.query_name = ss[1].strip()
			out.write(read)

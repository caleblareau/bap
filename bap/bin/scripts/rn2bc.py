#!/usr/bin/env python

import sys
import re
import pysam
from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process aligned bam files and generate the XB (Cell Barcode) tag"
opts = OptionParser(usage=usage)
opts.add_option("--input", help="Filename of the new .bam file to be generated")
opts.add_option("--output", help="Filename of new .bam file to be generated")
options, arguments = opts.parse_args()

bamname = options.input
outname = options.output

bam = pysam.AlignmentFile(bamname, "rb")
out = pysam.AlignmentFile(outname, "wb", template = bam)

# Loop over bam and extract the sequence 
for read in bam:
	name = read.query_name
	ss = name.split("_")
	read.query_name = ss[1]
	read.tags = read.tags + [('XB', ss[0])]
	out.write(read)
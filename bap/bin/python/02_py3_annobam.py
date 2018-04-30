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
opts.add_option("--sample", help="ID that will be added to each read in bam file")
opts.add_option("--tag", help="2 digit tag value for bam")
options, arguments = opts.parse_args()

bamname = options.input
outname = options.output
sampleID = options.sample
tag = options.tag

bam = pysam.AlignmentFile(bamname, "rb")
out = pysam.AlignmentFile(outname, "wb", template = bam)

# Loop over bam and extract the sequence 
for read in bam:
	read.tags = read.tags + [(tag, sampleID)]
	out.write(read)
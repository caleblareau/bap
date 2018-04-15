#!/usr/bin/env python

import sys
import re
import pysam
from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to change read tag"
opts = OptionParser(usage=usage)
opts.add_option("--input", help="Filename of the new .bam file to be generated")
opts.add_option("--blacklist", help="Filename of the new .bam file to be generated")
opts.add_option("--tag", help="tag")
opts.add_option("--output", help="Filename of new .bam file to be generated")
options, arguments = opts.parse_args()

bamname = options.input
outname = options.output
blacklistfile = options.blacklist
tag = options.tag

with open(blacklistfile) as file_handle:
    content = file_handle.readlines()
bl = [x.strip() for x in content] 
bl = set(bl)

bam = pysam.AlignmentFile(bamname, "rb")
out = pysam.AlignmentFile(outname, "wb", template = bam)

for read in bam:
    name = read.query_name
    if(name in bl):
        read.set_tag(tag, 4)
    out.write(read)

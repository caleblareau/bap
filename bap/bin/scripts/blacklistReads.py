
#!/usr/bin/env python

import sys
import re
import pysam
from optparse import OptionParser

opts = OptionParser()
usage = "usage: %prog [options] [inputs] Software to process WASP altered alignments"
opts = OptionParser(usage=usage)
opts.add_option("--input", help="Filename of the new .bam file to be generated")
opts.add_option("--output", help="Filename of new .txt file to be generated")
options, arguments = opts.parse_args()

bamname = options.input
outname = options.output

bam = pysam.AlignmentFile(bamname, "rb")
IDs = open(outname + ".bl.txt", 'w') 
log = open(outname + ".bl.log", 'w') 

# Loop over bam and extract whether or not the alignments are in agreement
good = 0
bad = 0

for read in bam:
    name = read.query_name
    ss = name.split("_-_")
    pos_orig = ss[1].strip().split(":")[1]
    pos2 = read.reference_start

    if(abs(float(pos_orig) - pos2) > 100):
        IDs.write(ss[0]+"\n")
        bad = bad + 1
    else:
        good = good +1

log.write("kept: " + str(good) +"\n" + "filtered: " + str(bad) + "\n")
log.close()

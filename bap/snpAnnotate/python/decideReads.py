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
usage = "usage: %prog [options] [inputs] Script to process re-aligned bam and determine read validity"
opts = OptionParser(usage=usage)

opts.add_option("--bam1", "-a", help="Name of the .bam file to parse for haplotype 1")
opts.add_option("--bam2", "-b", help="Name of the .bam file to parse for haplotype 2")
opts.add_option("--stats-file", "-s", help="Name of the stats file to parse")
opts.add_option("--out", "-o", help="Path to the output directory for these files")

options, arguments = opts.parse_args()

bam1 = options.bam1
bam2 = options.bam2
stats_file = options.stats_file
out = options.out

# Infer chromosome name from data file
chrn = str(re.sub("_stats.txt", "", os.path.basename(stats_file)))

def processOneBam(bamfile):
	bam = pysam.AlignmentFile(bamfile, "rb")
	Itr = bam.fetch(until_eof=True)

	nkeep = 0
	whitelist = []

	while(Itr):
		try:
			read1 = Itr.__next__()
			read2 = Itr.__next__()
	
			# Ensure that we are looking at the same reads
			while( not (read1.query_name == read2.query_name)):
				read1 = read2
				read2 = Itr.__next__()
			try:
				# Parse out the old mapping -- read 1
				name1 = read1.query_name
				ss1 = name1.split("_-_")
				vec1 = ss1[0].strip().split(":")
				ref1_old = vec1[0]
				ref1_bp1 = int(vec1[1])
				ref1_bp2 = int(vec1[2])
			
				name2 = read2.query_name
				ss2 = name2.split("_-_")
				vec2 = ss2[0].strip().split(":")
				ref2_old = vec2[0]
				ref2_bp1 = int(vec2[1])
				ref2_bp2 = int(vec2[2])
			
				# Find new positions
				ref1_new = read1.reference_name
				ref2_new = read2.reference_name
				pos1_new = read1.reference_start
				pos2_new = read2.reference_start
				qual1_new = read1.mapping_quality
				qual2_new = read2.mapping_quality
				
				# Given all of this...
				if(qual1_new > 20 and qual2_new > 20):
					if(ref1_new == ref1_old and ref2_new == ref2_old):
						if(min(abs(pos1_new - ref2_bp1), abs(pos1_new - ref1_bp1)) < 100 and
							min(abs(pos2_new - ref2_bp1), abs(pos2_new - ref1_bp1)) < 100 ):
							whitelist.append(ss1[1])
							nkeep += 1
			except ValueError:
				pass
			
		except StopIteration:
			bam.close()	
			break
	return(nkeep, whitelist)

hap1_nkeep, hap1_whitelist = processOneBam(bam1)
hap2_nkeep, hap2_whitelist = processOneBam(bam2)

# Update per-chromosome QC stats
handle_stat = open(stats_file, 'a')
handle_stat.write(chrn + "_Haplotype1kept_" + str(hap1_nkeep) + "\n")
handle_stat.write(chrn + "_Haplotype2kept_" + str(hap2_nkeep) + "\n")
handle_stat.close()

# Write the whitelists
with open(out + "/whitelist_hap1_" + chrn + ".txt", 'w') as t:
	for item in hap1_whitelist:
		t.write("%s\n" % item)
with open(out + "/whitelist_hap2_" + chrn + ".txt", 'w') as t:
	for item in hap2_whitelist:
		t.write("%s\n" % item)

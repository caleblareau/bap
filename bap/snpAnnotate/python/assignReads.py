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
usage = "usage: %prog [options] [inputs] Script to process bam files sorted by name "
opts = OptionParser(usage=usage)

opts.add_option("--input", "-i", help="Name of the .bam file to parse")
opts.add_option("--snps-file", "-s", help="Name of the SNPs file to parse")
opts.add_option("--out", "-o", help="Path to the output directory for these files")


options, arguments = opts.parse_args()

bamname = options.input
out = options.out
snps_file = options.snps_file

# Infer chromosome name from data file
chrn = str(re.sub(".tsv", "", re.sub("SNPs_", "", os.path.basename(snps_file))))

# Build dictionary of SNPs
snps = {}
with open(snps_file, 'r') as t:
	next(t)
	for line in t:
		four = line.strip().split("\t")
		snps[(str(four[0]), int(four[1]), '1')] = four[2]
		snps[(str(four[0]), int(four[1]), '2')] = four[3]


def getBaseAt(read, pos):
    """
    Extract nucleotide within a read at a given set of positions
    """
    nuc = []
    for p in pos:
        try:
            nuc.append(read.seq[p])
        except:
           pass
    return nuc



def getAllelicStatus(gpos, genotype, quality, snps):
    """
    For a given set of genomic position and assoctiated genotype, compare to a snp file and return a code status
    0 : unassigned - no snp information extracted from the read
    1 : genotype from REF genome is found
    2 : genotype from ALT genome is found
    3 : conflicting information extracted from the read

    """

    code = None
    g1_count = 0
    g2_count = 0
    l = len(genotype)

    for i in range(len(genotype)):
        if gpos[i] != None:
            if (str(chrn), int(gpos[i]), '1') in snps.keys() and (str(chrn), int(gpos[i]), '2') in snps.keys():
                if snps[(str(chrn), int(gpos[i]), '1')] == genotype[i] and quality[i] > 10:
                    g1_count+=1
                elif snps[(str(chrn), int(gpos[i]), '2')] == genotype[i] and quality[i] > 10:
                    g2_count+=1


    if g1_count > 0 and g2_count > 0:
        code = 3
    elif g1_count > 0 and g2_count == 0:
        code = 1
    elif g2_count > 0 and g1_count == 0:
        code = 2
    elif g1_count == 0 and g2_count == 0:
        code = 0

    return code

def replace_str_index(text,index=0,replacement=''):
    return '%s%s%s'%(text[:index],replacement,text[index+1:])

def WASP_PE(gpos, genotype, snps, tagval, read, readno, pos_write):
    """
    Do the WASP-style allele swap
    """
    
    name = read.query_name
    name = chrn + ":" + pos_write + "_-_" + name + ' ' + str(readno)

    # Read attributes
    r = read.query_sequence
    q = pysam.array_to_qualitystring(read.query_qualities)

    if(tagval == 1):
        for i in range(len(genotype)):
            if gpos[i] != None:
                if (str(chrn), int(gpos[i]), '1') in snps.keys() and (str(chrn), int(gpos[i]), '2') in snps.keys():
                    r = replace_str_index(r, i, snps[(str(chrn), int(gpos[i]), '2')])
    if(tagval == 2):
        for i in range(len(genotype)):
            if gpos[i] != None:
                if (str(chrn), int(gpos[i]), '1') in snps.keys() and (str(chrn), int(gpos[i]), '2') in snps.keys():
                    r = replace_str_index(r, i, snps[(str(chrn), int(gpos[i]), '1')])

    o = "@"+name+"\n"+r+"\n+\n"+q+"\n"
    return(o)

# Write to permuted fastq files
permuted_hap1_r1 = gzip.open(out + '/' + chrn + '_hap1.realign_1.fastq.gz', 'wb')
permuted_hap1_r2 = gzip.open(out + '/' + chrn + '_hap1.realign_2.fastq.gz', 'wb')

permuted_hap2_r1 = gzip.open(out + '/' + chrn + '_hap2.realign_1.fastq.gz', 'wb')
permuted_hap2_r2 = gzip.open(out + '/' + chrn + '_hap2.realign_2.fastq.gz', 'wb')

bam = pysam.AlignmentFile(bamname, "rb")
Itr = bam.fetch(until_eof=True)

n_hap1 = 0
n_hap2 = 0
n_discordant = 0
n_unassigned = 0

while(Itr):
	try:
		read1 = Itr.__next__()
		read2 = Itr.__next__()
	
		# Ensure that we are looking at the same reads
		while( not (read1.query_name == read2.query_name)):
			read1 = read2
			read2 = Itr.__next__()
	
		# Get genomic positions
		genomePos1 = read1.get_reference_positions(full_length=True)
		genomePos2 = read2.get_reference_positions(full_length=True)
	
		# Infer genotyped bases
		bases1 = getBaseAt(read1,range(0,len(genomePos1)-1))
		bases2 = getBaseAt(read2,range(0,len(genomePos2)-1))
	
		# Tag each read
		tagval1 = getAllelicStatus(genomePos1, bases1, read1.query_qualities, snps)
		tagval2 = getAllelicStatus(genomePos2, bases2, read2.query_qualities, snps)
	
		# Find discordance either between reads or within them
		if(tagval1 == 3 or tagval2 == 3 or (tagval1 + tagval2) == 3):
			n_discordant += 1
	
		# Throwing this in there to catch reads that are behaving oddly
		elif(len(genomePos1) == 0 or len(genomePos2) == 0):
			n_unassigned += 1
	
		elif(tagval1 == 1 or tagval2 == 1):
			pos_write = str(min(genomePos1[read1.query_alignment_start], genomePos2[read2.query_alignment_start])) + ":" + str(max(genomePos1[read1.query_alignment_start], genomePos2[read2.query_alignment_start]))
			n_hap1 += 1                                      # First int is for haplotype; second is for read number
			permuted_hap1_r1.write(str.encode(WASP_PE(genomePos1, bases1, snps, 1, read1, 1, pos_write))) 
			permuted_hap1_r2.write(str.encode(WASP_PE(genomePos2, bases2, snps, 1, read2, 2, pos_write))) 
	
		elif(tagval1 == 2 or tagval2 == 2):
			pos_write = str(min(genomePos1[read1.query_alignment_start], genomePos2[read2.query_alignment_start])) + ":" + str(max(genomePos1[read1.query_alignment_start], genomePos2[read2.query_alignment_start])) 
			n_hap2 += 1                                     # First int is for haplotype; second is for read number
			permuted_hap2_r1.write(str.encode(WASP_PE(genomePos1, bases1, snps, 2, read1, 1, pos_write))) 
			permuted_hap2_r2.write(str.encode(WASP_PE(genomePos2, bases2, snps, 2, read2, 2, pos_write))) 
		else:
			n_unassigned += 1
	except StopIteration:
		bam.close()	
		permuted_hap1_r1.close()
		permuted_hap1_r2.close()
		permuted_hap2_r1.close()
		permuted_hap2_r2.close()
		break

handle_stat = open(out + '/' + chrn + "_stats.txt", 'w')
handle_stat.write(chrn + "_Haplotype1_" + str(n_hap1) + "\n")
handle_stat.write(chrn + "_Haplotype2_" + str(n_hap2) + "\n")
handle_stat.write(chrn + "_Discordant_" + str(n_discordant) + "\n")
handle_stat.write(chrn + "_Unassigned_" + str(n_unassigned) + "\n")
handle_stat.close()


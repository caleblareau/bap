import itertools
import time
import shutil
import re
import os
import sys
import csv
import gzip
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator

# Define bead barcodes
def gettime(): 
	'''
	Matches `date` in unix
	'''
	return(time.strftime("%a ") + time.strftime("%b ") + time.strftime("%d ") + time.strftime("%X ") + 
		time.strftime("%Z ") + time.strftime("%Y")+ ": ")
		
script_dir = os.path.dirname(os.path.realpath(__file__))
barcodesfilepath = script_dir + '/whitelist/737K-cratac-v1.txt.gz'
with gzip.open(barcodesfilepath, "rt") as my_file:
	barcodesR = my_file.readlines()
barcodes = [barcode.rstrip() for barcode in barcodesR]
#print(gettime() + "Found and imported " + str(len(barcodes)) + " bead barcodes")	
global barcodes_set 
barcodes_set = set(barcodes)


# Reformat read for export
def formatRead(title, sequence, quality):
	return("@%s\n%s\n+\n%s\n" % (title, sequence, quality))


def debarcode_correct_trio(trio):
	"""
	Function that is called in parallel
	"""
	# Parse out inputs
	listRead1 = trio[0]; listRead2 = trio[1]; listRead3 = trio[2]

	# parameters to return
	fq1 = ""
	fq2 = ""

	# Grab attributes
	title1 = listRead1[0]; sequence1 = listRead1[1]; quality1 = listRead1[2]
	title2 = listRead2[0]; sequence2 = listRead2[1]; quality2 = listRead2[2]
	title3 = listRead3[0]; sequence3 = listRead3[1]; quality3 = listRead3[2]
	
	# Update sequence with reverse complement
	bio_r3 = str(Seq(sequence3).reverse_complement())
	sequence3 = bio_r3
	
	# update the quality
	quality3 = quality3[::-1]
	
	corrected_bead_barcode = correct_bead_barcode(sequence3)
	
	# Now adjust the scale barcode
	tn5_barcode = sequence2[0:8]
	corrected_tn5_barcode = correct_tn5_barcode(tn5_barcode)
	sequence2 = sequence2[8+19:]
	quality2  = quality2[8+19:]
	
	bad_tn5 = 0
	bad_bead = 0
	bad_both = 0

	if(corrected_tn5_barcode == "N"*8):
		bad_tn5 = 1
	if(corrected_bead_barcode == "N"*16):
		bad_bead = 1
		
	# Return the barcode with underscores + the biological sequence learned
	if(bad_bead + bad_tn5 == 2):
		bad_both = 1
	ok = 1 - bad_tn5 - bad_bead + bad_both
	
	if((ok == 1) & (len(sequence2) > 5)):
		cb = " CB:Z:" + corrected_bead_barcode + corrected_tn5_barcode
		ofq1 = formatRead(title1.split(" ")[0]+cb, sequence1, quality1)
		ofq2 = formatRead(title2.split(" ")[0]+cb, sequence2, quality2)
		return([[ofq1, ofq2], [bad_bead, bad_tn5, bad_both, ok]])
	else:
		return([["",""], [bad_bead, bad_tn5, bad_both, ok]])

#---------------------------
# Functions to deal with 10x barcodes
#---------------------------

#-----
# This code is modified from CellRanger-ATAC but I'm too lazy to factor in base qualities
# https://github.com/10XGenomics/cellranger-atac/blob/main/mro/atac/stages/processing/attach_bcs/__init__.py
#-----
DNA_ALPHABET = 'AGCT'
ALPHABET_MINUS = {char: {c for c in DNA_ALPHABET if c != char} for char in DNA_ALPHABET}
ALPHABET_MINUS['N'] = set(DNA_ALPHABET)
MAXDIST_CORRECT = 2

def gen_nearby_seqs(seq, test_me, maxdist=3):

	allowed_indices = [i for i in range(len(seq)) if seq[i] != 'N']
	required_indices = tuple([i for i in range(len(seq)) if seq[i] == 'N'])
	mindist = len(required_indices)
	if mindist > maxdist:
		return

	for dist in range(mindist + 1, maxdist + 1):
		for modified_indices in itertools.combinations(allowed_indices, dist - mindist):
			indices = set(modified_indices + required_indices)

			for substitutions in itertools.product(
					*[ALPHABET_MINUS[base] if i in indices else base for i, base in enumerate(seq)]):
				new_seq = ''.join(substitutions)
				if new_seq in test_me:
					yield new_seq


def correct_bead_barcode(seq, maxdist=2):
	if seq in barcodes_set:
		return seq
	for test_str in gen_nearby_seqs(seq, barcodes_set, maxdist):
		return(test_str)
	return "N"*16

# Define Tn5 barcodes
tn5_barcodes = [
	"GAACCGCG","AGGTTATA","TCATCCTT",
	"TGGCCGGT","CAATTAAC","ATAATGTG",
	"TCTGTTGG","CTCACCAA","TATTAGCT",
	"ATGTAAGT","GCACGGAC","GGTACCTT",
	"ATCCACTG","GCTTGTCA","CAAGCTAG",
	"TAAGTGGT","CGGACAAC","ATATGGAT",
	"GCTCATTG","ATCTGCCA","CTTGGTAT",
	"GATCTATC","AGCTCGCT","CGGAACTG"
	]
global set_tn5_barcodes
set_tn5_barcodes = set(tn5_barcodes)

def correct_tn5_barcode(tn5seq, maxdist=2):
	if tn5seq in set_tn5_barcodes:
		return tn5seq
	for test_str_tn5 in gen_nearby_seqs(tn5seq, tn5_barcodes, maxdist):
		return(test_str_tn5)
	return "N"*8

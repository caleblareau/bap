import os
import re
import regex
import sys
import gzip
from itertools import repeat

from fuzzywuzzy import fuzz
from fuzzywuzzy import process
from fuzzysearch import find_near_matches

def batch_iterator(iterator, batch_size):
	"""
	Returns lists of tuples of length batch_size.
	"""
	entry = True  # Make sure we loop once
	while entry:
		batch = []
		while len(batch) < batch_size:
			try:
				entry = iterator.__next__()
			except StopIteration:
				entry = None
			if entry is None:
				# End of file
				break
			batch.append(entry)
		if batch:
			yield batch


def chunk_writer_gzip(filename, what):
	'''
	Basic function to write a chunk of a fastq file
	to a gzipped file
	'''
	with gzip.open(filename, 'wt') as out_write:
				out_write.writelines(what)
	return(filename)			

def prove_barcode_simple(bc, valid_set):
	'''
	Function that takes a putative barcode and returns the nearest valid one
	'''
		
	if(bc in valid_set):
		return(bc)
	else:
		return("NA")

def prove_barcode(bc, valid_set, n_mismatch):
	'''
	Function that takes a putative barcode and returns the nearest valid one
	'''
		
	if(bc in valid_set):
		return(bc, "0")
	else:
		eo = process.extractOne(bc, valid_set)
		if(eo[1] >= (len(bc)-n_mismatch)/len(bc)*100): 
			return(eo[0], "1")
		else:
			return("N"*len(bc), "0")
			
def formatRead(title, sequence, quality):
	"""
	Takes three components of fastq file and stiches them together in a string
	"""
	return("@%s\n%s\n+\n%s\n" % (title, sequence, quality))

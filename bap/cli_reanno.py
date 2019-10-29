import click
import os
import os.path
import sys
import shutil
import yaml
import random
import string
import itertools
import time
import pysam

import csv
import re
from itertools import groupby
from .bapHelp import *
from pkg_resources import get_distribution

@click.command()
@click.version_option()

@click.option('--input', '-i', help='Input bam file.')
@click.option('--output', '-o', help='Output bam file.')

@click.option('--sep', '-s', default = "_", help='Separator for reannotation. Assume: {barcode}_{readname} (default delim = "_")')
@click.option('--tag', '-t', default="XB", help='Sam tag for barcode; by default, assume XB.\n\n')

def main(input, output, sep, tag):
	
	"""
	bap-reanno: Reannotate samples that were de-barcoded and aligned \n
	Caleb Lareau, clareau <at> broadinstitute <dot> org \n	
	"""
	
	__version__ = get_distribution('bap-atac').version
	script_dir = os.path.dirname(os.path.realpath(__file__))
	click.echo(gettime() + "Starting re-barcoding from bap pipeline v%s" % __version__)
	
	bam = pysam.AlignmentFile(input, "rb")
	out = pysam.AlignmentFile(output, "wb", template = bam)

	# Loop over bam and extract the sequence 
	for read in bam:
		name = read.query_name
		ss = name.split(sep)
		read.query_name = ss[1]
		read.tags = read.tags + [(tag, ss[0])]
		out.write(read)


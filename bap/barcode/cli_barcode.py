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

import csv
import re
from itertools import groupby
from ..bapHelp import *
from pkg_resources import get_distribution
from subprocess import call, check_call

@click.command()
@click.version_option()

@click.argument('mode', type=click.Choice(['v1.0', 'v2.0', 'v2.1', 'v2.1-multi']))

@click.option('--fastq1', '-a', help='Read 1 of the sequence; should contain the inline barcode in the sequence')
@click.option('--fastq2', '-b', help='Read 2 of the sequence; should contain the inline barcode in the sequence')
@click.option('--output', '-o', default="debarcode", help='Output prefix for processed fastq files. By default, a simple string in the execution directory.\n\n')

@click.option('--ncores', '-c', default=2, help='Number of cores to be used in parallel de-barcoding')
@click.option('--nreads', '-nr', default=5000000, help='Number of reads to be processed in a chunk (both for an output file unit and in parallel processing)')
@click.option('--nmismatches', '-nm', default=1, help='Number of mismatches to be tolerated when doing a section of a barcode matching')


def main(mode, fastq1, fastq2, output, ncores, nreads, nmismatches):
	
	"""
	bap-barcode: De-barcode samples from BioRad bead single cell atac \n
	Caleb Lareau, clareau <at> broadinstitute <dot> org \n
	
	mode = ['v1', 'v2', 'v2-multi'] for bead design\n
	"""
	
	__version__ = get_distribution('bap').version
	script_dir = os.path.dirname(os.path.realpath(__file__))
	click.echo(gettime() + "Starting de-barcoding from bap pipeline v%s" % __version__)
	
	# Parse user settings
	core_call1 = " --fastq1 " + fastq1 + " --fastq2 " + fastq2 + " --ncores " + str(ncores) 
	core_call2 = " --nreads " + str(nreads) + " --nmismatches " + str(nmismatches) + " --output " + output
	core_call = core_call1 + core_call2
	
	# Handle mode to handle the configuration and make the right system call
	if(mode == "v1.0"):
		cmd = 'python '+script_dir+'/modes/biorad_v1.py '
		earlier = " --constant1 " + "TAGCCATCGCATTGC" + " --constant2 " + "TACCTCTGAGCTGAA"
		later = " --nextera " + "TCGTCGGCAGCGTC" + " --me " + "AGATGTGTATAAGAGACAG"
	elif(mode == "v2.0"):
		cmd = 'python '+script_dir+'/modes/biorad_v2.py '
		earlier = " --constant1 " + "TATGCATGAC" + " --constant2 " + "AGTCACTGAG"
		later = " --nextera " + "TGGTAGAGAGGGTG" + " --me " + "AGATGTGTATAAGAGACAG"
	elif(mode == "v2.1"):
		cmd = 'python '+script_dir+'/modes/biorad_v2.py '
		earlier = " --constant1 " + "TATGCATGAC" + " --constant2 " + "AGTCACTGAG"
		later = " --nextera " + "TCGTCGGCAGCGTC" + " --me " + "AGATGTGTATAAGAGACAG"
	else:
		sys.exit(gettime() + "User-supplied mode %s not found!" % mode)
	
	# Assemble the final call
	sys_call = cmd + earlier + later + core_call
	os.system(sys_call)
		

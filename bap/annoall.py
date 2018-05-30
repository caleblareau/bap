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

from pkg_resources import get_distribution
from subprocess import call, check_call
from .bapHelp import *
from .bapProjectClass import *
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs

@click.command()
@click.version_option()

@click.argument('mode', type=click.Choice(['bam',]))

@click.option('--input', '-i', help='Input for bap; varies by which mode is specified; see documentation')
@click.option('--output', '-o', default="bap_out", help='Output directory for analysis; see documentation.')
@click.option('--name', '-n', default="default", help='Name for all of the output files (default: uses the .bam prefix)')

def main(mode, input, output, name,):
	
	"""
	annotateAlleles: Annotate alleles with maternal/paternal haplotype \n
	Caleb Lareau, clareau <at> broadinstitute <dot> org \n
	
	mode = ['bam']\n
	"""
	
	__version__ = get_distribution('bap').version
	script_dir = os.path.dirname(os.path.realpath(__file__))

	click.echo(gettime() + "Starting annotateAlleles from bap pipeline v%s" % __version__)
	
	click.echo(gettime() + "Complete.")
	

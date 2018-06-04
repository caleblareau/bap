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
import glob
import csv
import re
from itertools import groupby

from pkg_resources import get_distribution
from subprocess import call, check_call
from ..bapHelp import *
from .aaProjectClass import *
from ruamel import yaml
from ruamel.yaml.scalarstring import SingleQuotedScalarString as sqs

@click.command()
@click.version_option()

@click.option('--input', '-i', help='Bam file that will be annotated')
@click.option('--snp-table', '-s', help='Input for bap; varies by which mode is specified; see documentation')
@click.option('--output', '-o', default="annotateAlleles_out", help='Output directory for analysis; see documentation.')

@click.option('--name', '-n', default="default", help='Name for all of the output files (default: uses the .bam prefix)')
@click.option('--haplotype-tag', '-ht', default="GH", help='Two letter tag for bam file for phased reads')

@click.option('--ncores', '-c', default=2, help='Number of cores to be used in analysis')
@click.option('--keep-temp-files', '-z', is_flag=True, help='Keep all intermediate files.')

@click.option('--bwa-path', default = "", help='Path to bwa; by default, assumes that bwa is in PATH.')
@click.option('--bwa-index', '-bi', default = "", help='Path to the bwa (.fa) index; should be specified as if you were calling bwa (with file index prefix); should be the same as the original alignment.')


def main(input, snp_table, output,
	name, haplotype_tag,
	ncores, keep_temp_files,
	bwa_path, bwa_index):
	
	"""
	annotateAlleles: Annotate alleles with haplotypes \n
	Caleb Lareau, clareau <at> broadinstitute <dot> org \n
	
	mode = ['bam']\n
	"""
	
	__version__ = get_distribution('bap').version
	script_dir = os.path.dirname(os.path.realpath(__file__))
	print(script_dir)
	click.echo(gettime() + "Starting annotateAlleles from bap pipeline v%s" % __version__)
	
	# Determine chromosomes in bam file
	bam_chrs = []
	for l in pysam.idxstats(input).split('\n'):
		t = l.split("\t")
		if(len(t) > 3):
			if(float(t[2]) > 0):
				bam_chrs.append(t[0])
	
	# Make the project
	p = aaProject(script_dir, input, snp_table, output,
		name, haplotype_tag,
		ncores, keep_temp_files,
		bwa_path, bwa_index)
	
	# Make output folders
	of = output; fin = of + "/final"; temp = of + "/temp"; logs = of + "/logs"
	temp_split = temp + "/01_split"; temp_namesort = temp + "/02_namesort"
	temp_fastq_permuted = temp + "/03_fastq_permuted"; temp_whitelist = temp + "/04_whitelist"
	folders = [of, fin, temp, of + "/.internal", logs, logs + "/bwa", 
			temp_split, temp_namesort, temp_fastq_permuted, temp_whitelist,
			of + "/.internal/parseltongue", of + "/.internal/samples"]
	mkfolderout = [make_folder(x) for x in folders]
	
	# Make internal README files
	if not os.path.exists(of + "/.internal/README"):
		with open(of + "/.internal/README" , 'w') as outfile:
			outfile.write("This folder creates important (small) intermediate; don't modify it.\n\n")
	if not os.path.exists(of + "/.internal/parseltongue/README"):	
		with open(of + "/.internal/parseltongue/README" , 'w') as outfile:
			outfile.write("This folder creates intermediate output to be interpreted by Snakemake; don't modify it.\n\n")

	# Split discriminative SNPs by chr using a simple R call
	split_R = script_dir + "/R/01_splitSNPs.R"
	r_callSplit = " ".join(["Rscript", split_R, temp_split, snp_table])
	os.system(r_callSplit)
	
	# Determine chrs with snps
	snp_chrs = [re.sub(".tsv", "", re.sub(temp_split + "/SNPs_", "", l)) for l in glob.glob(temp_split + "/*.tsv")]
	chrs_go = [value for value in snp_chrs if value in bam_chrs]
	click.echo(gettime() + "Parsing reads from these chromosomes: ")
	click.echo(chrs_go)
	
	# Write them to a file
	with open(of + "/.internal/chrs.txt", 'w') as t:
		for item in chrs_go:
  			t.write("%s\n" % item)
	
	# Split bams by chromosome and sort by name
	filt_split_cmd = 'python ' +script_dir+'/python/splitNameBam.py --input '+p.bamfile + " --ncores " + str(ncores) + ' --chrfile ' + of + "/.internal/chrs.txt" + ' --out ' + of
	os.system(filt_split_cmd)
		
	# Let Snakemake process the chromosome files
	y_s = of + "/.internal/parseltongue/aa.object.yaml"
	with open(y_s, 'w') as yaml_file:
		yaml.dump(dict(p), yaml_file, default_flow_style=False, Dumper=yaml.RoundTripDumper)
			
	snakecmd_chr = 'snakemake --snakefile '+script_dir+'/Snakefile.snpAnnotate.txt --cores '+str(ncores)+' --config cfp="' + y_s + '" -T'
	os.system(snakecmd_chr)
	
	if keep_temp_files:
		click.echo(gettime() + "Temporary files not deleted since --keep-temp-files was specified.")
	else:
		byefolder = of
		shutil.rmtree(byefolder + "/.internal")
		shutil.rmtree(byefolder + "/temp")
		
		click.echo(gettime() + "Intermediate files successfully removed.")
	
	click.echo(gettime() + "Complete.")
	

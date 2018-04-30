#!/usr/bin/env python

from os.path import join
import os
import subprocess
import sys
import shutil
import pysam
import csv
from ruamel import yaml

configFile = sys.argv[1]
fastq1 = sys.argv[2]
fastq2 = sys.argv[3]
sample = sys.argv[4]

with open(configFile, 'r') as stream:
	config = yaml.load(stream, Loader=yaml.Loader)

# General
outdir = config["output"]
script_dir = config["script_dir"]
ncores = "4"
samtools = "samtools"

# Step 2
bwa = config["bwa"]
bwa_index = config["bwa_index"]

# Step 3
barcode_tag = config["barcode_tag"]

# 01 Trim using custom script
trim_py = script_dir + "/bin/python/01_py3_ATACtrim.py"
pycall = "python " + trim_py + " -a "+fastq1+" -b "+fastq2+" -s "+sample+" -o "+outdir+"/01_trimmed -q "+outdir+"/logs/trim"
os.system(pycall)

tfq1 = outdir + "/01_trimmed/" + sample + "_1.trim.fastq.gz"
tfq2 = outdir + "/01_trimmed/" + sample + "_2.trim.fastq.gz"

# 02 Align with bwa
bwalog = outdir+"/logs/bwa/" + sample + ".log"
sortedbam = outdir + '/02_aligned_reads/'+sample+'.st.bam'
bwacall = bwa + ' mem -t '+ncores+' '+ bwa_index +' -1 '+tfq1+' -2 '+tfq2+' 2> '+bwalog+' | ' + samtools + ' view -bS - | ' + samtools + ' sort -@ 4 - -o ' +sortedbam
print("Aligning data with bwa mem")
os.system(bwacall)

# 03 Process .bam files
anno_py = script_dir + "/bin/python/02_py3_annobam.py"
annotate1 = outdir + "/03_processed_reads/" + sample + ".st.anno.bam"

# pysam call
py2call = "python " + anno_py + " --input " +sortedbam+ " --output " +annotate1+ " --sample "+sample+ " --tag " +barcode_tag
os.system(py2call)
pysam.index(annotate1)



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
name = config["name"]
script_dir = config["script_dir"]
mode = config["mode"]

# Ned to update
keepchrs = ["chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22", "chrX"]
mitochrs = ['chrM', 'MT', 'humanM', 'mouseM', 'humanMT', 'mouseMT']
read_quality = "30"
ncores = "4"
# bigwig support

# Step 1
clipl = config["clipl"]
clipr = config["clipr"]
py_trim = str(config["py_trim"])
PEAT = config["PEAT"]
skip_fastqc = str(config["skip_fastqc"])
java = str(config["java"])

# Step 2
bowtie2 = config["bowtie2"]
samtools = config["samtools"]
bowtie2index = config["bowtie2_index"]

# Step 3
max_javamem = config["max_javamem"]
keep_duplicates = str(config["keep_duplicates"])
MarkDuplicatesCall = java + " -Xmx"+max_javamem+"  -jar " + script_dir + "/bin/picard/MarkDuplicates.jar"
CollectInsertCall = java + " -Xmx"+max_javamem+"  -jar " + script_dir + "/bin/picard/CollectInsertSizeMetrics.jar"

# Step 4
macs2 = config["macs2"]
tssFile = config["tssFile"]
macs2_genome_size = config["macs2_genome_size"]



# 01 Trim using custom script
trim_py = script_dir + "/bin/python/py3_ATACtrim.py"
pycall = "python " + trim_py + " -a "+fastq1+" -b "+fastq2+" -l "+clipl+" -r "+clipr+" -s "+sample+" -o "+outdir+"/01_trimmed -t hard -q "+outdir+"/logs/trim"
peatcall = PEAT + " paired -1 "+fastq1+" -2 "+fastq2+" -o "+ outdir+"/01_trimmed/"+sample+" -n "+ncores+" -l 20 -r 0.1 -a 0.1 -g 0.1 --out_gzip"
if not os.path.isfile(outdir+"/logs/trim/"+sample+".trim.log"):
	if(py_trim == "True"):
		print("running py trim")
		os.system(pycall)
	else:
		print("running py trim")
		os.system(pycall)
		#print("Running PEAT trimming")
		#os.system(peatcall)
		#os.system("mv "+outdir+"/01_trimmed/"+sample+"_report.txt "+outdir+"/logs/trim/"+sample+".trim.log")

tfq1 = outdir + "/01_trimmed/" + sample + "_1.trim.fastq.gz"
tfq2 = outdir + "/01_trimmed/" + sample + "_2.trim.fastq.gz"

# 01a fastqc
fastqc_path = shutil.which("fastqc")
if(skip_fastqc == "False" and str(fastqc_path) != "None"):
	fastqc_call = "fastqc -q -j "+ java + " " + tfq1 + " " + tfq2 + " -o "+outdir+"/logs/fastqc" 
	os.system(fastqc_call)
	
# 02 Align with bowtie2
bwt2log = outdir+"/logs/bowtie2/" + sample + ".log"
sortedbam = outdir + '/02_aligned_reads/'+sample+'.all.sorted.bam'
bwt2call = '(' + bowtie2 + ' -X 2000 -p '+ncores+' -x '+ bowtie2index +' --rg-id '+sample+' -1 '+tfq1+' -2 '+tfq2+' | ' + samtools + ' view -bS - | ' + samtools + ' sort -@ 2 - -o ' +sortedbam+') 2> ' + bwt2log 
if not os.path.isfile(sortedbam):
	print("Aligning data with Bowtie2")
	os.system(bwt2call)
	pysam.index(sortedbam)

# 02a mito-- extract if the user says to
mitobam = outdir + "/mito/" + sample + ".mito.bam"
mitocall = samtools + " view -b "+sortedbam+" -o "+mitobam+" " + " ".join(str(i) for i in mitochrs) + " && " + samtools + " index " + mitobam
if not os.path.isfile(mitobam):
	os.system(mitocall)

# 03 Process .bam files
temp1 = outdir + "/03_processed_reads/temp/" + sample + ".temp1.bam"
temp1call = samtools + " view -q " + read_quality + " -f 0x2 "+sortedbam+" -o "+ temp1
if not os.path.isfile(temp1):
	os.system(temp1call)
	pysam.index(temp1)

temp2 = outdir + "/03_processed_reads/temp/" + sample + ".temp2.bam"
temp2call = samtools + " view -b "+temp1+" -o "+temp2+" " + " ".join(str(i) for i in keepchrs)
if not os.path.isfile(temp2):
	os.system(temp2call)
	pysam.index(temp2)

if(mode == "bulk"):
	finalbam = outdir + "/final/bams/"+sample+".proatac.bam"
else:
	finalbam = outdir + "/03_processed_reads/bams/"+sample+".proatac.bam"	
	
rmlog = outdir + "/logs/picard/markdups/"+sample+".MarkDuplicates.log"
run_markdup = MarkDuplicatesCall + " INPUT="+temp2+" OUTPUT="+finalbam+" METRICS_FILE="+rmlog+" REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=SILENT"

if not os.path.isfile(finalbam):
	if(keep_duplicates == "True"):
		os.system("cp " + temp2 + " " + finalbam)
	else:
		os.system(run_markdup)
	pysam.index(finalbam)


##############
# GET QC STATS
##############

# Get Inserts
insertlog = outdir + "/logs/picard/inserts/"+sample+".inserts.log"
histofile = outdir + "/logs/picard/inserts/"+sample+".histogram.pdf"
run_insert = CollectInsertCall + " INPUT="+finalbam+" O="+insertlog+" H="+histofile+"  W=1000 VALIDATION_STRINGENCY=SILENT" 
os.system(run_insert)
Median_Insert=os.popen("grep -A1 'MEDIAN_INSERT_SIZE' "+insertlog+''' | grep -v "MEDIAN_INSERT_SIZE" | awk '{print $1}' ''').read().strip()
Mean_Insert=os.popen("grep -A1 'MEDIAN_INSERT_SIZE' "+insertlog+''' | grep -v "MEDIAN_INSERT_SIZE" | awk '{print $5}' ''').read().strip()
p80_Insert=os.popen("grep -A1 'MEDIAN_INSERT_SIZE' "+insertlog+''' | grep -v "MEDIAN_INSERT_SIZE" | awk '{print $16}' ''').read().strip()

# Do promoter overlap to get TSS rate
ptss = outdir + "/.internal/promoter.tss.bed"
TSSreads = os.popen(samtools + ' view -L ' + outdir + "/.internal/promoter.tss.bed " +finalbam+" | wc -l").read().strip()
ALLreads = os.popen(samtools + ' view '+finalbam+" | wc -l").read().strip()
TSSpercent = str(round(float(TSSreads)/float(ALLreads)*100,2))

# Get other existing summary statistics
Frags = os.popen(samtools + ' flagstat '+sortedbam+''' | head -1 | cut -d " " -f1 | awk '{print $1/2}' ''').read().strip()
Aligned_Reads=os.popen(samtools + ' view -b '+sortedbam+" "+" ".join(str(i) for i in mitochrs + keepchrs) +' | '+samtools+''' flagstat - | head -1 | cut -d " " -f1 | awk '{print $1/2}' ''').read().strip()
Aligned_noMT=os.popen(samtools + ' view -b '+sortedbam+" "+" ".join(str(i) for i in keepchrs) +' | '+samtools+''' flagstat - | head -1 | cut -d " " -f1 | awk '{print $1/2}' ''').read().strip()
MT_frags=os.popen(samtools + ' flagstat '+mitobam+''' | head -5 | tail -n 1 | cut -d " " -f 1 | awk '{print $1/2}' | awk -F. '{print $1}' ''').read().strip()
Dup_Rate=os.popen("grep -A1 'UNPAIRED_READS_EXAMINED' "+rmlog+''' | grep -v "UNPAIRED_READS_EXAMINED" | awk '{print $7*100}' ''').read().strip()
Lib_Size=os.popen("grep -A1 'UNPAIRED_READS_EXAMINED' "+rmlog+''' | grep -v "UNPAIRED_READS_EXAMINED" | awk '{print $8}' ''').read().strip()
Final_frags=os.popen(samtools+' flagstat '+finalbam+''' | head -1 | cut -d " " -f1 | awk '{print $1/2}' ''').read().strip()


# 04 Do macs2 on each sample only if bulk
# Need to build in support for bigwig here
if(mode == "bulk"):
	
	# Do MACS2
	macs2log = outdir + "/logs/macs2/"+sample+".peakcalls.log"
	macs2outdir = outdir + "/04_qc/macs2_each/"
	macs2call = "(" + macs2 + " callpeak -t "+finalbam+" -n " + macs2outdir + sample + " --nomodel --keep-dup all --call-summits -q 0.05 -g " + macs2_genome_size + ") 2> " + macs2log

	if not os.path.isfile(macs2outdir + sample + "_peaks.narrowPeak"):
		os.system(macs2call)
		os.system("mv " + macs2outdir + sample + "_peaks.xls " + outdir + "/logs/macs2/"+sample+"_peaks.xls")
		os.system("mv " + macs2outdir + sample + "_summits.bed " + outdir + "/final/summits/"+sample+"_summits.bed")
	
	# Build TSS Vector
	tssout = outdir + "/logs/tss/" + sample + ".tss.csv"
	Vvec_py = script_dir + "/bin/python/py3_makeVvec.py"
	vv_call = "python "+Vvec_py+" -a "+finalbam+" -b "+tssFile+" -e 2000 -p ends -o " + tssout+" -q "+outdir+"/final/plots/"+sample+".tss.png"
	if not os.path.isfile(tssout):
		os.system(vv_call)
		
# Build final output file
outitems = [ Frags,  TSSpercent,  Final_frags,  Frags,  Dup_Rate,  Lib_Size,  MT_frags,  Aligned_Reads,  Aligned_noMT,  Median_Insert,  Mean_Insert,  p80_Insert]
outnames = ['Frags','TSSpercent','Final_frags','Frags','Dup_Rate','Lib_Size','MT_frags','Aligned_Reads','Aligned_noMT','Median_Insert','Mean_Insert','p80_Insert']
with open(outdir+"/logs/samples/"+sample+".sampleQC.tsv", 'w') as outfile:
	outfile.write("Sample"+"\t"+ "\t".join(outnames)+"\n")
	outfile.write(sample+"\t"+ "\t".join(str(x) for x in outitems)+"\n")
# Stitching together components of bap


*Contact* [Caleb Lareau](mailto:clareau@broadinstitute.org)

## Step 0
```
This workflow assumes a python3 environment with several dependency packages.

Downstream use with bap requires other R packages as well. 

Further, we assume that samtools is in the environment.
```

## Step 1
Debarcode raw reads.

```
bap-barcode --help
```

## Step 2/3
Align and rebarcode with your favoriate aligner. Also, add the barcode to a sam flag instead of the read header. 

Run this for each of the split values, which you can do with a shell loop easily. An example usage for bwa is shown here

```
sh code/02_align_reannotate.sh testrun1-split001 test/output
# calls the 03_bamReannotate.py within the shell script which moves the barcode ID from the fastq header to a bam tag
```

## Step 4
Merge all of the `.bam` files. Index the merged file. 

```
samtools merge XXX
```

## Step 5

Run it through `bap`.

```
bap --help
```
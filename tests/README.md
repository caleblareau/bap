# Testing

In development, these commands are what is assumed to have been run


**bam** mode

### vanilla hg19 with peaks
```
bap bam -i data/test.small.bam -z -bt CB -ji 0.002 -em -pf data/test.small.peaks.bed 
```

### testing species mixing
```
 bap bam -i data/small_mix.bam -bt XB -ji 0.0001 -r hg19-mm10 -z --mapq 0 -bf 100 -cf 100
```


# Testing the pieces

*phasing*

```
annotateAlleles -i allele/mESC_test.bam -z --bwa-index /Volumes/dat/genomes/mm10_bwa/mm10.fa -s allele/CAST_129S_testSNPtable.tsv.gz 

```

```
python /Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/bin/python/markAllelicStatusCL.py -i data/test.small.bam -s data/t.vcf -t HP -r -o ann.bam
bowtie2 -x /Volumes/dat/genomes/hg19_bwt2/hg19 -U ann.realign.fastq.gz 
```

*c1fastq*
```
bap c1fastq --input c1fastqinput.csv -bi /Volumes/dat/genomes/hg19_bwa/hg19.fa -z
```

*step 11*
```
python ../bap/bin/python/11_quantBarcode_Filt.py --input data/small_hg19.bam --bedtools-genome ../bap/anno/bedtools/chrom_hg19.sizes --barcode-tag XB -n h -o data

```

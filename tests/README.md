# Testing

In development, these commands are what is assumed to have been run

**bam** mode

```
bap bam -i data/test.small.bam -z -bt CB -ji 0.002 
```

with mitochondria:
```
bap bam -i data/test.small.bam -z -bt CB -ji 0.002 -em
```

## testing species mixing
```
bap bam -i data/small_mix.bam -bt XB -ji 0.002 -r hg19-mm10 -z
```

## testing the phasing

```
python /Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/bin/python/markAllelicStatusCL.py -i data/test.small.bam -s data/t.vcf -t HP -r -o ann.bam

bowtie2 -x /Volumes/dat/genomes/hg19_bwt2/hg19 -U ann.realign.fastq.gz 

```

# Testing the C1 parsing
```
bap c1fastq --input c1fastqinput.csv -bi /Volumes/dat/genomes/hg19_bwa/hg19.fa -z
```


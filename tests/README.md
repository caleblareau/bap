# Testing

Some basic use cases for testing `bap`'s execution


### Basic comparison

```
time bap bam -i data/jaccardPairsForIGV.bam -bt XB -r hg19 -z -o bap

time bap2 bam -i data/jaccardPairsForIGV.bam -bt XB -r hg19 -z -o bap2
```

**Note:** as discussed elsewhere, there should be no reason to use `bap`; use `bap2` for markedly better performance

### Verify additional output for a species mix experiment or when peaks files are specified

```
time bap2 bam -i data/small_mix.bam -bt XB -ji 0.0001 -r hg19-mm10 -z --mapq 0 -bf 10 -o SM
```

### Test the ability to not merge when a prior is known
```
time bap2 bam -i data/jaccardPairsForIGV.bam -bt XB -r hg19 -z -o bap2 -pf data/test.small.peaks.bed -bp data/jaccardPairsTest_sep.tsv 
```

```
bap-barcode v2.1 -a fastq_br/biorad_v2_R1.fastq.gz -b fastq_br/biorad_v2_R2.fastq.gz -o test
```


#### Docker

Notes for CAL:
```
docker build -t caleblareau/bap bap

docker exec -it caleblareau/bap bash
```


<br><br>

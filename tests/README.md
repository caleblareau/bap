# Testing

Some basic use cases for testing `bap` and `bap2`


## NEW - bap2

```
 time bap2 bam -i data/small_mix.bam -bt XB -ji 0.0001 -r hg19-mm10 -z --mapq 0 -bf 10
```

### Comparison

```
time bap2 bam -i data/jaccardPairsForIGV.bam -bt XB -r hg19 -z -o bap2
time bap bam -i data/jaccardPairsForIGV.bam -bt XB -r hg19 -z -o bap
```

### Dedicated fragment generation

```
bap-frag -i data/jaccardPairsForIGV.bam -be for_frag/jaccardPairsForIGV.barcodeTranslate.tsv -o test_bap_frag -z

# Test snakemake
snakemake --snakefile /Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/bin/snake/Snakefile.bap_frags --cores 12 --config cfp="test_bap_frag/.internal/parseltongue/bap.object.bam.yaml" --stats test_bap_frag/logs/jaccardPairsForIGV.snakemake.stats

```

## debarcoding
```
bap-barcode v1.0 -a fastq_br/biorad_v1_R1.fastq.gz -b fastq_br/biorad_v1_R2.fastq.gz --nmismatches 1
bap-barcode v2.1 -a fastq_br/biorad_v2_R1.fastq.gz -b fastq_br/biorad_v2_R2.fastq.gz --nmismatches 1
bap-barcode v2.1-multi -a fastq_br/biorad_v2-multi_R1.fastq.gz -b fastq_br/biorad_v2-multi_R2.fastq.gz --nmismatches 1


bwa mem /Volumes/dat/genomes/hg19_bwa/hg19.fa debarcode-c001_1.fastq.gz debarcode-c001_2.fastq.gz | samtools view -bS - |  samtools sort -@ 4 - -o debarcode.bam
```

<br><br>

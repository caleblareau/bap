# Test Data

GM12878 from [Jason's original paper](http://www.nature.com/nature/journal/v523/n7561/abs/nature14590.html)

### Sample executions

```
proatac bulk -i data/fastq -o outtest -rg hg19 -bi /Volumes/dat/genomes/hg19_bwt2/hg19 -z

proatac bulk -i data/smallTable.txt -o small -rg hg19 -bi /Volumes/dat/genomes/hg19_bwt2/hg19 -z

proatac counts -i data/bam/RGID -o countsRGIDtest -pf data/test_bed.txt --by-rgid
proatac counts -i data/bam/individualSamples -o countstest -pf data/test_bed.txt

proatac check -i data/fastq -o checktest -rg hg19 -bi /Volumes/dat/genomes/hg19_bwt2/hg19

proatac single -i data/smallTable.txt -o small -rg hg19 -bi /Volumes/dat/genomes/hg19_bwt2/hg19 -z 

proatac summitsToPeaks -i data/summitfiles -o summittest -rg hg19
```


```
python ../proatac/bin/python/py3_makeVvec.py -a small/final/bams/SRR1779979.proatac.bam -b ../proatac/anno/TSS/hg19.refGene.TSS.bed -e 2000 -p ends -o a.txt


```
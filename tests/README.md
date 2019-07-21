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



<br><br>

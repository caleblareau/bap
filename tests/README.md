# Testing

Some basic use cases for testing `bap`'s execution


### Basic comparison

```
time bap bam -i data/jaccardPairsForIGV.bam -bt XB -r hg19 -z -o bap

time bap2 bam -i data/jaccardPairsForIGV.bam -bt XB -r hg19 -z -o bap2
```

**Note:** as discussed elsewhere, there should be no reason to use `bap`; use `bap2` for markedly better performance

### Verify additional output for a species mix experiment

```
time bap2 bam -i data/small_mix.bam -bt XB -ji 0.0001 -r hg19-mm10 -z --mapq 0 -bf 10 -o SM
```


#### Docker

Notes for CAL:
```
docker build -t caleblareau/bap bap

docker exec -it caleblareau/bap bash
```


<br><br>

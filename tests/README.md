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


## testing the phasing

```
 python /Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/bin/python/markAllelicStatus.py -i data/test.small.bam -s data/t.vcf -t HP -r -o ann2.bam
```


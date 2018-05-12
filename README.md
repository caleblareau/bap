# bap
**B**ead-based sc**A**TAC-seq data **P**rocessing

## About
This is an integrated `python` (>= 3.6) and `R` package that
performs analysis and identification of bead doublets and processes
data files for downstream analyses for SureCell single-cell ATAC-seq data. 

## Quality control statistics

[Descriptions about QC metrics](documentation/aboutQC.md)

## Docker

Notes for CAL:
```
docker build -t caleblareau/bap bap

docker exec -it caleblareau/bap bash
```

## Install

**Recommended:**
First, create a `python` virtual environment in some working directory to keep things tidy:

```
python3 -m venv venv3
source venv3/bin/activate
```

Next, install `bap`:

```
git clone https://github.com/caleblareau/bap.git
cd bap
pip3 install -e . 
```

Verify install:

```
bap --help
```


## Simple Test

Try a sample analysis on the test .bam file:

```
cd tests
bap bam -i data/test.small.bam -bt CB
```

Verify that there are files `bap_out/final`

#### Contact
Lead: [Caleb Lareau](mailto:clareau@broadinstitute.org)
<br><br>

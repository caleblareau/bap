# Testing `bap`

This directory holds the automated test suite (`pytest`) and the small data
fixtures it runs against.

## Quick start

Install the package with its test dependencies, then run the fast suite:

```bash
pip install -e .[test]

# Fast: unit + CLI smoke tests only (no external tools needed)
pytest -m "not integration"
```

The fast suite requires no bioinformatics tools and runs in under a second. It
is what CI runs on every push (see `.github/workflows/ci.yml`).

## What each test file covers

| File | Scope | External tools |
|------|-------|----------------|
| `test_helpers.py` | Unit tests for the pure functions in `bap/bapHelp.py` (sequence math, file helpers, tool discovery, the `run_cmd` subprocess wrapper) and `mitoChr`. | none |
| `test_cli_smoke.py` | In-process `click.testing.CliRunner` checks that every active CLI imports, renders `--help`, and reports `--version`; verifies `bap2 support` lists built-in genomes. | none |
| `test_integration.py` | End-to-end `bap2 bam` run against a bundled hg19 `.bam`, asserting the final `.bap.bam`, `.barcodeTranslate.tsv`, and `.fragments.tsv.gz` are produced. | samtools, bedtools, R, snakemake |

Shared fixtures (paths to the bundled data, `small_bam`) live in `conftest.py`.

## Running the integration tests

Integration tests are marked with `@pytest.mark.integration` and **skip
automatically** when `samtools`, `bedtools`, `R`, or `snakemake` are not on
`PATH`, so the default `pytest` run stays green on a bare machine.

To run them, install the external tools (plus the R packages the pipeline uses:
`Rsamtools`, `GenomicAlignments`, `GenomicRanges`, `dplyr`, `data.table`) and
then:

```bash
# Run everything, including integration tests
pytest

# Run only the integration tests
pytest -m integration
```

A convenient way to get the external tools is conda/mamba:

```bash
mamba install -c bioconda -c conda-forge samtools bedtools snakemake bioconductor-genomicalignments bioconductor-rsamtools r-dplyr r-data.table
```

## Test data

- `data/` — small indexed `.bam` files (`jaccardPairsForIGV.bam`,
  `small_mix.bam`, `test.small.bam`), a peaks bed, a whitelist, and a barcode
  prior table used by the pipeline tests.
- `fastq_br/` — small BioRad and Scale-ATAC FASTQs for barcode-parsing tests.
- `for_frag/` — pre-computed annotated fragment files used to exercise the
  fragment / adjacent-Tn5 code paths.

## Manual smoke tests

The commands below are handy for exercising the full pipeline by hand against
the bundled data (run from this `tests/` directory). They require the external
tools listed above.

### Basic run

```bash
bap2 bam -i data/jaccardPairsForIGV.bam -bt XB -r hg19 -z -o bap2
```

> **Note:** the legacy `bap` (v1) command is deprecated — use `bap2`.

### Species-mix / peaks-file output

```bash
bap2 bam -i data/small_mix.bam -bt XB -ji 0.0001 -r hg19-mm10 -z --mapq 0 -bf 10 -o SM
```

### Skip merging when a prior is known

```bash
bap2 bam -i data/jaccardPairsForIGV.bam -bt XB -r hg19 -z -o bap2 \
  -pf data/test.small.peaks.bed -bp data/jaccardPairsTest_sep.tsv
```

### BioRad barcode parsing

```bash
bap-barcode v2.1 -a fastq_br/biorad_v2_R1.fastq.gz -b fastq_br/biorad_v2_R2.fastq.gz -o test
```

### Scale-ATAC barcode parsing

```bash
bap-scale -f fastq_br/scale -s ScaleTest -o ScalePro
```

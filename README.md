# bap
**B**ead-based sc**A**TAC-seq data **P**rocessing

## About
This is an integrated `python` (>= 3.6) and `R` package that
performs analysis and identification of barcode multiplets and processes
data files for downstream analyses for droplet-based single-cell ATAC-seq data. 

## Documentation

[See the bap wiki page](https://github.com/caleblareau/bap/wiki) for details and FAQs about bap.

## Installation

```
pip install bap-atac          # from PyPI
pip install -e .              # from a local checkout (development)
```

Requires Python >= 3.9. The pipeline also shells out to `samtools`, `bedtools`,
`R`, and `snakemake`, which must be available on your `PATH`.

## Testing

Automated tests use `pytest`. Install the test extras and run the fast suite:

```
pip install -e .[test]
pytest -m "not integration"   # unit + CLI smoke tests, no external tools needed
pytest                        # also runs end-to-end tests when the external tools are present
```

See [`tests/README.md`](tests/README.md) for the full testing guide.




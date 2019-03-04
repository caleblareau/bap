
# What's the difference between bap and bap2?
A non-exhaustive list of what the differences are:

- `bap` was slow and relatively unstable system-to-system. `bap2` achieves the same goals
as `bap` but with considerably faster and is more stable both in terms of memory usage and
conceivable error handling. 
- `bap` counted each read (think `samtools idxstats`) when performing QC metrics, which leads
to inflated statistics. `bap2` counts each _fragment_. 
- `bap2` has been engineered to be compliant with CellRangerATAC output. 

# Single end reads?

In the current implementation, **bap** only supports paired-end reads. There's
some support on [Biostars](https://www.biostars.org/p/215988/) that explains
what one should do given single-end sequencing data. 

# I found a bug / error; what do I do? 

Please let us know if you find any errors/inconsistencies in the documentation or code
by filing a new [GitHub Issue](https://github.com/caleblareau/bap/issues).


<br><br>
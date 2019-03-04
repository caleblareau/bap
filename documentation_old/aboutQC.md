
# Current QC report diagnostics
[Caleb Lareau](mailto:clareau@broadinstitute.org)

## Summary summary
This file just computes the mean over all quality controls statistics described
below over barcodes that are defined with `isCell = TRUE`. The idea is that this
should represent the overall spectra of the experiment.

## Individual experiment QC statistics

**DropBarcode**: assigned per droplet (after computational pairing) using the 
previously described formulation. 

**Exp**: shorthand experiment ID for convenience.

**FRIP**: computed as the fraction of fragments overlapping chromatin accessibility
peaks for a given peak annotation (variable dataset to dataset).

**nBarcodesDroplet**: The number of bead barcodes combined to form the specific droplet barcode.

**filteredQ30**: The number of *reads* passing a quality alignment filter of 30.

**isDoublet**: A boolean from nBarcodesDroplet > 1.

**isCell**: defined computationally; presently using `uniqueNuclearFrags > 1000` and `FRIP > 0.5`

**tssPproportion**: proportion of fragments overlapping a TSS +/- 2kb for a given reference set of TSS genes (all Ref-Seq).

**meanInsertSize**: the mean of the insert sizes for each properly paired read per cell.

**medianInsertSize**: same but the median. 

**uniqueNuclearFrags**: number of de-duplicated fragments mapping to the nuclear genome.

**totalNuclearFrags**: number of fragments (including duplicates) mapping to the nuclear genome.

**totalMitoFrags**: number of fragments uniquely aligning to the mitochondrial genome. 

**duplicateProportion**: defined as `(totalNuclearFrags - uniqueNuclearFrags)/ totalNuclearFrags`

**librarySize**: estimated total number of unique reads that can be gathered per cell using the Lander-Waterman
computation, a specific case of the "Coupon Collector" problem in statistics.

<br><br>


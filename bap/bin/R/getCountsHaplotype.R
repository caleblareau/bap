library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(ggplot2)
library(BiocParallel)
library(chromVAR)
library(Matrix)
library(data.table)
library(SummarizedExperiment)

getCountsHaplotype <- function(bamfile, peaks, barcodeTag = "CB", haplotypeTag = "XG",
                               mapqFilter = 10, properPair = TRUE, minFragments = 1000){
  
  # Import Reads
  GA <- readGAlignments(bamfile, param = ScanBamParam(
    flag = scanBamFlag(isMinusStrand = FALSE, isProperPair = properPair),
    tag = c(barcodeTag, haplotypeTag), mapqFilter = mapqFilter))
  GA <- GA[!(as.character(seqnames(GA)) %in% c("chrY", "chrM"))]
  
  # Filter barcodes
  XBtable <- table(mcols(GA)[,barcodeTag])
  whichKeep <- names(XBtable)[as.numeric(XBtable) > minFragments]
  
  GAfilt <- GA[mcols(GA)[,barcodeTag] %in% whichKeep]
  rm(GA)
  grl <- split(GAfilt,mcols(GAfilt)[,barcodeTag], drop = TRUE)
  rm(GAfilt)
  
  # Get reads in peaks as a sparse-matrix-like matrix
  co <- bplapply(grl, function(gr1) {
    ugr <- unique(makeGRangesFromDataFrame(data.frame(
      seqnames = factor(as.character(seqnames(gr1))),
      start = start(gr1), 
      end = end(gr1),
      mcols(gr1)
    ), keep.extra.columns = TRUE))
    co0 <- countOverlaps(peaks, ugr)
    co1 <- countOverlaps(peaks, ugr[mcols(ugr)$XG == 1])
    co2 <- countOverlaps(peaks, ugr[mcols(ugr)$XG == 2])
    return(data.frame(
      peak = c(which(co0!=0), which(co1!=0),which(co2!=0)),
      count = c(co0[co0!=0], co1[co1!=0], co2[co2!=0]),
      genome = c(rep(0, sum(co0!=0)), rep(1, sum(co1!=0)), rep(2, sum(co2!=0)))
    ))
  })
  rm(grl)
  rb <- rbindlist(co, use.names = TRUE, idcol = TRUE)
  rm(co)
  id <- as.factor(rb[[1]])
  
  # Make matrices for counts and genome-specific instances
  m0 <- Matrix::sparseMatrix(i = c(rb[[2]][rb[[4]] == 0], length(peaks)), j = c(as.numeric(id)[rb[[4]] == 0],1), x = c(rb[[3]][rb[[4]] == 0], 0))
  m1 <- Matrix::sparseMatrix(i = c(rb[[2]][rb[[4]] == 1], length(peaks)), j = c(as.numeric(id)[rb[[4]] == 1],1), x = c(rb[[3]][rb[[4]] == 1], 0))
  m2 <- Matrix::sparseMatrix(i = c(rb[[2]][rb[[4]] == 2], length(peaks)), j = c(as.numeric(id)[rb[[4]] == 2],1), x = c(rb[[3]][rb[[4]] == 2], 0))
  colnames(m0) <- levels(id); colnames(m1) <- levels(id); colnames(m2) <- levels(id)
  
  # Make Summarized Experiment
  SE <- SummarizedExperiment(
    rowRanges = peaks,
    assays = list("counts" = m0, "counts_genome1" = m1, "counts_genome2" = m2),
    colData = DataFrame(barcode = colnames(m0))
  )
  SE <- SE[,colSums(m0) > minFragments]
  return(SE)
}

# Import Peaks
peakfile <- "GM12878_consensus.fixedwidthpeaks.bed"
peaks <- sort(sortSeqlevels(suppressWarnings(chromVAR::getPeaks(peakfile))))

bams <- list.files(".", full.names = TRUE, pattern = "*ready2.bam.allele.bam")
#bamfile <- "N711-Exp31-Sample9.ready2.bam.allele.bam"

# Make SE
lo <- lapply(bams, function(bamfile){
  SE <- getCountsHaplotype(bamfile, peaks)
  SE <- filterPeaks(SE)
  
  SEx <- SE[as.character(seqnames(SE)) == "chrX",]
  df <- data.frame(
    full = colSums(assays(SEx)[["counts"]]),
    genome1 = colSums(assays(SEx)[["counts_genome1"]]),
    genome2 = colSums(assays(SEx)[["counts_genome2"]])
  )
})
saveRDS(lo, file = "GM12878_phasedCounts.rds")
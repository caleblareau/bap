library(Rsamtools)
library(GenomicAlignments)
library(GenomicRanges)
library(ggplot2)
library(BiocParallel)
library(chromVAR)
library(Matrix)
library(data.table)
library(SummarizedExperiment)
library(dplyr)

# Function that returns a data.frame of implicated barcode pairs and relevant summary statistics
# Trying a new "jaccard_approx_frag" that is very fast but not techincally a true Jaccard index

findDoubles <- function(bamfile, barcodeTag = "CB", jaccard_frag_cutoff = 0.02, compute_everything = TRUE,
                        mapqFilter = 10, properPair = TRUE, minFragments = 1000){
  
  # Import Reads
  GA <- readGAlignments(bamfile, param = ScanBamParam(
    flag = scanBamFlag(isProperPair = properPair),
    tag = c(barcodeTag), mapqFilter = mapqFilter))
  GA <- GA[!(as.character(seqnames(GA)) %in% c("chrY", "chrM"))]
  
  # Filter barcodes
  XBtable <- table(mcols(GA)[,barcodeTag])
  whichKeep <- names(XBtable)[as.numeric(XBtable) > minFragments]
  nKept <- as.numeric(XBtable)[as.numeric(XBtable) > minFragments]
  names(nKept) <- whichKeep
  GAfilt <- GA[mcols(GA)[,barcodeTag] %in% whichKeep]
  rm(GA); rm(XBtable)
  
  # Make a couple of handy smaller objects
  barcode <- mcols(GAfilt)[,barcodeTag]
  
  # Find exact fragments
  ov <- findOverlaps(GRanges(GAfilt), GRanges(GAfilt), type = "equal")
  rm(gr)
  
  # Make a dataframe of all combinations that have fragments overlapping
  hugeDF <-
    data.frame(
      bc1 = barcode[queryHits(ov)[ queryHits(ov) !=  subjectHits(ov)]],
      bc2 = barcode[subjectHits(ov)[ queryHits(ov) !=  subjectHits(ov)]],
      stringsAsFactors = FALSE
    )
  rm(ov)
  
  # One dplyr command to save the day -- does everything 
  # Computes an "approximate jaccard frag" as it can be > 1 in certain weird duplicate settings
  hugeDF %>% filter(bc1 != bc2) %>%  # filter out reads that map to the same barcode
    mutate(barc1 = ifelse(bc1 > bc2, bc1, bc2),
           barc2 = ifelse(bc1 > bc2, bc2, bc1)) %>%  # switch to respect alphabetical order if necessary
    select(-one_of(c("bc1", "bc2"))) %>% # drop non-ordered barcode columns
    group_by(barc1, barc2) %>% summarise(n_both = n()) %>%  # group and count
    mutate(n_barc1 = nKept[barc1], n_barc2 = nKept[barc2]) %>% # add the baseline number of reads
    mutate(jaccard_approx_frag = (n_both)/(n_barc1 + n_barc2- n_both)) %>% # compute the Jaccard index and output
    filter(jaccard_approx_frag > jaccard_frag_cutoff) %>% arrange(desc(jaccard_approx_frag)) %>% as.data.frame() -> implicatedPairs # filter based on 
  rm(hugeDF)
  
  # Only compute the exact jaccard_frag and
  # jaccard_bp metrics if user specicies to do so
  if(compute_everything){
    
    # still slow... function to compute jaccard_bp
    jaccard_bp_frag_exact <- function(gr_a, gr_b) {
      
      # per base-pair a la Samtools
      intersection <- sum(width(intersect(gr_a, gr_b)))
      union <- sum(width(union(gr_a, gr_b)))
      j_bp <- intersection/union
      
      ov <- findOverlaps(gr_a, gr_b, type = "equal")
      j_frag <- length(ov)/(length(gr_a) + length(gr_b) - length(ov))
      
      return(c(j_bp, j_frag))
    }
    
    # loop back through pairs and compute jaccard_bp
    GAsmall <- GAfilt[barcode %in% unique(c(as.character(implicatedPairs[,"barc1"]),
                                            as.character(implicatedPairs[,"barc2"])))]
    
    # Compute the exact metrics
    twoOut <-  sapply(1:dim(implicatedPairs)[1], function(i){
      out2 <- jaccard_bp_frag_exact(
        GRanges(GAsmall[mcols(GAsmall)[,barcodeTag] == as.character(implicatedPairs[i,"barc1"])]),
        GRanges(GAsmall[mcols(GAsmall)[,barcodeTag] == as.character(implicatedPairs[i,"barc2"])])
      )
      out2
    })
    implicatedPairs$jaccard_bp <- twoOut[1,]
    implicatedPairs$jaccard_frag<- twoOut[2,]
    rm(GAsmall)
  }
  
  rm(barcode)
  return(implicatedPairs)
}

# Execute
bamfile <- "N711-Exp31-Sample9.ready2.bam.allele.bam"
implicatedPairs <- findDoubles(bamfile)


p1 <- ggplot(implicatedPairs, aes(x = jaccard_frag, y = jaccard_bp)) + geom_point() +
  theme_bw() + labs(x = "Jaccard Frag", y = "Jaccard BP")

p2 <- ggplot(implicatedPairs, aes(x = jaccard_frag, y = jaccard_approx_frag)) + geom_point() +
  theme_bw() + labs(x = "Jaccard Frag", y = "Jaccard Approx. Frag")

cowplot::ggsave(cowplot::plot_grid(p1, p2, nrow = 1), filename = "compare2jaccard.pdf", 
                width = 8, height = 4)

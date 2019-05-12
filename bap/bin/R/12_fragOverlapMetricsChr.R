options(warn=-1)

suppressMessages(suppressWarnings(library(Rsamtools)))
suppressMessages(suppressWarnings(library(GenomicAlignments)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))

"%ni%" <- Negate("%in%")

# This function / script for a given .bam file saves a .rds file
# For each pair of barcodes given how many times they share an identical fragment in this 
# given chromosome

# TO DO:
# parameterize mapq and properpair

args <- commandArgs(trailingOnly = TRUE)

if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}
bamfile <- args[i+1]
barcodeTag <- args[i+2]
barcodeQuantFile <- args[i+3]
blacklistRegionsFile <- args[i+4]

# Don't execute-- strictly for testing
if(FALSE){
  base <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests"
  bamfile <- paste0(base, "/", "bap_out/temp/filt_split/test.small.chr21.bam")
  barcodeTag <- "CB"
  barcodeQuantFile <- paste0(base, "/", "bap_out/final/test.small.barcodequants.csv")
  blacklistRegionsFile <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/anno/blacklist/hg19.full.blacklist.bed"
}


# Import blacklist file
bldf <- data.frame(fread(blacklistRegionsFile))
blgr <- makeGRangesFromDataFrame(bldf, seqnames.field = "V1", start.field = "V2", end.field = "V3")

# Establish the barcode counts normally
bq <- data.frame(fread(barcodeQuantFile, header = TRUE, sep = ","))
keepBarcodesGlobal <- as.character(bq[,1])
rm(bq)

findDoubles_df <- function(bamfile, barcodeTag, mapqFilter = 0, properPair = TRUE){
  
  rdsOut <- gsub(".bam", "_overlapCount.rds", gsub("/filt_split/", "/frag_overlap/", bamfile))
  
  # Import Reads
  GA <- readGAlignments(bamfile, param = ScanBamParam(
    flag = scanBamFlag(isProperPair = properPair),
    tag = c(barcodeTag), mapqFilter = mapqFilter))
  
  # Filter 1 for eligible barcodes
  GAfilt1 <- GA[mcols(GA)[,barcodeTag] %in% keepBarcodesGlobal]
  rm(GA)
  
  # Filter 2 for blacklist
  ov_bl <- findOverlaps(GAfilt1, blgr)
  GAfilt <- GAfilt1[1:length(GAfilt1) %ni% queryHits(ov_bl)]
  rm(GAfilt1)
  
  # Pull out barcode for retained
  barcode <- mcols(GAfilt)[,barcodeTag]
  
  # Find exact fragments
  ov <- findOverlaps(GRanges(GAfilt), GRanges(GAfilt), type = "equal")
  
  # Determine baseline numbers
  XBtable <- table(mcols(GAfilt)[,barcodeTag])
  whichKeep <- names(XBtable)
  nKept <- as.numeric(XBtable); names(nKept) <- whichKeep
  rm(XBtable)
  
  # Make a dataframe of all combinations that have fragments overlapping
  hugeDF <-
    data.frame(
      bc1 = barcode[queryHits(ov)[ queryHits(ov) !=  subjectHits(ov)]],
      bc2 = barcode[subjectHits(ov)[ queryHits(ov) !=  subjectHits(ov)]],
      stringsAsFactors = FALSE
    )
  rm(ov)
  
  # One dplyr command to save the day -- determine occurences of overlapping reads
  hugeDF %>% filter(bc1 != bc2) %>%  # filter out reads that map to the same barcode
    mutate(barc1 = ifelse(bc1 > bc2, bc1, bc2),
           barc2 = ifelse(bc1 > bc2, bc2, bc1)) %>%  # switch to respect alphabetical order if necessary
    select(-one_of(c("bc1", "bc2"))) %>% # drop non-ordered barcode columns
    group_by(barc1, barc2) %>% summarise(n_both = n()) %>%  # group and count
    mutate(n_barc1 = nKept[barc1], n_barc2 = nKept[barc2]) %>% as.data.frame() -> implicatedPairs # keep baseline numbers
  rm(hugeDF)
  rm(barcode)
  saveRDS(implicatedPairs, file = rdsOut)
}

findDoubles_df(bamfile, barcodeTag)

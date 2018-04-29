suppressMessages(suppressWarnings(library(Rsamtools)))
suppressMessages(suppressWarnings(library(GenomicAlignments)))
suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(tools)))

"%ni%" <- Negate("%in%")

options(scipen=999)

# Working with the raw .bam file to generate QC statistics 
# Per barcode 

#-----------------
# Helper Functions
#-----------------

# See https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/DuplicationMetrics.java#L115
estimateLibrarySize <- function(nTotal, nUnique){
  
  f <- function(x, c, n) {
    return(c / x - 1 + exp(-n / x))
  }
  
  m = 1
  M = 100
  
  nDuplicates <- nTotal - nUnique
  
  # Checks
  stopifnot(nTotal > 0)
  stopifnot(nDuplicates >0)
  stopifnot(nUnique > 0)
  
  if (nUnique >= nTotal | (f(m * nUnique, nUnique, nTotal) < 0)) {
    stop("Error: invalid inputs")
  }
  
  while (f(M * nUnique, nUnique, nTotal) > 0) {
    M <- M*10.0
  }
  
  for(i in seq(0, 40)) {
    r <-  (m + M) / 2.0
    u <- f(r * nUnique, nUnique, nTotal);
    if (u == 0) {
      break
    } else if (u > 0) {
      m = r
    } else if (u < 0) {
      M = r
    }
  }
  
  return(round(nUnique * (m + M) / 2.0))
}

estimateDuplicateRate <- function(nTotal, nUnique){
  round((nTotal-nUnique)/nTotal, 3)
}

#-----------------
# Parse arguments
#-----------------

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}


bamFile <- args[i+1] # directory of .rds files
barcodeTranslateFile <- args[i+2] # file path to the number of barcodes for each observed barcode
barcodeQuantsFile <- args[i+3] # filepath to write the implicated barcode pairs
tssFile <- args[i+4]
tag <- args[i+5]
blacklistFile <- args[i+6]

# For devel only
if(FALSE){
  bamFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.bap.bam"
  barcodeTranslateFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.barcodeTranslate.tsv"
  barcodeQuantsFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.barcodequants.csv"
  tssFile <-  "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/anno/TSS/hg19.refGene.TSS.bed"
  tag <- "DB"
  blacklistFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/anno/blacklist/hg19.full.blacklist.bed"
}

# Parse .bam file
GA <- readGAlignments(bamFile, param = ScanBamParam(
    flag = scanBamFlag(isMinusStrand = FALSE, isProperPair = TRUE),
    tag = c(tag), mapqFilter = 0, what="isize"))

# Set up TSS file for scoring
df <- data.frame(fread(tssFile))
df$V2 <- df$V2 - 2000
df$V3 <- df$V3 + 2000


# Deal with TSS enrichment %
promoter <- makeGRangesFromDataFrame(df, seqnames.field = "V1", start.field = "V2", end.field = "V3")
ov <- findOverlaps(promoter, GA)
df <- data.frame(mcols(GA))
df$TSS <- as.numeric(1:dim(df)[1] %in% subjectHits(ov))
rm(ov); rm(GA); rm(promoter)
colnames(df) <- c("isize", "DropBarcode", "TSS")
sbdf <- df %>% group_by(DropBarcode) %>% summarise(tssPproportion = round(mean(TSS), 3),
                                                   meanInsertSize = round(mean(isize)),
                                                   medianInsertSize = round(median(isize))) %>% data.frame()
rm(df)

# Import old barcode stats
mss <- data.frame(merge(fread(barcodeTranslateFile), fread(barcodeQuantsFile), by.x = "BeadBarcode", by.y = "Barcode"))
msss <- mss %>% group_by(DropBarcode) %>% summarise(uniqueNuclearFrags = sum(UniqueNuclear), 
                                                    totalNuclearFrags = sum(TotalNuclear),
                                                    totalMitoFrags = sum(TotalMito)) %>% data.frame()

msss %>% mutate(duplicateProportion = estimateDuplicateRate(totalNuclearFrags, uniqueNuclearFrags),
                librarySize = estimateLibrarySize(totalNuclearFrags, uniqueNuclearFrags)) -> msss
qcStats <- merge(sbdf, msss)

write.table(qcStats,
            file = gsub(".barcodequants.csv$", ".QCstats.csv", barcodeQuantsFile),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")


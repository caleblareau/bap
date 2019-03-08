options(warn=-1)

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
  
  nDuplicates <- (nTotal - nUnique) + 1 # charity to handle only unique reads observed
  
  if (nUnique > nTotal | (f(m * nUnique, nUnique, nTotal) < 0) | nUnique < 0 | nTotal < 0 | nDuplicates < 0) {
    message("Library size returns 0 -- invalid inputs; check this cell more closely")
    return(0)
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
dropbarcodeTag <- args[i+5]
blacklistFile <- args[i+6]
peakFile <- args[i+7]
speciesMix <- args[i+8]
oneone <- args[i+9]

# For devel only
if(FALSE){
  bamFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.bap.bam"
  barcodeTranslateFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.barcodeTranslate.tsv"
  barcodeQuantsFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.barcodequants.csv"
  tssFile <-  "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/anno/TSS/hg19.refGene.TSS.bed"
  dropbarcodeTag <- "DB"
  blacklistFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/anno/blacklist/hg19.full.blacklist.bed"
  peakFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/data/test.small.peaks.bed"
  speciesMix <- "no"
  oneone <- "yes"
}

# We only need the read names if we have a valid peak file
import.names <- peakFile != "none"

# Parse .bam file
GA <- readGAlignments(bamFile, use.names = import.names, param = ScanBamParam(
  flag = scanBamFlag(isProperPair = TRUE),
  tag = c(dropbarcodeTag), mapqFilter = 0, what="isize"))

# Set up TSS file for scoring
tssdf <- data.frame(fread(tssFile))
tssdf$V2 <- tssdf$V2 - 2000
tssdf$V3 <- tssdf$V3 + 2000

# Deal with TSS enrichment %
promoter <- makeGRangesFromDataFrame(tssdf, seqnames.field = "V1", start.field = "V2", end.field = "V3")
ovTSS <- findOverlaps(promoter, GA)
df <- data.frame(mcols(GA))

df$TSS <- as.numeric(1:dim(df)[1] %in% subjectHits(ovTSS))
rm(ovTSS); rm(promoter)

# Deal with FRIP if we have a peak file
if(peakFile != "none"){
  peakdf <- data.frame(fread(peakFile))[,c(1,2,3)]
  colnames(peakdf) <- c("V1", "V2", "V3")
  peaks <- makeGRangesFromDataFrame(peakdf, seqnames.field = "V1", start.field = "V2", end.field = "V3")
  ovPEAK <- findOverlaps(peaks, GA)
  df$Peak <- as.numeric(1:dim(df)[1] %in% subjectHits(ovPEAK))
  rm(peakdf)
} else{
  df$Peak <- 0
}


if(speciesMix == "no"){
  colnames(df) <- c("isize", "DropBarcode", "TSS", "Peak")
  sbdf <- df %>% group_by(DropBarcode) %>% summarise(FRIP = round(mean(Peak), 3),
                                                     tssProportion = round(mean(TSS), 3),
                                                     meanInsertSize = round(mean(abs(isize))),
                                                     medianInsertSize = round(median(abs(isize)))) %>% data.frame()
} else {
  df$mouse <- as.numeric(substr(as.character(seqnames(GA)),1,1) == "m")
  df$human <- as.numeric(substr(as.character(seqnames(GA)),1,1) == "h")
  colnames(df) <- c("isize", "DropBarcode", "TSS", "Peak", "mouse", "human")
  sbdf <- df %>% group_by(DropBarcode) %>% summarise(FRIP = round(mean(Peak), 3),
                                                     tssProportion = round(mean(TSS), 3),
                                                     meanInsertSize = round(mean(abs(isize))),
                                                     medianInsertSize = round(median(abs(isize))),
                                                     nHumanReads = sum(human),
                                                     nMouseReads = sum(mouse)) %>% data.frame()
}
rm(df)

# Import old barcode stats and merge with what has been computed
mss <- data.frame(merge(fread(barcodeTranslateFile), fread(barcodeQuantsFile), by.x = "BeadBarcode", by.y = "Barcode"))

# Also subtract out the total number of overlap reads (divided by 2 sincel double counting) to get proper estimate of # 
# of unique nuclear reads. This technically will be a slight UNDERestimate in cases of 3 or more drops merging where
# all 3 share the same read, but we assume that this is a weird edge case that infrequently happens and thus is not
# of the utmost importance in handling here. 
msss <- mss %>% group_by(DropBarcode) %>% summarise(uniqueNuclearFrags = round(sum(UniqueNuclear) - (sum(OverlapReads)/2)), 
                                                    totalNuclearFrags = sum(TotalNuclear),
                                                    totalMitoFrags = sum(TotalMito),
                                                    totalNCfilteredFrags = sum(TotalNC)) %>% data.frame()

msss$duplicateProportion <- estimateDuplicateRate(msss$totalNuclearFrags, msss$uniqueNuclearFrags)
msss$librarySize <- sapply(1:dim(msss)[1], function(i){
  estimateLibrarySize(msss[i,"totalNuclearFrags"],msss[i,"uniqueNuclearFrags"])})
qcStats <- merge(sbdf, msss)

# add barcodes back if we need to (specified with the one to one option)
if(oneone == "yes"){
  translateDF <- read.table(barcodeTranslateFile, header = TRUE, stringsAsFactors = FALSE)
  qcStatsDF <- qcStats
  tvec <- translateDF[,1]; names(tvec) <- translateDF[,2]
  qcStats <- qcStatsDF
  qcStats$OriginalBarcode <- tvec[as.character(qcStats$DropBarcode)] 
}

if(peakFile != "none"){
  # recall ovPEAK from above 

  id <- factor(mcols(GA)[,dropbarcodeTag], levels = as.character(qcStats$DropBarcode))
  countdf <- data.frame(peak = queryHits(ovPEAK), sample = as.numeric(id)[subjectHits(ovPEAK)], read = names(GA)[subjectHits(ovPEAK)]) %>%
    distinct() %>%  # by filtering on distinct read / peak / sample trios, ensure that PE reads that overlap peak aren't double counted
    select(-one_of("read")) %>% 
    group_by(peak,sample) %>% summarise(count = n()) %>% data.matrix()
  #rm(ovPEAK)
  
  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks)), j =c(countdf[,2], 1), x = c(countdf[,3],0))
  colnames(m) <- levels(id)
  SE <- SummarizedExperiment(
    rowRanges = peaks, 
    assays = list(counts = m),
    colData = qcStats
  )
  saveRDS(SE, gsub(".barcodequants.csv$", ".SE.rds", barcodeQuantsFile))
}

rm(GA)

write.table(qcStats,
            file = gsub(".barcodequants.csv$", ".QCstats.csv", barcodeQuantsFile),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")

# Reformat the Barcode Translate file to remove the spurious additional column that was needed
write.table(data.frame(fread(barcodeTranslateFile))[,c(1,2)],
            file = barcodeTranslateFile, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")

options(warn=-1)

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

fragsFile <- args[i+1] 
tssFile <- args[i+2]
peakFile <- args[i+3]
simpleQCfile <- args[i+4] 
speciesMix <- args[i+5]
oneone <- args[i+6]
barcodeTranslateFile <- args[i+7]

# For dev only
if(FALSE){
  fragsFile <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/final/jaccardPairsForIGV.fragments.tsv.gz"
  tssFile <-  "~/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/anno/TSS/hg19.refGene.TSS.bed"
  peakFile <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/data/test.small.peaks.bed"
  simpleQCfile <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/final/jaccardPairsForIGV.basicQC.tsv"
  speciesMix <- "no"
  oneone <- "no"
  barcodeTranslateFile <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/final/jaccardPairsForIGV.barcodeTranslate.tsv"
}

# Import fragments
frag_df <- fread(fragsFile, header = FALSE)
if(dim(frag_df)[2] == 4){
  colnames(frag_df) = c("chr", "start", "end", "cell_barcode")
} else {
  colnames(frag_df) = c("chr", "start", "end", "cell_barcode", "read_name")
}
frag_gr <- makeGRangesFromDataFrame(frag_df)

# Set up TSS file for scoring
tssdf <- data.frame(fread(tssFile))
tssdf$V2 <- tssdf$V2 - 1000
tssdf$V3 <- tssdf$V3 + 1000

# Deal with TSS enrichment %
promoter <- makeGRangesFromDataFrame(tssdf, seqnames.field = "V1", start.field = "V2", end.field = "V3")
ovTSS <- findOverlaps(promoter, frag_gr)

frag_df$TSS <- as.numeric(1:dim(frag_df)[1] %in% subjectHits(ovTSS))
rm(ovTSS); rm(promoter)

# Deal with FRIP if we have a peak file
if(peakFile != "none"){
  suppressMessages(suppressWarnings(library(SummarizedExperiment)))
  peakdf <- data.frame(fread(peakFile))[,c(1,2,3)]
  colnames(peakdf) <- c("V1", "V2", "V3")
  peaks <- makeGRangesFromDataFrame(peakdf, seqnames.field = "V1", start.field = "V2", end.field = "V3")
  ovPEAK <- findOverlaps(peaks, frag_gr)
  frag_df$Peak <- as.numeric(1:dim(frag_df)[1] %in% subjectHits(ovPEAK))
  rm(peakdf) # Keep ovPEAK around for a while in case we build a SE
} else{
  frag_df$Peak <- 0
}

# Summarize frag attributes
frag_df$insert_size <- frag_df$end - frag_df$start
frag_df_summary <- frag_df[, list(meanInsertSize = mean(insert_size),
                                  medianInsertSize = median(insert_size),
                                  tssProportion = mean(TSS),
                                  FRIP = mean(Peak)), by = cell_barcode]

#Import existing QC stats
qc_df <- fread(simpleQCfile, header = TRUE)
qc_df$duplicateProportion <- estimateDuplicateRate(qc_df$totalNuclearFrags, qc_df$uniqueNuclearFrags)
qc_df$librarySize <- sapply(1:dim(qc_df)[1], function(i){
  estimateLibrarySize(as.numeric(qc_df[i,"totalNuclearFrags"]),as.numeric(qc_df[i,"uniqueNuclearFrags"]))})
qcStats <- merge(qc_df, frag_df_summary)

# Make it look pretty
qcStats$meanInsertSize <- round(qcStats$meanInsertSize, 1)
qcStats$duplicateProportion <- round(qcStats$duplicateProportion, 4)
qcStats$tssProportion <- round(qcStats$tssProportion, 4)
qcStats$FRIP <- round(qcStats$FRIP, 4)
colnames(qcStats)[1] <- "DropBarcode"
qcStats <- data.frame(qcStats)

# add barcodes back if we need to (specified with the one to one option)
if(oneone == "yes"){
  translateDF <- read.table(barcodeTranslateFile, header = TRUE, stringsAsFactors = FALSE)
  tvec <- translateDF[,1]; names(tvec) <- translateDF[,2]
  qcStats$OriginalBarcode <- tvec[as.character(qcStats$DropBarcode)] 
}

# Create a Summarized Experiment if we have a valid peak file
if(peakFile != "none"){
  
  # recall ovPEAK from above 
  id <- factor(frag_df[["cell_barcode"]], levels = as.character(qcStats$DropBarcode))
  countdf <- data.frame(peak = queryHits(ovPEAK), sample = as.numeric(id)[subjectHits(ovPEAK)]) %>%
    group_by(peak,sample) %>% summarise(count = n()) %>% data.matrix()
  #rm(ovPEAK)
  
  m <- Matrix::sparseMatrix(i = c(countdf[,1], length(peaks)), j= c(countdf[,2], 1), x = c(countdf[,3],0))
  colnames(m) <- levels(id)
  SE <- SummarizedExperiment(
    rowRanges = peaks, 
    assays = list(counts = m),
    colData = qcStats
  )
  saveRDS(SE, gsub(".basicQC.tsv$", ".SE.rds", simpleQCfile))
}

# Either way, export QC stats
write.table(qcStats,
            file = gsub(".basicQC.tsv$", ".QCstats.csv", simpleQCfile),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")


options(warn=-1)

suppressMessages(suppressWarnings(library(tools)))

args <- commandArgs(trailingOnly = TRUE)

if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}

# barcodeTranslateFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.barcodeTranslate.tsv"
barcodeTranslateFile <- args[i+1] # directory of .rds files
qcStatsFile <- gsub(".barcodeTranslate.tsv$", ".QCstats.csv", barcodeTranslateFile)

translateDF <- read.table(barcodeTranslateFile, header = TRUE, stringsAsFactors = FALSE)
qcStatsDF <- read.table(qcStatsFile, header = TRUE, sep = ",")
tvec <- translateDF[,1]; names(tvec) <- translateDF[,2]
qcStats <- qcStatsDF
qcStats$DropBarcode <- tvec[as.character(qcStats$DropBarcode)]

write.table(qcStats,
            file = qcStatsFile,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
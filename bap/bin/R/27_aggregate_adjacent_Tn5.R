options(warn=-1)

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))

"%ni%" <- Negate("%in%")

options(scipen=999)

# This function / script for processed chromosome files to produce consensus barcode doublets
# as well as summary and QC metrics and visuals

substrRight <- function(x, n = 6){
  substr(x, nchar(x)-n+1, nchar(x))
}

# I/O
args <- commandArgs(trailingOnly = FALSE)

# Grab knee call functions from figuring out where this script lives and then changing the file name
needle <- "--file="
match <- grep(needle, args)
source(normalizePath(gsub("27_aggregate_adjacent_Tn5.R", "00_knee_CL.R", gsub(needle, "", args[match]))))

# Import parameters using logic from the end3
nn <- length(args)
rdsDir <- args[nn-1] # directory of .rds files
outputRDS <- args[nn] # directory of .rds files

# For devel only
if(FALSE){
  rdsDir <- "~/dat/Research/Boston/lareau_dev/bap/tests/bap2/temp/frag_overlap"
  outputRDS <- "~/dat/Research/Boston/lareau_dev/bap/tests/bap2/final/jaccardPairsForIGV.adjacentTn5.rds"
}

rdsFiles <- list.files(rdsDir, full.names = TRUE, pattern = "_adjacentBpCount.rds$")
rdsFiles <- rdsFiles[sapply(lapply(rdsFiles,file.info), function(x) x$size) > 0]
lapply(rdsFiles, readRDS) %>%
  rbindlist(fill = TRUE) -> inputDF

saveRDS(
  inputDF[,.(total = sum(count)), by=list(barc1, barc2, distance_computed)],
  outputRDS)

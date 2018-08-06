options(warn=-1)

suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(tools)))

"%ni%" <- Negate("%in%")

options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}

out_directory <- args[i+1] # directory of final output
keep_temp_files <- args[i+2] # whether or not to keep temporary files
barcodeQuantsFile <- args[i+3]

# For devel only
if(FALSE){
  out_directory <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out"
  keep_temp_files <= "True"
  barcodeQuantsFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.barcodequants.csv"
}

possible_ncRead <- list.files(paste0(out_directory, "/temp/", "filt_split"), full.names = TRUE, pattern = "*_ncRead.tsv$")
if(length(possible_ncRead) > 0){
  
  # Filter out files that have zero size
  possible_ncRead <- possible_ncRead[file.info(possible_ncRead)$size > 0]
  
  lapply(possible_ncRead, function(file){
    data.frame(table(fread(file)[[2]]))
  }) %>% rbindlist() %>% data.frame() -> count_ncDF
  
  # Aggregate the number of counts of NC reads across chromosomes
  count_ncDF$Var1 <- as.numeric(as.character(count_ncDF$Var1))
  count_ncDF %>% group_by(Var1) %>% summarise(Total = sum(Freq)) %>% arrange(Var1) %>% data.frame() -> out_count
  colnames(out_count) <- c("NC_value", "NumberOfReads")
  
  write.table(out_count,
              file = gsub(".barcodequants.csv$", ".NCsumstats.tsv", barcodeQuantsFile),
              quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")
  
  rm(count_ncDF)
  
  # Possibly remove files at this point
  if(keep_temp_files != "True"){
    toke <- lapply(possible_ncRead, file.remove)
  }
}

# Augment quants file
quants <- data.frame(fread(barcodeQuantsFile))
quants$TotalNC <- 0

possible_ncBarcode <- list.files(paste0(out_directory, "/temp/", "filt_split"), full.names = TRUE, pattern = "*_ncBarcode.tsv$")
if(length(possible_ncBarcode) > 0){
  
  # Filter out files that have zero size
  possible_ncBarcode <- possible_ncBarcode[file.info(possible_ncBarcode)$size > 0]
  
  # Import across all 
  lapply(possible_ncBarcode, function(file){
    data.frame(fread(file))
  }) %>% rbindlist() %>% data.frame()  %>% group_by(V1) %>% summarise(totalNC = sum(V2)) %>% data.frame() -> count_bcDF
  vec <- as.numeric(count_bcDF$totalNC); names(vec) <- as.character(count_bcDF$V1)
  
  # Update values in quants if present here
  quants$TotalNC <- sapply(1:dim(quants), function(j){
    bc <- as.character(quants[j,"Barcode"])
    ifelse(bc %in% names(vec), vec[bc], 0)
  })
  
  if(keep_temp_files != "True"){
    toke <- lapply(possible_ncBarcode, file.remove)
  }
}

# Write back out with new column
write.table(quants,
            file = barcodeQuantsFile,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")



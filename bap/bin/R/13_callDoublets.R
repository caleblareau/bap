library(data.table)
library(dplyr)
library(ggplot2)
library(tools)
"%ni%" <- Negate("%in%")

options(scipen=999)

# This function / script for processed chromosome files to produce consensus barcode doublets
# as well as summary and QC metrics and visuals

# TO DO:
# actually call doublets -- make another additional text file
# produce QC plot
# make barcode key/value

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}

rdsDir <- args[i+1] # directory of .rds files
nbcin <- args[i+2] # file path to the number of barcodes for each observed barcode
tblOut <- args[i+3] # filepath to write the implicated barcode pairs
min_jaccard_frag <- as.numeric(args[i+4])
name <- args[i+5] #name prefix for file naming convention
one_to_one <- args[i+6] #arguement for keeping bead : drop conversion 1 to 1
one_to_one <- one_to_one == "True" # Fix to R boolean


# For devel only
if(FALSE){
  rdsDir <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/frag_overlap"
  nbcin <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.barcodequants.csv"
  tblOut <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.implicatedBarcodes.csv"
  name <- "test.small"
  min_jaccard_frag <- 0.005
  one_to_one <- FALSE
}
rdsFiles <- list.files(rdsDir, full.names = TRUE)


# Import the number of barcodes per bead
nBC <- data.frame(fread(nbcin, header = TRUE)) %>% arrange(desc(UniqueNuclear)); colnames(nBC) <- c("BeadBarcode", "UniqueFragCount", "TotalFragCount", "TotalMitoCount")
nBC_keep <- nBC; nBC_keep$DropBarcode <- ""
nBCv <- nBC$UniqueFragCount; names(nBCv) <- nBC$BeadBarcode

# Implicate overlapping fragments
lapply(rdsFiles, readRDS) -> o
o[sapply(o, is.data.frame)] -> o
o %>%
  rbindlist(fill = TRUE) %>% as.data.frame() %>%
  group_by(barc1, barc2) %>%
  summarise(N_both = sum(n_both)/2) %>% # overlaps called twice
  filter(N_both > 1) %>% mutate(N_barc1 = nBCv[barc1], N_barc2 = nBCv[barc2]) %>% 
  mutate(jaccard_frag = round((N_both)/(N_barc1 + N_barc2 - N_both + 1),4)) %>% filter(jaccard_frag > min_jaccard_frag) %>% 
  arrange(desc(jaccard_frag)) %>% data.frame() -> ovdf

write.table(ovdf, file = tblOut,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")

# Guess at how wide we need to make the barcodes to handle leading zeros
guess <- ceiling(log10(dim(nBC)[1]))

idx <- 1

# Loop through and eat up barcodes
while(dim(nBC)[1] > 0){
  barcode <- as.character(nBC[1,1])
  barcode_combine <- barcode
  
  # Find friends that are similarly implicated and append from Barcode 1
  friendsRow1 <- which(barcode ==  ovdf[,"barc1", drop = TRUE])
  if(length(friendsRow1) > 0){
    friends1 <- as.character(ovdf[friendsRow1,"barc2"])
    ovdf <- ovdf[1:dim(ovdf)[1] %ni% friendsRow1,]
    barcode_combine <- c(barcode_combine, friends1)
  }
  
  # Find friends that are similarly implicated and append from Barcode 2
  friendsRow2 <- which(barcode ==  ovdf[,"barc2", drop = TRUE])
  if(length(friendsRow2) > 0){
    friends2 <- as.character(ovdf[friendsRow2,"barc1"])
    ovdf <- ovdf[1:dim(ovdf)[1] %ni% friendsRow2,]
    barcode_combine <- c(barcode_combine, friends2)
  }
  
  # If user species one to one, then only remove that one barcode
  if(one_to_one){
    barcode_combine <- barcode
  }
  
  # Make a drop barcode and save our progress
  dropBarcode <- paste0(name, "_BC", formatC(idx, width=guess, flag="0", digits = 20), "_N", sprintf("%02d", length(barcode_combine)))
  
  nBC_keep[nBC_keep$BeadBarcode %in% barcode_combine, "DropBarcode"] <- dropBarcode
  idx <- idx + 1
  
  # Remove barcodes that we've dealt with
  nBC <- nBC[nBC$BeadBarcode %ni% barcode_combine,]
}

write.table(nBC_keep[,c("BeadBarcode", "DropBarcode")],
            file = gsub(".implicatedBarcodes.csv$", ".barcodeTranslate.tsv", tblOut),
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")


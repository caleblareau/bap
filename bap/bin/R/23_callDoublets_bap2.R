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
source(normalizePath(gsub("23_callDoublets_bap2.R", "00_knee_CL.R", gsub(needle, "", args[match]))))

# Import parameters using logic from the end3
nn <- length(args)
rdsDir <- args[nn-8] # directory of .rds files
n_bc_file <- args[nn-7] # file for the number of reads supporting each barcode
hq_bc_file <- args[nn-6] # file for the HQ barcodes that were nominated
tblOut <- args[nn-5] # filepath to write the implicated barcode pairs
min_jaccard_frag <- as.numeric(args[nn-4])
name <- args[nn-3] #name prefix for file naming convention
one_to_one <- args[nn-2] #arguement for keeping bead : drop conversion 1 to 1
barcoded_tn5 <- args[nn-1]
barcode_prior_file <- args[nn]

# Fix to R boolean
one_to_one <- one_to_one == "True" 
barcoded_tn5 <- barcoded_tn5 == "True" 

# Replace the .gz convention
tblOut <- gsub(".gz$", "", tblOut)

# For devel only
if(FALSE){
  rdsDir <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/temp/frag_overlap"
  n_bc_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/knee/jaccardPairsForIGV.barcodeQuantSimple.csv"
  hq_bc_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/final/jaccardPairsForIGV.HQbeads.tsv"
  tblOut <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/final/o.tsv"
  min_jaccard_frag <- 0.005
  name <- "x"
  one_to_one <- FALSE
  barcoded_tn5 <- FALSE
  barcode_prior_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/data/jaccardPairsTest_sep.tsv"
  
}

rdsFiles <- list.files(rdsDir, full.names = TRUE, pattern = "_overlapCount.rds$")
rdsFiles <- rdsFiles[sapply(lapply(rdsFiles,file.info), function(x) x$size) > 0]


lapply(rdsFiles, readRDS) %>%
  rbindlist(fill = TRUE) -> inputDF

# Only consider merging when Tn5 is the same
if(barcoded_tn5){
  tn5_1 = substrRight(inputDF[["barc1"]])
  tn5_2 = substrRight(inputDF[["barc2"]])
  inputDF <- inputDF[tn5_1 == tn5_2]
}


# Import number of barcodes
valid_barcodes <- fread(hq_bc_file, col.names = c("bc"), header = FALSE)[["bc"]]
nBC <- fread(n_bc_file, col.names = c("BeadBarcode", "count"), sep = ",") %>%
  data.frame() %>% filter(BeadBarcode %in% valid_barcodes) %>% arrange(desc(count))
count_vec <- nBC$count*2; names(count_vec) <- as.character(nBC$BeadBarcode)

sum_dt <- inputDF[, .(N_both = sum(n_both)), by = list(barc1, barc2)] 
sum_dt$N_barc1 <- count_vec[sum_dt$barc1]
sum_dt$N_barc2 <- count_vec[sum_dt$barc2]

data.frame(sum_dt) %>% # fixed the previous divided by two in the upstream script (22) for overall accuracy
  mutate(jaccard_frag = round((N_both)/(N_barc1 + N_barc2 - N_both + 0.05), 5)) %>% 
  filter(jaccard_frag > 0) %>% 
  arrange(desc(jaccard_frag)) %>% data.frame() -> ovdf


# Call knee if we need to
if(min_jaccard_frag == 0){
  message("Computing jaccard index for bead merging via a knee call--")
  two <- get_density_threshold_CL(head(ovdf$jaccard_frag, 1000000), "jaccard", logTransform = TRUE)
  min_jaccard_frag <- two[1]
  called_jaccard_frag <- two[2]
  
  # Write out what gets called
  write(paste0("jaccard_threshold_nosafety,", as.character(called_jaccard_frag)),
        gsub("/final/", "/knee/", gsub(".implicatedBarcodes.csv$", ".bapParams.csv", tblOut)), append = TRUE)
  
}

# Append to bap parameters
write(paste0("jaccard_threshold,", as.character(min_jaccard_frag)),
      gsub("/final/", "/knee/", gsub(".implicatedBarcodes.csv$", ".bapParams.csv", tblOut)), append = TRUE)

ovdf$merged <- ovdf$jaccard_frag > min_jaccard_frag

# Merge barcodes above the knee unless there's a prior reason why we shouldn't
if(barcode_prior_file != "none"){
  
  # Import and assign each 
  bp_df <- fread(barcode_prior_file, header = FALSE) 
  feature_vec <- as.character(bp_df[["V2"]]); names(feature_vec) <- as.character(bp_df[["V1"]])
  merged_df <- ovdf %>% filter(merged)
  priorf1 <- feature_vec[as.character(merged_df[["barc1"]])] %>% unname()
  priorf2 <- feature_vec[as.character(merged_df[["barc2"]])] %>% unname()
  
  # Find conflicts,  missing values, and 
  conflicts <- priorf1 != priorf2
  sum_is_na <- sum(is.na(conflicts))
  sum_conflicts = sum(conflicts, na.rm = TRUE)
  sum_valid = sum(!conflicts, na.rm = TRUE)
  
  # Write out conflicts and filter if we observe them
  if(sum_conflicts > 0){
    conflictFile = gsub("/final/", "/knee/", gsub(".implicatedBarcodes.csv$", ".barcode_prior_conflicts.csv", tblOut))
    write.table(merged_df[which(conflicts),], file = conflictFile, row.names = FALSE, col.names = TRUE, sep = ",", quote = FALSE)
    ovdf <- ovdf[1:dim(ovdf)[1] %ni% which(conflicts),]
  }
  
  # Either way, write out statistics
  out_stat_prior_df <- data.frame(
    what = c("Valid merges (stil merged)", "Merges with 1 or more NA (still merged)", "Conflicted merges (not merged)"),
    count = c(sum_valid, sum_is_na, sum_conflicts)
  )
  statFile = gsub("/final/", "/knee/", gsub(".implicatedBarcodes.csv$", ".barcode_prior_stats.csv", tblOut))
  write.table(out_stat_prior_df, file = statFile, row.names = FALSE, col.names = FALSE, sep = ",", quote = FALSE)
  
} 


# Export the implicated barcodes
write.table(ovdf, file = tblOut,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
system(paste0("gzip ", tblOut))

# Now filter based on the min_jaccard_frag
ovdf %>% filter(merged) %>% data.frame() -> ovdf

# Guess at how wide we need to make the barcodes to handle leading zeros
guess <- ceiling(log10(dim(nBC)[1]))

nBC_keep <- nBC; nBC_keep$DropBarcode <- ""; nBC_keep$OverlapReads <- ""

# Loop through and eat up barcodes
idx <- 1
while(dim(nBC)[1] > 0){
  barcode <- as.character(nBC[1,1])
  barcode_combine <- barcode
  OverlapReads = 0
  
  # Find friends that are similarly implicated and append from Barcode 1
  friendsRow1 <- which(barcode ==  ovdf[,"barc1", drop = TRUE])
  if(length(friendsRow1) > 0){
    friends1 <- as.character(ovdf[friendsRow1,"barc2"])
    #OverlapReads = OverlapReads + sum(ovdf[1:dim(ovdf)[1] %in% friendsRow1, "N_both"])
    #ovdf <- ovdf[1:dim(ovdf)[1] %ni% friendsRow1,]
    barcode_combine <- c(barcode_combine, friends1)
  }
  
  # Find friends that are similarly implicated and append from Barcode 2
  friendsRow2 <- which(barcode ==  ovdf[,"barc2", drop = TRUE])
  if(length(friendsRow2) > 0){
    friends2 <- as.character(ovdf[friendsRow2,"barc1"])
    #OverlapReads = OverlapReads + sum(ovdf[1:dim(ovdf)[1] %in% friendsRow2, "N_both"])
    #ovdf <- ovdf[1:dim(ovdf)[1] %ni% friendsRow2,]
    barcode_combine <- c(barcode_combine, friends2)
  }
  
  # If user species one to one, then only remove that one barcode
  if(one_to_one) barcode_combine <- barcode
  
  # Make a drop barcode and save our progress
  if(!barcoded_tn5){
    dropBarcode <- paste0(name, "_BC", formatC(idx, width=guess, flag="0", digits = 20), "_N", sprintf("%02d", length(barcode_combine)))
  } else {
    dropBarcode <- paste0(name, "_Tn5-", substrRight(barcode_combine), 
                          "_BC", formatC(idx, width=guess, flag="0", digits = 20),
                          "_N", sprintf("%02d", length(barcode_combine)))
  }
  
  # Annotate with new values
  nBC_keep[nBC_keep$BeadBarcode %in% barcode_combine, "DropBarcode"] <- dropBarcode
  #nBC_keep[nBC_keep$BeadBarcode %in% barcode_combine, "OverlapReads"] <- OverlapReads
  
  idx <- idx + 1
  
  # Remove barcodes that we've dealt with
  nBC <- nBC[nBC$BeadBarcode %ni% barcode_combine,]
}

write.table(nBC_keep[,c("BeadBarcode", "DropBarcode")],
            file = gsub(".implicatedBarcodes.csv$", ".barcodeTranslate.tsv", tblOut),
            quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")


# Finally, collate the NC values
tsvFiles <- list.files(rdsDir, full.names = TRUE, pattern = "_ncCount.tsv$")
tsvFiles <- tsvFiles[sapply(lapply(tsvFiles,file.info), function(x) x$size) > 0]

lapply(tsvFiles, read.table, header = FALSE) %>%
  rbindlist(fill = TRUE) %>% data.frame() -> inputDF
colnames(inputDF) <- c("NC_value", "i")
inputDF %>% group_by(NC_value) %>% 
  summarize(NumberOfFragments = sum(i)) -> NCtsv_out_df

nc_file_out <- gsub(".implicatedBarcodes.csv$", ".NCsumstats.tsv", tblOut)
write.table(NCtsv_out_df, file = nc_file_out, 
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")



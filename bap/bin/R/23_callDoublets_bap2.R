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
rdsDir <- args[nn-7] # directory of .rds files
n_bc_file <- args[nn-6] # file for the number of reads supporting each barcode
hq_bc_file <- args[nn-5] # file for the HQ barcodes that were nominated
tblOut <- args[nn-4] # filepath to write the implicated barcode pairs
min_jaccard_frag <- as.numeric(args[nn-3])
name <- args[nn-2] #name prefix for file naming convention
one_to_one <- args[nn-1] #arguement for keeping bead : drop conversion 1 to 1
barcoded_tn5 <- args[nn]

# Fix to R boolean
one_to_one <- one_to_one == "True" 
barcoded_tn5 <- barcoded_tn5 == "True" 

# Replace the .gz convention
tblOut <- gsub(".gz$", "", tblOut)

# For devel only
if(FALSE){
  rdsDir <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/frag_overlap"
  n_bc_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/knee/jaccardPairsForIGVbarcodeQuantsSimple.csv"
  tblOut <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/o.tsv"
  min_jaccard_frag <- 0.005
  name <- "x"
  one_to_one <- FALSE
  barcoded_tn5 <- FALSE
}

if(FALSE){
  rdsDir <- "/data/aryee/caleb/biorad/mouse_brain/N704_Exp119_sample4_S1_2b2a2p/temp/frag_overlap/"
  n_bc_file <- "/data/aryee/caleb/biorad/mouse_brain/N704_Exp119_sample4_S1_2b2a2p/knee/N704_Exp119_sample4_S1barcodeQuantsSimple.csv"
  hq_bc_file <- "/data/aryee/caleb/biorad/mouse_brain/N704_Exp119_sample4_S1_2b2a2p/final/N704_Exp119_sample4_S1.HQbeads.tsv"
  tblOut <- "/data/aryee/caleb/biorad/mouse_brain/N704_Exp119_sample4_S1_2b2a2p/final/o.tsv"
  min_jaccard_frag <- 0.005
  name <- "x"
  one_to_one <- FALSE
  barcoded_tn5 <- FALSE
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

inputDF[, .(N_barc1 = sum(n_barc1), N_barc2 = sum(n_barc2), N_both = sum(n_both)), by = list(barc1, barc2)] %>%
  data.frame() %>% # fixed the previous divided by two in the upstream script (22) for overall accuracy
  mutate(jaccard_frag = round((N_both)/(N_barc1 + N_barc2 - N_both + 1),4)) %>% 
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

# Export the implicated barcodes
ovdf$merged <- ovdf$jaccard_frag > min_jaccard_frag
write.table(ovdf, file = tblOut,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")
system(paste0("gzip ", tblOut))

# Now filter based on the min_jaccard_frag
ovdf %>% filter(jaccard_frag > min_jaccard_frag) %>% data.frame() -> ovdf

# Import number of barcodes
valid_barcodes <- fread(hq_bc_file, col.names = c("bc"))[["bc"]]
nBC <- fread(n_bc_file, col.names = c("BeadBarcode", "count"), sep = ",") %>%
  data.frame() %>% filter(BeadBarcode %in% valid_barcodes) %>% arrange(desc(count))

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


library(data.table)
library(BuenColors)
library(dplyr)

"%ni%" <- Negate("%in%")
substrRight <- function(x, n)  substr(x, nchar(x)-n+1, nchar(x))

# Import barcode quants and pairwise jaccard frag overlap files
ovdf_all <- data.frame(fread(paste0("zcat < ", "N714_Exp68_sample11_S1.implicatedBarcodes.csv.gz")))
nBC <- data.frame(fread("N714_Exp68_sample11_S1.barcodequants.csv"))


# Guess at how wide we need to make the barcodes to handle leading zeros
guess <- ceiling(log10(dim(nBC)[1]))

# Function to take these two data frames and a parameterized jaccard
# and returns the proportion of doublets or some other merge attribute
#
# TO DO: run a knee-call over this to make the denominator smaller /
# only have 'cells' that pass knee
determineMergeProportion <- function(min_jaccard_threshold){
  
  # Print what threshold we are looking at
  print(min_jaccard_threshold)
  
  # Filter paired observations such that they have at least
  # a sufficient jaccard overlap propriton
  ovdf_all %>% filter(jaccard_frag > min_jaccard_threshold) %>% data.frame() -> ovdf
  
  # Initialize other new data framges
  nBC_keep <- nBC; nBC_keep$DropBarcode <- ""; nBC_keep$OverlapReads <- ""
  
  # Loop through and eat up barcodes
  idx <- 1
  name <- "TestMerge"
  while(dim(nBC)[1] > 0){
    
    # Set up a new barcode
    barcode <- as.character(nBC[1,1])
    barcode_combine <- barcode
    OverlapReads = 0
    
    # Find friends that are similarly implicated and append from Barcode 1
    friendsRow1 <- which(barcode ==  ovdf[,"barc1", drop = TRUE])
    if(length(friendsRow1) > 0){
      friends1 <- as.character(ovdf[friendsRow1,"barc2"])
      OverlapReads = OverlapReads + sum(ovdf[1:dim(ovdf)[1] %in% friendsRow1, "N_both"])
      ovdf <- ovdf[1:dim(ovdf)[1] %ni% friendsRow1,]
      barcode_combine <- c(barcode_combine, friends1)
    }
    
    # Find friends that are similarly implicated and append from Barcode 2
    friendsRow2 <- which(barcode ==  ovdf[,"barc2", drop = TRUE])
    if(length(friendsRow2) > 0){
      friends2 <- as.character(ovdf[friendsRow2,"barc1"])
      OverlapReads = OverlapReads + sum(ovdf[1:dim(ovdf)[1] %in% friendsRow2, "N_both"])
      ovdf <- ovdf[1:dim(ovdf)[1] %ni% friendsRow2,]
      barcode_combine <- c(barcode_combine, friends2)
    }
    
    # Make a drop barcode and save our progress
    dropBarcode <- paste0(name, "_BC", formatC(idx, width=guess, flag="0", digits = 20), "_N", sprintf("%02d", length(barcode_combine)))
    
    # Annotate with new values
    nBC_keep[nBC_keep$Barcode %in% barcode_combine, "DropBarcode"] <- dropBarcode
    nBC_keep[nBC_keep$Barcode %in% barcode_combine, "OverlapReads"] <- OverlapReads
    
    idx <- idx + 1
    
    # Remove barcodes that we've dealt with
    nBC <- nBC[nBC$Barcode %ni% barcode_combine,]
  }
  
  # Get summary statistics per bead barcode
  merged_summary_statistics <- nBC_keep %>% group_by(DropBarcode) %>%
    summarise(uniqueNuclearFrags = (sum(UniqueNuclear) - (sum(as.numeric(OverlapReads))/2)), 
              totalNuclearFrags = sum(TotalNuclear),
              totalMitoFrags = sum(TotalMito),
              totalNCfilteredFrags = sum(TotalNC)) %>% data.frame()
  
  merged_summary_statistics$ncount <- as.numeric(substrRight(merged_summary_statistics$DropBarcode, 2))
  
  # PLUG and play here... function could return a vector of the names 
  keep_drop_barcodes <- as.character(merged_summary_statistics$DropBarcode)
  # keep_drop_barcodes <- knee_call(merged_summary_statistics$DropBarcode, merged_summary_statistics$uniqueNuclearFrags)
  
  # Filter out cells that pass adaptive knee filter and have less than 6 unique bead barcoes
  merged_summary_statistics_filt <- merged_summary_statistics %>% filter(ncount <= 6 & DropBarcode %in% keep_drop_barcodes)
  
  # Return a single number of the merge rate-- could slightly change this metric
  return(merged_summary_statistics_filt %>% summarize(Percent = 100*sum(ncount > 1)/n()) %>% as.numeric())
  
}

# Determine the merge proportion over a few values of jaccard_index
plot_df <- data.frame(
  jaccard_index = c(seq(0.001, 0.01, 0.002),
                    0.015, 0.02, 0.02)
)
plot_df$merge_proportion <- sapply(plot_df$jaccard_index, determineMergeProportion)

# Plot the function
ggplot(plot_df, aes(x = jaccard_index, y = merge_proportion)) + geom_point() + geom_line() +
  pretty_plot() + L_border() + labs(x = "Jaccard Fragment Index", y = "% of drops that were merged")





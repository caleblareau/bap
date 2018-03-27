library(data.table)
library(dplyr)
library(ggplot2)
library(tools)

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

rdsDir <- args[i+1]
tblOut <- args[i+2]

# For devel only
if(FALSE){
  rdsDir <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/frag_overlap"
  tblOut <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.implicatedBarcodes.csv"
}
rdsFiles <- list.files(rdsDir, full.names = TRUE)

lapply(rdsFiles, readRDS) %>%
  rbindlist() %>% as.data.frame() %>%
  group_by(barc1, barc2) %>%
  summarise(N_both = sum(n_both)/2, N_barc1 = sum(n_barc1), N_barc2 = sum(n_barc2)) %>% # overlaps called twice
  filter(N_both > 1) %>%
  mutate(jaccard_frag = round((N_both)/(N_barc1 + N_barc2 - N_both + 1),4)) %>%
  arrange(desc(jaccard_frag)) -> ovdf

write.table(ovdf, file = tblOut,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = ",")

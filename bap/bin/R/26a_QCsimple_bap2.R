options(warn=-1)

suppressMessages(suppressWarnings(library(GenomicRanges)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(tools)))

"%ni%" <- Negate("%in%")

options(scipen=999)

#-----------------
# Parse arguments
#-----------------

args <- commandArgs(trailingOnly = TRUE)
if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}

tempQcIn <- args[i+1]
mitoIn <- args[i+2]
firstQCout <- args[i+3]
species_mix <- args[i+4]

# For devel only
if(FALSE){
  tempQcIn <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/temp/jaccardPairsForIGV.basicQC-temp.tsv"
  mitoIn <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/temp/filt_split/jaccardPairsForIGV.chrM_frag.sumstats.tsv"
  firstQCout <- "/Users/clareau/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/final/jaccardPairsForIGV.basicQC.tsv"
  species_mix <- "none"
}

nuc_in <- fread(tempQcIn, col.names = c("cell_barcode", "n_total", "n_unique", "chromosome"))
mito_in <- fread(mitoIn, col.names = c("cell_barcode", "n_total", "n_unique", "chromosome"))

# Summarize counts
# No species mixing
nuc_sum <- nuc_in %>% group_by(cell_barcode) %>%
  summarize(totalNuclearFrags = sum(n_total), 
            uniqueNuclearFrags = sum(n_unique))

mito_sum <- mito_in %>% group_by(cell_barcode) %>%
  summarize(totalMitoFrags = sum(n_total), 
            uniqueMitoFrags = sum(n_unique))
basic_qc <- left_join(nuc_sum, mito_sum, by = "cell_barcode")

if(species_mix != "none"){
  
  nuc_in$mouse <- as.numeric(substr(nuc_in$chromosome,1,1) == "m")
  nuc_in$human <- as.numeric(substr(nuc_in$chromosome,1,1) == "h")
  
  human_df <- nuc_in %>% filter(human == 1) %>%
    group_by(cell_barcode) %>%
    summarize(totalHumanFrags = sum(n_total), 
              uniqueHumanFrags = sum(n_unique))
  
  mouse_df <- nuc_in %>% filter(mouse == 1) %>%
    group_by(cell_barcode) %>%
    summarize(totalMouseFrags = sum(n_total), 
              uniqueMouseFrags = sum(n_unique))
  QCstats <- left_join(basic_qc, human_df, by = "cell_barcode") %>%
    left_join(mouse_df, by = "cell_barcode")
  
} else {
  QCstats <- basic_qc
}

QCstats[is.na(QCstats)] <- 0
write.table(QCstats,
            file = firstQCout,
            quote = FALSE, row.names = FALSE, col.names = TRUE, sep = "\t")

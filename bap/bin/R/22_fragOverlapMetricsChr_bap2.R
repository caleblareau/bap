options(warn=-1)

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(GenomicRanges)))

"%ni%" <- Negate("%in%")

# This function / script for a given annotated chromosome fragments file saves a .rds file
# for each pair of barcodes given how many times they share an identical fragment in this 
# given chromosome


args <- commandArgs(trailingOnly = TRUE)

if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}
nc_threshold <- as.numeric(args[i+1])
hq_beads_file <- args[i+2]
anno_bedpe_file <- args[i+3]
out_frag_rds_file <- args[i+4]
out_nc_count_file <- args[i+5]
regularize_threshold <- as.numeric(args[i+6])

# Don't execute-- strictly for testing
if(FALSE){
  nc_threshold <- 6
  hq_beads_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/final/jaccardPairsForIGV.HQbeads.tsv"
  anno_bedpe_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/temp/filt_split/jaccardPairsForIGV.chr12.frag.bedpe.annotated.tsv.gz"
  out_frag_rds_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/temp/frag_overlap/jaccardPairsForIGV.chr12_overlapCount.rds"
  out_nc_count_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/temp/frag_overlap/jaccardPairsForIGV.chr12_ncCount.tsv"
}


# Import fragments
frags <- fread(cmd = paste0("zcat < ", anno_bedpe_file), col.names = c("read_name", "chr", "start", "end", "bead_barcode"))
HQbeads <- fread(hq_beads_file, col.names = "beads", header = FALSE)[[1]]

# Filter 1 for eligible barcodes
frags_filt1 <- frags[bead_barcode %in% HQbeads] 

# Quantify NC + export
frags_filt1[, `:=` (n_distinct_barcodes = .N), by = list(start, end)]
write.table(frags_filt1 %>% group_by(n_distinct_barcodes) %>% summarise(count = n()), 
            file = out_nc_count_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Filter out high NC values + first pass PCR duplicate removal at the bead level
nc_value <- frags_filt1[n_distinct_barcodes <= nc_threshold] 
frags_filt2 <- nc_value[, .(PCRdupCount = .N), by = list(chr, start, end, bead_barcode)]

# Determine baseline numbers
barcodes <- frags_filt2[["bead_barcode"]]
XBtable <- table(barcodes)
whichKeep <- names(XBtable)
nKept <- as.numeric(XBtable); names(nKept) <- whichKeep
rm(XBtable)

# Pull out barcode for retained fragments
# CONSIDER: padding by 1bp to fix old behavior
make_insert_overlap <- function(element_in_table){
  
  # Find find overlaps of Tn5 insertions
  inserts_df <- data.table(
    chr = c(frags_filt2[["chr"]]),
    start = c(frags_filt2[[element_in_table]]),
    end = c(frags_filt2[[element_in_table]]) ,
    bead_barcode = c(frags_filt2[["bead_barcode"]])
  )
  GAfilt <- makeGRangesFromDataFrame(inserts_df)
  ov <- findOverlaps(GAfilt, GAfilt, type = "equal")

  # Make a dataframe of all combinations that have fragments overlapping
  bc1 = barcodes[queryHits(ov)[ queryHits(ov) !=  subjectHits(ov)]]
  bc2 = barcodes[subjectHits(ov)[ queryHits(ov) !=  subjectHits(ov)]]
  boo_barc <- bc1 != bc2  # filter out reads that map to the same barcode
  bc1 <- bc1[boo_barc]
  bc2 <- bc2[boo_barc]
  rm(ov)
  
  hugeDF <- data.table(
    barc1 = ifelse(bc1 > bc2, bc1, bc2),
    barc2 = ifelse(bc1 > bc2, bc2, bc1)
  )
  hugeDF
}

#Double to consider left and right inserts
hugeDF <- rbind(
  make_insert_overlap("start"),
  make_insert_overlap("end")
)

# Break up previous massive dplyr call for speed in data.table
implicatedPairs <- hugeDF[, .(n_both = .N/2), by = list(barc1, barc2)]
implicatedPairs <- implicatedPairs[n_both >= regularize_threshold]
  
# Export
saveRDS(implicatedPairs, file = out_frag_rds_file)


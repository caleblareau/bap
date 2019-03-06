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

# Don't execute-- strictly for testing
if(FALSE){
  nc_threshold <- 6
  hq_beads_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/jaccardPairsForIGV.HQbeads.tsv"
  anno_bedpe_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/filt_split/jaccardPairsForIGV.chr12.frag.bedpe.annotated.tsv.gz"
  out_frag_rds_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/frag_overlap/jaccardPairsForIGV.chr12_overlapCount.rds"
  out_nc_count_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/frag_overlap/jaccardPairsForIGV.chr12_ncCount.tsv"
}

if(FALSE){
  nc_threshold <- 6
  base <- "/data/aryee/caleb/biorad/mouse_brain/N729_Exp110_sample8_combined_S1_2b2a2p/temp/filt_split/"
  anno_bedpe_file <- paste0(base, "/", "N729_Exp110_sample8_combined_S1.chr1.frag.bedpe.annotated.tsv.gz")
  hq_beads_file <- paste0("/data/aryee/caleb/biorad/mouse_brain/N729_Exp110_sample8_combined_S1_2b2a2p/final/N729_Exp110_sample8_combined_S1.HQbeads.tsv")
  out_frag_rds_file <- "/data/aryee/caleb/biorad/mouse_brain/N729_Exp110_sample8_combined_S1_2b2a2p/temp/filt_split/N729_Exp110_sample8_combined_S1.chr1_overlapCount.rds"
  out_nc_count_file <- "/data/aryee/caleb/biorad/mouse_brain/N729_Exp110_sample8_combined_S1_2b2a2p/temp/filt_split/N729_Exp110_sample8_combined_S1.chr1_ncCount.tsv"
  
}


# Import fragments
frags <- fread(cmd = paste0("zcat < ", anno_bedpe_file), col.names = c("read_name", "chr", "start", "end", "bead_barcode"))
HQbeads <- fread(hq_beads_file, col.names = "beads")[[1]]

# Filter 1 for eligible barcodes
frags_filt1 <- frags[bead_barcode %in% HQbeads] 

# Quantify NC + export
frags_filt1[, `:=` (n_distinct_barcodes = .N), by = list(start, end)]
write.table(frags_filt1 %>% group_by(n_distinct_barcodes) %>% summarise(count = n()), 
            file = out_nc_count_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Filter out high NC values + first pass PCR duplicate removal at the bead level
nc_value <- frags_filt1[n_distinct_barcodes <= nc_threshold] 
frags_filt2 <- nc_value[, .(PCRdupCount = .N), by = list(chr, start, end, bead_barcode)]

# Pull out barcode for retained
barcodes <- frags_filt2[["bead_barcode"]]
GAfilt <- makeGRangesFromDataFrame(frags_filt2)

# Find exact fragments
ov <- findOverlaps(GRanges(GAfilt), GRanges(GAfilt), type = "equal")

# Determine baseline numbers
XBtable <- table(barcodes)
whichKeep <- names(XBtable)
nKept <- as.numeric(XBtable); names(nKept) <- whichKeep
rm(XBtable)

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

# Break up previous massive dplyr call for speed in data.table
implicatedPairs <- hugeDF[, .(n_both = .N/2), by = list(barc1, barc2)]
implicatedPairs$n_barc1 <- nKept[implicatedPairs[["barc1"]]]
implicatedPairs$n_barc2 <- nKept[implicatedPairs[["barc2"]]]

# Export
saveRDS(implicatedPairs, file = out_frag_rds_file)



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
bin_resolution <- as.numeric(args[i+1])
nc_threshold <- as.numeric(args[i+2])
hq_beads_file <- args[i+3]
anno_bedpe_file <- args[i+4]
out_frag_rds_file <- args[i+5]
out_nc_count_file <- args[i+6]

# Don't execute-- strictly for testing
if(FALSE){
  bin_resolution <- 5000
  nc_threshold <- 1000
  hq_beads_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/final/jaccardPairsForIGV.HQbeads.tsv"
  anno_bedpe_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/temp/filt_split/jaccardPairsForIGV.chr12.frag.bedpe.annotated.tsv.gz"
  out_frag_rds_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/temp/frag_overlap/jaccardPairsForIGV.chr12_overlapCount.rds"
  out_nc_count_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap2/temp/frag_overlap/jaccardPairsForIGV.chr12_ncCount.tsv"
}


# Import fragments
frags <- fread(cmd = paste0("zcat < ", anno_bedpe_file), col.names = c("read_name", "chr", "start", "end", "bead_barcode"))
HQbeads <- fread(hq_beads_file, col.names = "beads")[[1]]

# Modify frag format for binning
frags$mid <- round((frags$start + frags$end)/2)
frags$start <- round(frags$mid/bin_resolution)*bin_resolution
frags$end <- frags$start + bin_resolution

# Filter 1 for eligible barcodes
frags_filt1 <- frags[bead_barcode %in% HQbeads] 

# Quantify NC + export
frags_filt1[, `:=` (n_distinct_barcodes = .N), by = list(start, end)]
write.table(frags_filt1 %>% group_by(n_distinct_barcodes) %>% summarise(count = n()), 
            file = out_nc_count_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

# Filter out high NC values + first pass PCR duplicate removal at the bead level
nc_value <- frags_filt1[n_distinct_barcodes <= nc_threshold] 
frags_filt2 <- nc_value[, .(PCRdupCount = .N), by = list(chr, start, end, bead_barcode)]

barcodes <- frags_filt2[["bead_barcode"]]
GAfilt <- makeGRangesFromDataFrame(frags_filt2)

# Find find overlaps of Tn5 insertions
ov <- findOverlaps(GAfilt, GAfilt, type = "equal")

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


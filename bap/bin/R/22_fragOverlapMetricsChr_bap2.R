options(warn=-1)

suppressMessages(suppressWarnings(library(Rsamtools)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))

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


# Import fragments
frags <- fread(cmd = paste0("zcat < ", anno_bedpe_file), col.names = c("chr", "start", "end", "read_name", "bead_barcode"))
HQbeads <- fread(hq_beads_file, col.names = "beads")[[1]]

# Filter 1 for eligible barcodes
frags_filt1 <- frags %>% filter(bead_barcode %in% HQbeads)
rm(frags)

# Quantify NC + export
nc_value <- frags_filt1 %>%
  group_by(chr, start, end) %>% mutate(n_distinct_barcodes = n_distinct(bead_barcode))

write.table(nc_value %>% group_by(n_distinct_barcodes) %>% summarise(count = n()), 
            file = out_nc_count_file, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
rm(frags_filt1)

# Filter out high NC values + first pass PCR duplicate removal at the bead level.
frags_filt2 <- nc_value %>%
  filter(n_distinct_barcodes <= nc_threshold) %>%
  group_by(chr, start, end, bead_barcode) %>% summarize(PCRdupCount = n()); rm(nc_value)

# Pull out barcode for retained
barcodes <- frags_filt2 %>% pull(bead_barcode) %>% as.character()
GAfilt <- makeGRangesFromDataFrame(frags_filt2)

# Find exact fragments
ov <- findOverlaps(GRanges(GAfilt), GRanges(GAfilt), type = "equal")

# Determine baseline numbers
XBtable <- table(barcodes)
whichKeep <- names(XBtable)
nKept <- as.numeric(XBtable); names(nKept) <- whichKeep
rm(XBtable)

# Make a dataframe of all combinations that have fragments overlapping
hugeDF <-
  data.frame(
    bc1 = barcodes[queryHits(ov)[ queryHits(ov) !=  subjectHits(ov)]],
    bc2 = barcodes[subjectHits(ov)[ queryHits(ov) !=  subjectHits(ov)]],
    stringsAsFactors = FALSE
  )
rm(ov)

# One dplyr command to save the day -- determine occurences of overlapping reads
hugeDF %>% filter(bc1 != bc2) %>%  # filter out reads that map to the same barcode
  mutate(barc1 = ifelse(bc1 > bc2, bc1, bc2),
         barc2 = ifelse(bc1 > bc2, bc2, bc1)) %>%  # switch to respect alphabetical order if necessary
  select(-one_of(c("bc1", "bc2"))) %>% # drop non-ordered barcode columns
  group_by(barc1, barc2) %>% summarise(n_both = n()/2) %>%  # group and count; divide by 2 to solve double counting (symmetry in overlap)
  mutate(n_barc1 = nKept[barc1], n_barc2 = nKept[barc2]) %>% as.data.frame() -> implicatedPairs # keep baseline numbers
rm(hugeDF)
rm(barcodes)

# Export
saveRDS(implicatedPairs, file = out_frag_rds_file)



options(warn=-1)
options(datatable.fread.input.cmd.message=FALSE)

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(GenomicRanges)))

"%ni%" <- Negate("%in%")

options(scipen=999)

# This script processes fragments and annotates them with the bead barcode ID
# based on a separate dictionary file

substrRight <- function(x, n = 6){
  substr(x, nchar(x)-n+1, nchar(x))
}

# I/O
args <- commandArgs(trailingOnly = FALSE)
nn <- length(args)

# Import parameters using logic from the end
blacklist_file <- args[nn-4]
frag_bedpe_file <- args[nn-3]
read_bead_file <- args[nn-2] 
annotated_out_file <- args[nn-1] 
unique_count_out_file <- args[nn]

# For devel only
if(FALSE){
  blacklist_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/bap/anno/blacklist/hg19-mm10.full.blacklist.bed"
  frag_bedpe_file <-"~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/filt_split/small_mix.hg19_chr1.frag.bedpe.gz"
  read_bead_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/filt_split/small_mix.hg19_chr1.read_bead.tsv.gz"
  annotated_out_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/filt_split/small_mix.hg19_chr1.frag.bedpe.annotated.tsv"
  unique_count_out_file <- "~/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/temp/filt_split/small_mix.hg19_chr1.bead_counts.tsv"
}

if(FALSE){
  
  base <- "/data/aryee/caleb/biorad/mouse_brain/N729_Exp110_sample8_combined_S1_2b2a2p/temp/filt_split/"
  blacklist_file <- "/data/aryee/caleb/pythondev/bap/bap/anno/blacklist/mm10.full.blacklist.bed"
  frag_bedpe_file <- paste0(base, "/", "N729_Exp110_sample8_combined_S1.chrX.frag.bedpe.gz")
  read_bead_file <-  paste0(base, "/", "N729_Exp110_sample8_combined_S1.chrX.read_bead.tsv.gz")
  annotated_out_file <- paste0(base, "/", "N729_Exp110_sample8_combined_S1.chrX.frag.bedpe.annotated.tsv.gz")
  unique_count_out_file <- paste0(base, "/", "N729_Exp110_sample8_combined_S1.chrX.bead_counts.tsv")
  
}

# Import frags and annotate with bead
frags <- fread(cmd = paste0("zcat < ", frag_bedpe_file), col.names = c("chr", "start", "end", "read_name")) 
bead_read <- fread(cmd = paste0("zcat < ", read_bead_file), col.names = c("read_name", "bead_id"))%>% na.omit() %>% unique()
mdf <- merge(frags, bead_read, by= "read_name") %>% na.omit()

# Filter for fragments overlapping the blacklist
bl <- fread(blacklist_file, col.names = c("chr", "start", "end")) %>% data.frame() %>% makeGRangesFromDataFrame()
ov_bl <- findOverlaps(bl, makeGRangesFromDataFrame(mdf))
blacklist_reads <- mdf$read_name[subjectHits(ov_bl)]
mdf <- mdf[read_name %ni% blacklist_reads]

# NOTE: the majority of the missing read names are instances where the MAPQ filter was not surpassed.

# Quantify the number of unique fragments per barcode
pcr_dup_df <- mdf[, .(count = .N), by = list(chr, start, end, bead_id)]
out_bead_quant <- pcr_dup_df[,.(nUnique = .N), by = bead_id]

# Export tables
write.table(mdf, file = annotated_out_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(out_bead_quant, file = unique_count_out_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



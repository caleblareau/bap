options(warn=-1)
options(datatable.fread.input.cmd.message=FALSE)

suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(GenomicRanges)))

"%ni%" <- Negate("%in%")

options(scipen=999)

# This script takes fragments and cleans them up (PCR deduplication; sort by coordinate)

substrRight <- function(x, n = 6){
  substr(x, nchar(x)-n+1, nchar(x))
}

# I/O
args <- commandArgs(trailingOnly = FALSE)
nn <- length(args)

# Import parameters using logic from the end
blacklist_file <- args[nn-2]
frag_bedpe_file <- args[nn-1]
out_file <- args[nn] 

# Import frags and annotate with bead
frags <- fread(cmd = paste0("cat < ", frag_bedpe_file), col.names = c("chr", "start", "end", "read_name")) 

# Filter for fragments overlapping the blacklist
bl <- fread(paste0("cat < ", blacklist_file), col.names = c("chr", "start", "end")) %>% data.frame() %>% makeGRangesFromDataFrame()
ov_bl <- findOverlaps(bl, makeGRangesFromDataFrame(frags))
blacklist_reads <- frags$read_name[subjectHits(ov_bl)]
frags <- frags[read_name %ni% blacklist_reads]

# NOTE: the majority of the missing read names are instances where the MAPQ filter was not surpassed.

# Perform PCR deduplication
pcr_dedup_df <- frags[, .(count = .N), by = list(chr, start, end)]
pcr_dedup_df <- pcr_dedup_df[order(start)]

# Export tables
write.table(pcr_dedup_df %>% dplyr::select(c("chr", "start", "end")),
            file = out_file, row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")



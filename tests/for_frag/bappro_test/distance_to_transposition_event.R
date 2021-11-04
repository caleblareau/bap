library(data.table)
library(GenomicRanges)
library(dplyr)

#Import insertion events
bap_pairs <- fread("jaccardPairsForIGV.implicatedBarcodes.csv") %>% 
  mutate(bc_pair = ifelse(barc1 < barc2, paste0(barc1, "_", barc2), paste0(barc2, "_", barc1))) %>%
  filter(merged) %>% pull(bc_pair)
lapply(list.files(".", pattern = ".gz"),fread) %>% rbindlist() -> d

data.frame(
  chr = d$V2,
  start = d$V3,
  end = d$V3 +1, 
  barcode = d$V5
) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE) -> gr1

data.frame(
  chr = d$V2,
  start = d$V4,
  end = d$V4 +1, 
  barcode = d$V5
) %>% makeGRangesFromDataFrame(keep.extra.columns = TRUE) -> gr2

# Get distances to close transposition events and summarize things
grd <- distanceToNearest(gr1, gr2)

dfdf <- data.frame(
  distance = grd@elementMetadata@listData$distance,
  idx1 = gr1@ranges@start[grd@from],
  idx2 = gr2@ranges@start[grd@to],
  bc1 = gr1@elementMetadata$barcode[grd@from],
  bc2 = gr2@elementMetadata$barcode[grd@to]
) %>% mutate(bc_pair = ifelse(bc1 < bc2, paste0(bc1, "_", bc2), paste0(bc2, "_", bc1))) %>% 
  mutate(dist_CL = idx2 - idx1) %>%
  mutate(bc_pair_anno = case_when(
    bc_pair %in% bap_pairs ~ "bap_pair", 
    bc1 == bc2 ~ "identical_bead_barcode", 
    TRUE ~ "other"
  )) %>%
  mutate(bc_pair_anno_detail = case_when(
    bc_pair %in% bap_pairs ~ bc_pair, 
    bc1 == bc2 ~ "z_identical_bead_barcode", 
    TRUE ~ "z_other"
  ))

df <- dfdf %>% 
  filter(dist_CL < 31 & dist_CL > -31) %>%
  group_by(
    dist_CL, bc_pair_anno
  ) %>% summarize(count = n())

# Highlight a few important boys
df <- df %>% mutate(color = case_when(dist_CL == "1" ~ "1", 
                                      dist_CL == "7" ~ "7", 
                                      dist_CL == "9" ~ "9", 
                                      TRUE ~ "other"))

ggplot(df, aes(x = dist_CL, y = count, fill = color)) +
  geom_bar(stat = "identity", color = "black") + 
  facet_wrap(~bc_pair_anno) +
  scale_fill_manual(values = c("green4", "dodgerblue2", "firebrick", "lightgrey")) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Distance between transposition events", y = "count", fill = "highlights")

# Characterize individual pairs
df2 <- dfdf %>% 
  filter(dist_CL < 31 & dist_CL > -31) %>%
  group_by(
    dist_CL, bc_pair_anno_detail
  ) %>% summarize(count = n())
df2 <- df2 %>% mutate(color = case_when(dist_CL == "1" ~ "1", 
                                      dist_CL == "7" ~ "7", 
                                      dist_CL == "9" ~ "9", 
                                      TRUE ~ "other"))

ggplot(df2, aes(x = dist_CL, y = count, fill = color)) +
  geom_bar(stat = "identity", color = "black") + 
  facet_wrap(~bc_pair_anno_detail, scales = "free_y") +
  scale_fill_manual(values = c("green4", "dodgerblue2", "firebrick", "lightgrey")) +
  theme_bw() + theme(legend.position = "bottom") +
  labs(x = "Distance between transposition events", y = "count", fill = "highlights")

options(warn=-1)

suppressMessages(suppressWarnings(library(ggplot2)))
suppressMessages(suppressWarnings(library(dplyr)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(tools)))

"%ni%" <- Negate("%in%")

# Helper functions
get_ticks <- function(data, col_name) {
  n1 <- floor(log10(range(data[col_name])))
  if (n1[1] == -Inf) n1[1] = 0
  pow <- seq(n1[1], n1[2]+1)
  ticks <- as.vector(sapply(pow, function(p) (c(1,5)*10^p)))
  return(ticks)
}

args <- commandArgs(trailingOnly = TRUE)

if(file_path_sans_ext(basename(args[1])) == "R"){
  i <- 2
} else { # Rscript
  i <- 0
}
bapParamsFile <- args[i+1]
beadBarcodesFile <- args[i+2]
jaccardFragsFile <- args[i+3]

# Don't execute-- strictly for testing
if(FALSE){
  bapParamsFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/knee/test.small.bapParams.csv"
  beadBarcodesFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/knee/test.small.barcodeQuantSimple.csv"
  jaccardFragsFile <- "/Volumes/dat/Research/BuenrostroResearch/lareau_dev/bap/tests/bap_out/final/test.small.implicatedBarcodes.csv.gz"
  
}

# Make output files
beadBarcodesKneePlot <- gsub(".barcodeQuantSimple.csv", ".beadBarcodeKnee.pdf", beadBarcodesFile)
jaccardKneePlot <- gsub(".barcodeQuantSimple.csv", ".jaccardOverlapKnee.pdf", beadBarcodesFile)

stopifnot(file.exists(bapParamsFile))

# Parse knee-called values
bapParams <- read.table(bapParamsFile, sep = ",")
bead_threshold <- bapParams %>% filter(V1 == "bead_threshold") %>% pull(V2) %>% as.character() %>% as.numeric()
jaccard_threshold <- bapParams %>% filter(V1 == "jaccard_threshold") %>% pull(V2)%>% as.character() %>% as.numeric()
bead_threshold <- ifelse(is.na(bead_threshold), 0, bead_threshold)

# Make bead barcodes knee plot
if(file.exists(beadBarcodesFile)){
  
  # Import bead counts
  bead_vec <- fread(beadBarcodesFile)[["V2"]]
  
  # Make a data frame for plotting with ticks
  barcode_counts_mut <- data.frame(count = bead_vec) %>%
    arrange(desc(count)) %>%
    mutate(rank = 1:n(), passKnee = count > bead_threshold)
  ticks_at_y <- get_ticks(barcode_counts_mut, 'count')
  ticks_at_x <- get_ticks(barcode_counts_mut, 'rank')
  x_intercept_bead <- sum(barcode_counts_mut$passKnee)
  
  barcode_counts_plot <- barcode_counts_mut %>%
    ggplot(aes(x=rank, y=count)) +
    scale_y_log10(breaks = ticks_at_y, labels = as.integer(ticks_at_y)) +
    scale_x_log10(breaks = ticks_at_x, labels = as.integer(ticks_at_x)) +
    xlab("Barcode in rank-descending order") +
    ylab("Reads per barcode") +
    theme(axis.text.x = element_text(angle=90)) +
    geom_point(aes(color = passKnee)) +
    scale_color_manual(values = c("black", "dodgerblue3")) +
    theme(legend.position = c(0.02, 0.1)) +
    labs(color = "Pass Knee") +
    geom_vline(xintercept = x_intercept_bead, color = "dodgerblue4") +
    geom_hline(yintercept = bead_threshold, color = "dodgerblue4")
  
  # Export to PDF for vectorized stuff
  ggsave(barcode_counts_plot, file = beadBarcodesKneePlot,
         width = 6, height = 6, useDingbats=FALSE)
  
  # Also save as a png for good measure
  ggsave(barcode_counts_plot, file = gsub(".pdf$", ".png", beadBarcodesKneePlot),
         width = 6, height = 6)
  
  #--------------------------
  # Additional QC plots
  #--------------------------
  
  big_c <- ceiling(max(log10(barcode_counts_mut$count)))
  
  # Density Plot
  barcode_counts_mut %>%
    filter(count > 50) %>% 
    ggplot() +
    geom_density(aes(x = log10(count+1))) +
    ylab("Density") +
    scale_x_continuous(limits = c(0,big_c)) +
    xlab("Count per barcode log10") +
    geom_vline(xintercept = log10(x_intercept_bead), colour= "dodgerblue4") ->
    density_plot
  
  beadBarcodesKneeDensityPlot <- gsub(".barcodeQuantSimple.csv", ".beadBarcodeKneeDensity.pdf", beadBarcodesFile)
  ggsave(density_plot, file = beadBarcodesKneeDensityPlot, useDingbats=FALSE, 
         width = 6, height = 6)
  
  # Knee Plot
  barcode_counts_mut %>%
    mutate(cum_sum = cumsum(count)) %>%
    ggplot() +
    geom_line(aes(x=rank, y=cum_sum)) +
    ylab("Cumulative read count") +
    xlab("Barcode rank") +
    geom_vline(xintercept = x_intercept_bead, colour= "dodgerblue4") ->
    knee_plot
  
  beadBarcodesKneeCurvePlot <- gsub(".barcodeQuantSimple.csv", ".beadBarcodeKneeCurve.pdf", beadBarcodesFile)
  ggsave(knee_plot, file = beadBarcodesKneeCurvePlot,
         width = 6, height = 6,  useDingbats=FALSE)
  
}

# Make jaccard overlaps knee plot
if(file.exists(jaccardFragsFile)){
  
  # Import jaccard overlap measures
  jaccard_frag_vec <- fread(cmd = paste0("zcat < ", jaccardFragsFile, " | head -1000000"))[["jaccard_frag"]]
  
  # Make a data frame for plotting with ticks
  jaccard_mut <- data.frame(jaccard_frag = jaccard_frag_vec) %>%
    arrange(desc(jaccard_frag)) %>%
    mutate(rank = 1:n(), passKnee = jaccard_frag > jaccard_threshold)
  ticks_at_y <- get_ticks(jaccard_mut, 'jaccard_frag')
  ticks_at_x <- get_ticks(jaccard_mut, 'rank')
  x_intercept_frag <- sum(jaccard_mut$passKnee)
  
  jaccard_plot <- jaccard_mut %>%
    ggplot(aes(x=rank, y=jaccard_frag)) +
    scale_y_log10(breaks = ticks_at_y, labels = as.numeric(ticks_at_y)) +
    scale_x_log10(breaks = ticks_at_x, labels = as.integer(ticks_at_x)) +
    xlab("bap overlap score in rank-descending order") +
    ylab("bap overlap score per barcode pair") +
    theme(axis.text.x = element_text(angle=90)) +
    geom_point(aes(color = passKnee)) +
    scale_color_manual(values = c("black", "dodgerblue3")) +
    theme(legend.position = c(0.02, 0.1)) +
    labs(color = "Pass Knee") +
    geom_vline(xintercept = x_intercept_frag, color = "dodgerblue4") +
    geom_hline(yintercept = jaccard_threshold, color = "dodgerblue4")
  
  # Export to PDF for vectorized stuff
  ggsave(jaccard_plot, file = jaccardKneePlot,
         width = 6, height = 6, useDingbats=FALSE)
  
  # Also save as a png for good measure
  ggsave(jaccard_plot, file = gsub(".pdf$", ".png", jaccardKneePlot),
         width = 6, height = 6)
  
}
statusUpdateTxt <- gsub(".barcodeQuantSimple.csv", ".kneesPlotted.txt", beadBarcodesFile)
write.table(data.frame("yes"), 
            file = statusUpdateTxt, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)

options(warn=-1)

suppressMessages(suppressWarnings(library(tools)))
suppressMessages(suppressWarnings(library(data.table)))
suppressMessages(suppressWarnings(library(dplyr)))

# Simple function to determine the mode of a vector
# See here: https://stackoverflow.com/questions/2547402/is-there-a-built-in-function-for-finding-the-mode
get_mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}


get_local_minima_CL <- function(x) {
  
  # Find the indicies where we go from negative diff's to positive
  # if this returns less then 0 then we are descending
  y <- diff(c(Inf, x)) < 0L
  
  # rle to get the number of trues (downwards steps) and falses (upwards)
  # cumsum will then the index at which each of these inflections happens
  y <- cumsum(rle(y)$lengths)
  
  # If we only get three values, then this becomes a problem-- 
  # the every/other TRUE,FALSE,TRUE will only keep the extreme values
  # and then remove the only true minimum
  if(length(y) > 3){
    y <- y[seq.int(1L, length(y), 2L)]
  }
  
  # Seth's modification for removing duplicated elements at the beginning
  if (x[[1]] == x[[2]]) {
    y <- y[-1]
  }
  y
}

# Complicated function converted by CL
get_density_threshold_CL <- function(count_vector=NULL, type, logTransform = TRUE) {
  
  # Initialize using some reasonable value and filter anything below
  threshold <- min(get_mode(count_vector)) * 0.001
  filtered_tbl <- data.frame(count = count_vector[count_vector > threshold]) 
  
  # Parameterize the log transformation to work with non-count data
  # May get confusing evtually but for now, seemingly a decent hack
  if(logTransform){
    filtered_tbl$log_counts <- log10(filtered_tbl$count)
  } else{
    filtered_tbl$log_counts <- filtered_tbl$count
  }
  
  # Calculate the density using a gaussian kernel
  xx_values = 10000
  vector_density <- density(filtered_tbl$log_counts, bw=0.1, kernel="gaussian", n=xx_values,
                            from=min(filtered_tbl$log_counts),to=max(filtered_tbl$log_counts))
  
  local_mins <- get_local_minima_CL(vector_density$y)
  
  # If a minima was called at the very start or end of the distribution, remove it
  local_mins <- local_mins[!(local_mins == 1 | local_mins == length(vector_density$y))]
  
  # Find the first local min that meets the following criteria, starting right to left
  local_min <- Find(function(x) {
    
    # Make sure that the selected min includes at least 20% of barcodes
    # and that the difference between the min and the max differences by some appreciable amount
    # both in terms of absolute difference (Caleb changed 0.5 to 0.05) AND relative difference
    abs_difference <- ifelse(logTransform, 0.5, 0.05)
    
    return (x >= (0.2 * xx_values) & ((max(filtered_tbl$log_counts) - vector_density$x[x]) > abs_difference |
                                        (vector_density$x[x] < max(filtered_tbl$log_counts) / 2))) 
    
  }, rev(local_mins))
  
  if (!is.null(local_min)) {
    if(logTransform){
      threshold <- 10^vector_density$x[local_min]
    } else {
      threshold <- vector_density$x[local_min]
    }
    message("Setting knee threshold to: ", threshold)
    
  } else {
    
    # Null local min-- take a best guess
    if(logTransform){
      message("No reliable knee found-- setting threshold to 0")
      threshold <- 0
    } else {
      message("No reliable knee found-- setting threshold to 0")
      threshold <- 0
    }
    local_min <- 1
    local_mins <- 1
  }
  
  safety <- 0
  
  # Safe guard for Jaccard Index failure
  if(type == "jaccard" & (threshold > 0.5 | threshold < 0.000001)){
    message("No reliable knee found-- setting threshold to 0.005")
    safety <- 0.005
  } 
  
  # Safe guard for knee counts failure
    if(type == "bead" & (threshold > 100000 | threshold < 100)){
    message("No reliable knee found-- setting threshold to 500")
    safety <- 500
  } 
  
  # Safety is with the guard rails; threshold is what the knee calls
  safety <- ifelse(safety > 0, safety, threshold)
  
  return(c(safety, threshold))
}


# See https://github.com/broadinstitute/picard/blob/master/src/main/java/picard/sam/DuplicationMetrics.java#L115
estimateLibrarySize <- function(nTotal, nUnique){
  
  f <- function(x, c, n) {
    return(c / x - 1 + exp(-n / x))
  }
  
  m = 1
  M = 100
  
  nDuplicates <- nTotal - nUnique
  
  # Checks
  stopifnot(nTotal > 0)
  stopifnot(nDuplicates >0)
  stopifnot(nUnique > 0)
  
  if (nUnique >= nTotal | f(m * nUnique, nUnique, nTotal) < 0) {
    stop("Error: invalid inputs")
  }
  
  while (f(M * nUnique, nUnique, nTotal) > 0) {
    M <- M*10.0
  }
  
  for(i in seq(0, 40)) {
    r <-  (m + M) / 2.0
    u <- f(r * nUnique, nUnique, nTotal);
    if (u == 0) {
      break
    } else if (u > 0) {
      m = r
    } else if (u < 0) {
      M = r
    }
  }
  
  return((nUnique * (m + M) / 2.0))
}

estimateLibrarySize(1000, 800)

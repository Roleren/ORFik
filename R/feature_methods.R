#' Subset GRanges to get desired frame. GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#'
#’ @export
#' @param x A tiled to size of 1 GRanges object
#' @param frame A numeric indicating which frame to extract
#' @return GRanges object reduced to only first frame
#' @examples
#' #subset_to_frame(x, 1)
#'
subset_to_frame <- function(x, frame){
  if(as.vector(strand(x) == "+")[1]){
    x[seq(frame, length(x), 3)]
  }else{
    x[seq(length(x) + 1 - frame, 1, -3)]
  }
}


#' Subset GRanges to get stop codons. GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#'
#' @param x A tiled to size of 1 GRanges object
#' @return GRanges object reduced to only stop codon
#' @export
#' @examples
#' #subset_to_stop(x)
#'
subset_to_stop <- function(x){
  if(as.vector(strand(x))[1] == "+"){
    x[c(length(x) - 3, length(x) - 4, length(x) - 5)]
  } else {
    x[c(4, 5, 6)]
  }
}


#' Subset GRanges to get coverage. GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#'
#’
#' @param cov A coverage object from coverage()
#' @param y GRanges object for which coverage should be extracted
#' @return numeric vector of coverage of input GRanges object
#' @export
#' @examples
#' #subset_coverage(x)
#'
subset_coverage <- function(cov, y) {
  cov1 <- cov[[as.vector(seqnames(y)[1])]]
  return(as.vector(cov1[ranges(y)]))
}


#' Calucalte Entropy value of input reads. Based on ...
#'
#' @export
#' @param countsOver A numeric vector containing reads over given ORF
#' @return numeric value of entropy
#' @examples
#' calculateEntropy(c(30,15,5,45,23,4,23,4,1,14,7,3,24,6,2,10,6,2))
#'
calculateEntropy <- function(countsOver) {
  reg_count <- 0
  reg_len <- 3
  N <- sum(countsOver)
  L <- length(countsOver)
  if(N <= L){
    reg_len <- floor(L/N)
  }
  reg_count <- L/reg_len
  Hx <- 0
  MHx <- 0
  for(h in 0:(reg_count-1)){
    which_reads <- (1 + h * reg_len):(h * reg_len + reg_len)
    Xi <- sum(countsOver[which_reads]) / N
    Hx <- Hx + ifelse(is.nan(Xi * log2(Xi)), 0, Xi * log2(Xi))
    MHx <- MHx + 1/reg_count * log2(1/reg_count)
  }
  return(Hx/MHx)
}


#' Calculate Coverage of the input reads.
#' What is percent of the vector covered by the values higher than zero?
#'
#' @export
#' @param countsOver A numerc vector
#' @return numeric value in between 0 and 1
#'
calculateCoverage <- function(countsOver) {
  return(sum(countsOver >= 1)/length(countsOver))
}

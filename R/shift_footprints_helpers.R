#' Shift ribo-seq reads using cigar string
#'
#' @param cigar the cigar of the reads
#' @param shift the shift as integer
#' @param is_plus_strand logical
#' @return the shifted read
#'
parseCigar <- function(cigar, shift, is_plus_strand) {
  c_signs <- unlist(explodeCigarOps(cigar))
  c_counts <- unlist(explodeCigarOpLengths(cigar))

  i <- ifelse(is_plus_strand, 0L, length(c_signs) + 1L)
  increment <- ifelse(is_plus_strand, 1L, -1L)
  limit <- 0L
  plusShift <- 0L
  while (shift >= limit) {
    i = i + increment
    if (c_signs[i] == "M") {
      limit = limit + c_counts[i]
    } else if (c_signs[i] == "N" || c_signs[i] == "D") {
      plusShift = plusShift + c_counts[i]
    } else if (c_signs[i] == "I") {
      plusShift = plusShift - c_counts[i]
    } else {
      warning(paste0("Not supported sign:", c_signs[i]))
    }
  }
  shift = shift + plusShift
  return(shift)
}

#' Find if there is periodicity in the vector
#'
#' Checks if there is a periodicity and if the periodicity is 3.
#'
#' It uses Fourier transform for finding periodic vectors
#' @param x (numeric) Vector of values to detect periodicity of 3 like in
#' RiboSeq data.
#' @return a logical, if it is periodic.
#' @importFrom stats fft spec.pgram
#'
isPeriodic <- function(x) {
  if (sum(x) == 0) return(FALSE)
  amplitudes <- abs(fft(x))
  amp <- amplitudes[2 : (length(amplitudes) / 2 + 1)]
  periods <- 1 / spec.pgram(x = x, plot = FALSE)$freq
  return((periods[which.max(amp)] > 2.9) & (periods[which.max(amp)] < 3.1))
}

#' Get the offset for specific RiboSeq read width
#' @param x a vector with count per position to analyse, assumes the zero
#'  is in the middle + 1 (position 0)
#' @param feature (character) either "start" or "stop"
#' @return a single numeric offset
#' @family pshifting
#'
changePointAnalysis <- function(x, feature = "start") {
  meta <- x[seq.int(40L)]
  pos <- -(length(x)/2):(length(x)/2 - 1)
  if (feature == "start") {
    means <- c()
    for (j in seq.int(15L, 35L)) {
      m <- mean(meta[seq.int(j, 40L)]) - mean(meta[seq.int(j - 1L)])
      means <- c(means, m)
    }
    shift <- which.max(abs(means)) + 14
    offset <- pos[shift]
  } else if (feature == "stop") {
    shift <- which.max(meta)
    offset <- pos[shift] + 6
  }
  return(offset)
}

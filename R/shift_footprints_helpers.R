#' Shift ribo-seq reads using cigar string
#'
#' Example if you want to change a read of length 20, by +12.
#' You need to account for gaps etc, this is done using the
#' cigar string of the read.
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

#' Find if there is a periodicity of 3 in the vector
#'
#' It uses Fourier transform for finding periodic vectors on the transcript
#' normalized counts over all CDS TIS regions from -30 to 29, where TIS is
#' 0.\cr
#' Checks if there is a periodicity and if the periodicity is 3,
#' more precisely between 2.9 and 3.1.
#'
#' Transcript normalized means per CDS TIS region, count reads per position,
#' divide that number per position by the total of that transcript, then sum
#' up these numbers per position for all transcripts.
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
#'
#' Creates sliding windows of transcript normalized counts per position
#' and check which window has most in upstream window vs downstream window.
#' Pick the position with highest absolute value maximum of the window difference.
#' Checks windows with split sites between positions -16 to +4, where 0 is TIS.
#'
#' Transcript normalized means per CDS TIS region, count reads per position,
#' divide that number per position by the total of that transcript, then sum
#' up these numbers per position for all transcripts.
#' @param x a vector with count per position to analyse, assumes the zero
#' position (TIS) is in the middle + 1 (position 0). Default it is size 60,
#' from -30 to 29 in p-shifting
#' @param feature (character) either "start" or "stop"
#' @param max.pos integer, default 40L, subset x to go from index 1 to max.pos,
#'  if tail is not relevant.
#' @param interval integer vector, default seq.int(15L, 35L).
#'  Seperation points for upstream and downstream windows
#' @return a single numeric offset, -12 would mean p-site is 12 bases upstream
#' @family pshifting
#'
changePointAnalysis <- function(x, feature = "start", max.pos = 40L,
                                interval = seq.int(15L, 35L)) {
  if (max.pos > length(x)) stop("Can not subset max.pos > length of x")
  if (!all(interval %in% seq.int(x)))
    stop("interval vector must be subset of x indices!")
  meta <- x[seq.int(max.pos)] # First 40 or other specified
  pos <- -(length(x)/2):(length(x)/2 - 1) # The positions
  if (feature == "start") {
    means <- c()
    # upstream window vs downstream window, Check counts in area: pos -16 to +4
    for (j in interval) {
      down <- meta[seq.int(j, max.pos)]
      only.per.third <- down[pos[seq.int(j, max.pos)] %% 3 == pos[j] %% 3]
      downstream <- mean(only.per.third) # down window
      upstream <- mean(meta[seq.int(j - 2L)]) # up window (-2 is skip 1 base, less noise)
      m <- downstream - upstream
      means <- c(means, m)
    }
    shift <- which.max(abs(means)) + min(interval) - 1 # +14 since we start on 15 in interval
    offset <- pos[shift]
  } else if (feature == "stop") {
    shift <- which.max(meta)
    offset <- pos[shift] + 6
  }
  return(offset)
}

#' Convert percentage to ratio of 1
#'
#' 50 -> 0.5 etc, if length cds > minimum.cds
#' @param top_tx numeric
#' @param cds GRangesList object
#' @param minimum.cds numeric, default 1000
#' @return numeric, as ratio of 1
percentage_to_ratio <- function(top_tx, cds, minimum.cds = 1000) {
  if (!is.numeric(top_tx) |top_tx < 0 | top_tx > 100)
    stop("top_tx must be numeric between 0 and 100 %")
  top_tx <- if (length(cds) > minimum.cds) { # as ratio of 1
    (100 - top_tx) / 100
  } else 0
  return(top_tx)
}

#' Pre shifting plot analysis
#'
#' For internal use only!
#' @param rw a data.table of position, score and fraction (read length)
#' of either TIS or TES (translation end site, around 3' UTR)
#' @param heatmap a logical or character string, default FALSE.
#' If TRUE, will plot heatmap of
#' raw reads before p-shifting to console, to see if shifts given make sense.
#' You can also set a filepath to save the file there.
#' @param region a character string, default "start of CDS"
#' @return invisible(NULL)
footprints.analysis <- function(rw, heatmap, region = "start of CDS") {
  region <- paste("Position relative to", region)
  if (heatmap != FALSE) {
    gg <- coverageHeatMap(rw, scoring = "transcriptNormalized",
                          xlab = region)
    plot(gg)
    if (is.character(heatmap)) {
      ggsave(heatmap, gg)
    }
  }
  return(invisible(NULL))
}

#' Calculate metaplot coverage of reads around input GRangesList object.
#'
#' Sums up coverage over set of GRanges objects that.
#' @param x GRanges object of your reads.
#' You should resize them beforehand to width of 1 to focus on
#' 5' ends of footprints.
#' @param windows GRanges object of your CDSs start or stop postions. Its width
#' has to be even number as we will assume in the middle is position zero which
#' is included in the downstream window.
#' @return A data.frame with average counts (avg_counts) of reads mapped to
#' positions (position) specified in windows along with frame (frame).
#' @export
#' @importFrom BiocGenerics Reduce
#' @examples
#' windows <- GenomicRanges::GRangesList(
#'   GenomicRanges::GRanges(seqnames = "chr1",
#'                          ranges = IRanges::IRanges(c(50, 100), c(80, 200)),
#'                          strand = "-"))
#' x <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges =  IRanges::IRanges(c(100, 180), c(200, 300)),
#'   strand = "-")
#' metaWindow(x, windows)
#'
metaWindow <- function(x, windows) {

  window_size <- unique(sum(width(windows)))
  if (length(window_size) != 1) {
    stop("All input 'windows' should have the same sum(width(windows))")
  }
  if (window_size %% 2 != 0) stop("Width of the window has to be even number.")
  window_size <- (window_size)/2
  cvg <- coverageByTranscript(x, windows)
  cvg <- Reduce(`+`, cvg) / length(cvg)

  hitMap <- data.frame(avg_counts = as.vector(cvg),
                       position = -window_size:(window_size - 1),
                       frame = c(rev(rep_len(3:1, window_size)),
                                 rep_len(seq(3), window_size)))
  return(hitMap)
}


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

  i = ifelse(is_plus_strand, 0L, length(c_signs) + 1L)
  increment = ifelse(is_plus_strand, 1L, -1L)
  limit = 0L
  plusShift = 0L
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


#' Get the transcripts that have minimum lengths of leaders and cds.
#'
#' Filter transcripts to those whohave 5' UTR, CDS, 3' UTR of
#' some lengths, pick the longest per gene.
#' @param txdb a TxDb object from gtf
#' @param minFiveUTR (integer) minimum bp for 5' UTR during filtering for the
#' transcripts
#' @param minCDS (integer) minimum bp for CDS during filtering for the
#' transcripts
#' @param minThreeUTR (integer) minimum bp for 3' UTR during filtering for the
#' transcripts
#' @return a character vector of valid tramscript names
#' @export
#' @examples
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
#' txNames <- txNamesWithLeaders(txdb)
#'
txNamesWithLeaders <- function(txdb, minFiveUTR = 30L,
                               minCDS = 150L, minThreeUTR = 30L) {
  if(!is(txdb, "TxDb")) stop("txdb must be a TxDb object")

  tx <- data.table::setDT(
    GenomicFeatures::transcriptLengths(
      txdb, with.cds_len = TRUE, with.utr5_len = TRUE, with.utr3_len = TRUE))
  tx <- tx[tx$utr5_len >= minFiveUTR & tx$cds_len >= minCDS &
             tx$utr3_len >= minThreeUTR, ]
  gene_id <- cds_len <- NULL
  data.table::setorder(tx, gene_id, -cds_len)
  tx <- tx[!duplicated(tx$gene_id), ]
  tx <- tx[!is.na(tx$gene_id)]
  return(tx$tx_name)
}


#' Get Start and Stop codon within specified windows over CDS.
#'
#' For each cds in `txdb` object, filtered by `txNames`,
#' get a window around start and stop codons within `window_size`
#' downstream and upstream of the codon.
#' @param txdb a txdb object of annotations
#' @param txNames a character vector of the transcript names to use
#' @param start (logical) whether to include start codons
#' @param stop (logical) whether to incude stop codons
#' @param window_size (integer) size of the window to extract upstream of the
#' start/stop codon and downstream
#' @param cds a GRangesList with cds, a speedup if you already have them
#'  loaded.
#' @return a list with two slots "starts" and "stops", each contains a
#' GRangesList of windows around start and stop codons for the transcripts of
#' interest
#' @export
#' @examples
#' gtf_file <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' txdb <- GenomicFeatures::makeTxDbFromGFF(gtf_file, format = "gtf")
#' txNames <- txNamesWithLeaders(txdb)
#' getStartStopWindows(txdb, txNames)
#'
getStartStopWindows <- function(
  txdb, txNames, start = TRUE, stop = TRUE, window_size = 30L, cds = NULL) {

  if (is.null(cds)) {
    cds <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)[txNames]
  } else {
    cds <- cds[txNames]
  }
  cdsTiled <- tile1(cds, matchNaming = FALSE)

  if (start) {
    fiveUTRs <- GenomicFeatures::fiveUTRsByTranscript(
      txdb, use.names=TRUE)[txNames]
    fiveUTRs <- removeMetaCols(fiveUTRs)
    cdsHead <- heads(cdsTiled, window_size)
    fiveTiled <- tile1(fiveUTRs, matchNaming = FALSE)
    fiveTail <- tails(fiveTiled, window_size)
    partitions <- c(fiveTail, cdsHead)
    merged <- unlist(partitions, use.names = FALSE)
    start <- reduceKeepAttr(split(merged,
      names(partitions)[groupings(partitions)]))
  } else {
    start <- NULL
  }

  if (stop) {
    threeUTRs <- GenomicFeatures::threeUTRsByTranscript(
      txdb, use.names=TRUE)[txNames]
    threeUTRs <- removeMetaCols(threeUTRs)
    cdsTail <- tails(cdsTiled, window_size)
    threeHead <- downstreamN(threeUTRs, window_size)
    partitions <- c(cdsTail, threeHead)
    merged <- unlist(partitions, use.names = FALSE)
    stop <- reduceKeepAttr(split(merged,
      names(partitions)[groupings(partitions)]))
  } else {
    stop <- NULL
  }
  ss <- list(start, stop)
  names(ss) <- c("starts", "stops")
  return(ss)
}


#' Find if there is periodicity in the vector
#'
#' @param x (numeric) Vector of values to detect periodicity of 3 like in
#' RiboSeq data.
#' @return a logical, if it is periodic.
#' @importFrom stats fft spec.pgram
#'
isPeriodic <- function(x) {
  amplitudes <- abs(fft(x))
  amp <- amplitudes[2:(length(amplitudes)/2+1)]
  periods <- 1/spec.pgram(x = x, plot = FALSE)$freq
  return((periods[which.max(amp)] > 2.9) & (periods[which.max(amp)] < 3.1))
}


#' Get the offset for specific RiboSeq read width
#' @param x a vector with points to analyse, assumes the zero is in the
#' middle + 1
#' @param feature (character) either "start" or "stop"
#' @return a single numeric offset
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

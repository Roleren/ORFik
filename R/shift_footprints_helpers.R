#' Compute coverage for every GRangesList subset.
#'
#' This is similar to \code{\link[GenomicFeatures]{coverageByTranscript}}, but
#' it adds: automatic sorting of the windows, fix for some rare cases when
#' subsetting fails on minus/plus strands and security that subseting of windows
#' will always return values (zeros) istead of out of bounds error.
#'
#' Minus strand is already flipped so that the most 5' position on the window
#' is the first position in the returned Rle.
#' @param x the cigar of the reads
#' @param windows (GRangesList) of transcripts or CDS or other ranges that will
#' be subseting coverage of \code{x}
#' @param ignore.strand (logical) Whether to consider all reads to be "*".
#' @return (RleList) of positional counts of \code{x} ranges overlapping each
#' consecutive position of the elements of \code{windows}
#' @examples
#' cds <- GenomicRanges::GRangesList(
#'   GenomicRanges::GRanges(seqnames = "chr1",
#'                          ranges = IRanges::IRanges(100, 200),
#'                          strand = "+"))
#' reads <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges =  IRanges::IRanges(c(100, 150), c(110, 160)),
#'   strand = "+")
#' ORFik:::coverageByWindow(reads, cds)
#'
coverageByWindow <- function(x, windows, ignore.strand = FALSE) {
  if (!methods::is(windows, "GRangesList")) {
    windows <- try(exonsBy(windows, by="tx", use.names=TRUE),
                   silent=TRUE)
    if (methods::is(windows, "try-error"))
      stop("failed to extract the exon ranges ", "from 'windows' with ",
           "exonsBy(windows, by=\"tx\", use.names=TRUE)")
  }
  seqinfo(x) <- merge(seqinfo(x), seqinfo(windows))
  if (!S4Vectors::isTRUEorFALSE(ignore.strand))
    stop("'ignore.strand' must be TRUE or FALSE")

  ## Fix seqlengths so that windows will always be able to subset on x coverage
  tDT <- data.table::as.data.table(windows)
  tDT <- tDT[, list(seqlengths = max(end)), by = "seqnames"]
  xDT <- data.table::as.data.table(x)
  xDT <- xDT[, list(seqlengths = max(end)), by = "seqnames"]
  cDT <- rbind(xDT, tDT)
  cDT <- cDT[, list(seqlengths = max(seqlengths)), by = "seqnames"]
  seqlengths(x) <- cDT$seqlengths[match(seqlevels(x), cDT$seqnames)]

  # sort windows
  if (is.null(names(windows))) names(windows) <- seq_along(windows)
  windows <- sortPerGroup(windows)

  ## 1) Compute unique exons ('uex').

  ex <- unlist(windows, use.names = FALSE)
  ## We could simply do 'uex <- unique(ex)' here but we're going to need
  ## 'sm' and 'is_unique' later to compute the "reverse index" so we compute
  ## them now and use them to extract the unique exons. That way we hash
  ## 'ex' only once (the expensive operation).
  sm <- selfmatch(ex)  # uses a hash table internally
  is_unique <- sm == seq_along(sm)
  uex2ex <- which(is_unique)  # index of unique exons
  uex <- ex[uex2ex]  # unique exons

  ## 2) Compute coverage for each unique exon ('uex_cvg').

  #There doesn't seem to be much benefit in doing this.
  #x <- subsetByOverlaps(x, windows, ignore.strand=TRUE)
  if (ignore.strand) {
    cvg <- coverage(x)
    uex_cvg <- cvg[uex]  # parallel to 'uex'
  } else {
    x1 <- x[strand(x) %in% c("+", "*")]
    x2 <- x[strand(x) %in% c("-", "*")]
    cvg1 <- coverage(x1)
    cvg2 <- coverage(x2)
    is_plus_ex <- strand(uex) == "+"
    is_minus_ex <- strand(uex) == "-"
    if (!identical(is_plus_ex, !is_minus_ex))
      stop("'windows' has exons on the * strand. ",
           "This is not supported at the moment.")
    uex_cvg <- RleList(as.list(rep(0, length(uex))))
    uex_cvg[is_plus_ex] <- cvg1[uex[is_plus_ex]]
    uex_cvg[is_minus_ex] <- cvg2[uex[is_minus_ex]]
  }

  ## 3) Flip coverage for exons on minus strand.

  ## It feels like this is not as fast as it could be (the bottleneck being
  ## subsetting an Rle object which needs to be revisited at some point).
  uex_cvg <- IRanges::revElements(uex_cvg, strand(uex) == "-")

  ## 4) Compute coverage by original exon ('ex_cvg').

  ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]  # reverse index
  #stopifnot(identical(ex2uex[uex2ex], seq_along(uex2ex)))  # sanity
  #stopifnot(identical(ex2uex[sm], ex2uex))  # sanity
  #stopifnot(all(uex[ex2uex] == ex))  # sanity

  ex_cvg <- uex_cvg[ex2uex]  # parallel go 'ex'

  ## 5) Compute coverage of each transcript by concatenating coverage of its
  ##    exons.

  ans <- IRanges:::regroupBySupergroup(ex_cvg, windows)

  ## 6) Propagate 'mcols(windows)'.

  mcols(ans) <- mcols(windows)
  ans
}


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
  cvg <- coverageByWindow(x, windows)
  cvg <- Reduce(`+`, cvg) / length(cvg)

  hitMap <- data.frame(avg_counts = as.vector(cvg),
                       position = -window_size:(window_size - 1),
                       frame = c(rev(rep_len(3:1, window_size)),
                                 rep_len(1:3, window_size)))
  return(hitMap)
}


#' Shift ribo-seq reads using cigar string
#'
#' @param cigar the cigar of the reads
#' @param shift the shift as numeric
#' @param is_plus_strand logical
#' @return the shifted read
#'
parseCigar <- function(cigar, shift, is_plus_strand) {
  c_signs <- unlist(explodeCigarOps(cigar))
  c_counts <- unlist(explodeCigarOpLengths(cigar))

  i = ifelse(is_plus_strand, 0, length(c_signs) + 1)
  increment = ifelse(is_plus_strand, 1, -1)
  limit = 0
  plusShift = 0
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
#'
txNamesWithLeaders <- function(txdb, minFiveUTR = 30L,
                               minCDS = 150L, minThreeUTR = 30L) {
  if(class(txdb) != "TxDb") stop("txdb must be a TxDb object")

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
#' For each cds in \code{txdb} object, filtered by \code{txNames},
#' get a window around start and stop codons within \code{window_size}
#' downstream and upstream of the codon.
#' @param txdb a txdb object of annotations
#' @param txNames a character vector of the transcript names to use
#' @param start (logical) whether to include start codons
#' @param stop (logical) whether to incude stop codons
#' @param window_size (integer) size of the window to extract upstream of the
#' start/stop codon and downstream
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
  txdb, txNames, start = TRUE, stop = TRUE, window_size = 30L) {
  cds <- GenomicFeatures::cdsBy(txdb, by = "tx", use.names = TRUE)[txNames]
  cdsTiled <- tile1(groupGRangesBy(unlist(cds, use.names = TRUE)))

  if (start) {
    fiveUTRs <- GenomicFeatures::fiveUTRsByTranscript(
      txdb, use.names=TRUE)[txNames]
    cdsHead <- phead(cdsTiled, window_size)
    fiveTiled <- tile1(groupGRangesBy(unlist(fiveUTRs, use.names = TRUE)))
    fiveTail <- ptail(fiveTiled, window_size)
    merged <- unlist(c(fiveTail, cdsHead), use.names = FALSE)
    start <- reduce(split(merged, names(merged)))
  } else {
    start <- NULL
  }

  if (stop) {
    threeUTRs <- GenomicFeatures::threeUTRsByTranscript(
      txdb, use.names=TRUE)[txNames]
    cdsTail <- ptail(cdsTiled, window_size)
    threeTiled <- tile1(groupGRangesBy(unlist(threeUTRs, use.names = TRUE)))
    threeHead <- phead(threeTiled, window_size)
    merged <- unlist(c(cdsTail, threeHead), use.names = FALSE)
    stop <- reduce(split(merged, names(merged)))
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
  periods <- 1/spec.pgram(x = x, plot = F)$freq
  return((periods[which.max(amp)] > 2.9) & (periods[which.max(amp)] < 3.1))
}


#' Get the offset for specific RiboSeq read width
#' @param x a vector with points to analyse, assumes the zero is in the
#' middle + 1
#' @param feature (character) either "start" or "stop"
#' @return a single numeric offset
#'
changePointAnalysis <- function(x, feature = "start") {
  meta <- x[1:40]
  pos <- -(length(x)/2):(length(x)/2 - 1)
  if (feature == "start") {
    means <- c()
    for (j in 15:35) {
      m <- mean(meta[j:40]) - mean(meta[1:(j-1)])
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

#' Get a binned coverage window per transcript
#'
#' Per transcript (or other regions), bin them all to windowSize (default 100),
#' and make a data.table, rows are positions, useful for plotting with ORFik
#' and ggplot2.
#'
#' NOTE: All ranges with smaller width than windowSize, will of course be
#' removed. What is the 100th position on a 1 width object ?
#' @param txdb a TxDb object or a path to gtf/gff/db file.
#' @param reads GRanges or GAlignment of reads
#' @param splitIn3 a logical(TRUE), split window in 3 (leader, cds, trailer)
#' @param windowSize an integer (100), size of windows (columns). All genes with
#' region smaller than this size are filter out for metacoverage.
#' @param fraction a character (1), info on reads (which read length,
#' or which type (RNA seq)) (row names)
#' @inheritParams coveragePerTiling
#' @param BPPARAM how many cores/threads to use? default: bpparam()
#' @return a data.table with columns position, score
#' @importFrom GenomicFeatures fiveUTRsByTranscript cdsBy threeUTRsByTranscript
#' exonsBy
#' @keywords internal
windowPerTranscript <- function(txdb, reads, splitIn3 = TRUE,
                                windowSize = 100, fraction = "1",
                                weight = "score", drop.zero.dt = FALSE,
                                BPPARAM = bpparam()) {
  # Load data
  txdb <- loadTxdb(txdb)

  if (splitIn3) {
    # Lets take all valid transcripts, with size restrictions:
    # leader > 100 bases, cds > 100 bases, trailer > 100 bases
    txNames <- filterTranscripts(txdb, windowSize, windowSize,
                                 windowSize) # valid transcripts
    leaders <- loadRegion(txdb, "leaders")[txNames]
    cds <- loadRegion(txdb, "cds")[txNames]
    trailers <- loadRegion(txdb, "trailers")[txNames]
    txCov <- splitIn3Tx(leaders, cds, trailers, reads, windowSize, fraction,
                        weight, is.sorted = TRUE, drop.zero.dt = drop.zero.dt,
                        BPPARAM = BPPARAM)

  } else {
    tx <- loadRegion(txdb, "tx")
    txNames <- widthPerGroup(tx) >= windowSize
    if (!any(txNames)) stop(paste0("no valid transcripts with length",
                                   windowSize))
    tx <- tx[txNames]
    txCov <- scaledWindowPositions(tx, reads, windowSize, weight = weight,
                                   is.sorted = TRUE, drop.zero.dt = drop.zero.dt)
    txCov[, `:=` (fraction = fraction, feature = "transcript")]
  }
  txCov[] # for print
  return(txCov)
}

#' Create binned coverage of transcripts, split into the 3 parts.
#'
#' The 3 parts  of transcripts are the leaders, the cds' and trailers.
#' Per transcript part, bin them all to windowSize (default 100),
#' and make a data.table, rows are positions, useful for plotting with ORFik
#' and ggplot2.
#' @param leaders a \code{\link{GRangesList}} of leaders (5' UTRs)
#' @param cds a \code{\link{GRangesList}} of coding sequences
#' @param trailers a \code{\link{GRangesList}} of trailers (3' UTRs)
#' @param is.sorted logical (FALSE), is grl sorted. That is + strand groups in
#' increasing ranges (1,2,3), and - strand groups in decreasing ranges (3,2,1)
#' @inheritParams windowPerTranscript
#' @return a data.table with columns position, score
#' @keywords internal
splitIn3Tx <- function(leaders, cds, trailers, reads, windowSize = 100,
                       fraction = "1", weight = "score",
                       is.sorted = FALSE, drop.zero.dt = FALSE,
                       BPPARAM = BiocParallel::SerialParam()) {
  features <- c("leaders", "cds", "trailers")
  features <- features[c(length(leaders) > 0, length(cds) > 0, length(trailers) > 0)]
  txCov <- bplapply(features, FUN = function(feature, reads, leaders, cds,
                                    trailers, windowSize, weight, drop.zero.dt) {
    cov <- scaledWindowPositions(get(feature, mode = "S4"), reads, windowSize,
                                 weight = weight, is.sorted = is.sorted,
                                 drop.zero.dt = drop.zero.dt)
    cov[, `:=` (feature = feature)]
  },  reads = reads, leaders = leaders, cds = cds,trailers = trailers,
  windowSize = windowSize, drop.zero.dt = drop.zero.dt, weight = weight,
  BPPARAM = BPPARAM)

  txCov <- rbindlist(txCov)
  txCov[, `:=` (fraction = fraction)]
  return(txCov)
}

#' Calculate meta-coverage of reads around input GRanges/List object.
#'
#' Sums up coverage over set of GRanges objects as a meta representation.
#' @param x GRanges/GAlignment object of your reads.
#' Remember to resize them beforehand to width of 1 to focus on
#' 5' ends of footprints etc, if that is wanted.
#' @param windows GRangesList or GRanges of your ranges
#' @param scoring a character, default: "sum", one of
#' (zscore, transcriptNormalized, mean, median, sum, sumLength, NULL),
#' see ?coverageScorings for info and more alternatives.
#' @param withFrames a logical (TRUE), return positions with the 3 frames,
#' relative to zeroPosition. zeroPosition is frame 0.
#' @param zeroPosition an integer DEFAULT (NULL), the point if all windows
#' are equal size, that should be set to position 0. Like leaders and
#' cds combination, then 0 is the TIS and -1 is last base in leader. NOTE!:
#' if not all windows have equal width, this will be ignored. If all have
#' equal width and zeroPosition is NULL, it is set to as.integer(width / 2).
#' @param scaleTo an integer (100), if windows have different size,
#'  a meta window can not directly be created, since a meta window must
#'  have equal size for all windows. Rescale (bin) all windows to scaleTo.
#'  i.e c(1,2,3) -> size 2 -> coverage of position c(1, mean(2,3)) etc.
#' @param fraction a character/integer (NULL), the fraction i.e
#' (27) for read length 27, or ("LSU") for large sub-unit TCP-seq.
#' @param feature a character string, info on region. Usually either
#' gene name, transcript part like cds, leader, or CpG motifs etc.
#' @param forceUniqueEven, a logical (TRUE), if TRUE; require that all windows
#' are of same width and even. To avoid bugs. FALSE if score is NULL.
#' @param forceRescale logical, default TRUE. If TRUE, if
#' \code{unique(widthPerGroup(windows))} has length > 1, it will force all
#' windows to width of the \code{scaleTo} argument, making a binned meta
#' coverage.
#' @inheritParams coveragePerTiling
#' @param append.zeroes logical, default FALSE. If TRUE and drop.zero.dt
#' is TRUE and all windows have equal length, it will add back 0 values after transformation.
#' Sometimes needed for correct plots, if TRUE, will call abort if not all
#' windows are equal length!
#' @return A data.table with scored counts (score) of
#' reads mapped to positions (position) specified in windows along with
#' frame (frame) per gene (genes) per library (fraction) per transcript region
#' (feature). Column that does not apply is not given, but position and (score/count)
#' is always returned.
#' @family coverage
#' @export
#' @importFrom BiocGenerics Reduce
#' @examples
#' library(GenomicRanges)
#' windows <- GRangesList(GRanges("chr1", IRanges(c(50, 100), c(80, 200)),
#'                                "-"))
#' x <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges =  IRanges::IRanges(c(100, 180), c(200, 300)),
#'   strand = "-")
#' metaWindow(x, windows, withFrames = FALSE)
#'
metaWindow <- function(x, windows, scoring = "sum", withFrames = FALSE,
                       zeroPosition = NULL, scaleTo = 100,
                       fraction = NULL, feature = NULL,
                       forceUniqueEven = !is.null(scoring),
                       forceRescale = TRUE,
                       weight = "score", drop.zero.dt = FALSE,
                       append.zeroes = FALSE) {
  window_size <- unique(widthPerGroup(windows))
  if (!is.null(zeroPosition) & !is.numeric(zeroPosition))
    stop("zeroPosition must be numeric / integer if defined!")
  if (length(window_size) != 1 & forceUniqueEven)
    stop("All input 'windows' should have the same sum(width(windows)),
          when forceUniqueEven is TRUE")
  if (forceUniqueEven & any(window_size %% 2 != 0))
    stop("Width of the window has to be even number, when forceUniqueEven
          is TRUE")
  if ((length(window_size) != 1) & append.zeroes) {
    stop("Can not append zeroes when windows are of different size, run without optimization!")
  }
  if (!is.null(scoring)) {
    if ((drop.zero.dt | append.zeroes) & (scoring %in% c("periodic", "periodicv"))) {
      stop("append.zeroes / drop.zero.dt optimization with score as periodic is not allowed!")
    }
  }

  if ((length(window_size) != 1) & forceRescale) {
    hitMap <- scaledWindowPositions(windows, x, scaleTo, weight = weight,
                                    drop.zero.dt = drop.zero.dt)
  } else {
    hitMap <- coveragePerTiling(windows, x, is.sorted = TRUE,
                                keep.names = FALSE, as.data.table = TRUE,
                                withFrames = FALSE, weight = weight,
                                drop.zero.dt = drop.zero.dt)
    zeroPosition <- ifelse(is.null(zeroPosition), as.integer(window_size / 2),
                           as.integer(zeroPosition))
    hitMap[, position := position - (zeroPosition + 1L) ]
  }
  # Special "if", for debug periodicity
  if (!is.null(scoring)) {
    if ((scoring == "periodicv") & !is.null(fraction)) {
      hitMap[, fraction := rep(fraction, nrow(hitMap))]
    }
  }

  hitMap <- coverageScorings(hitMap, scoring, copy.dt = FALSE)
  if (append.zeroes) { # Add back zero values now
    hitMap <- appendZeroes(hitMap, max.pos = window_size - zeroPosition -1,
                           min.pos = -zeroPosition, fractions = NULL)
  }

  if (withFrames & length(window_size) == 1) {
    hitMap[, frame := c(rev(rep_len(seq.int(2L, 0L), zeroPosition)),
                        rep_len(seq.int(0L, 2L), window_size - zeroPosition))]
  }
  if (!is.null(fraction)) {
    hitMap[, fraction := rep(fraction, nrow(hitMap))]
  }
  if (!is.null(feature)) {
    hitMap[, feature := rep(feature, nrow(hitMap))]
  }

  hitMap[] # for print
  return(hitMap)
}

#' Scale (bin) windows to a meta window of given size
#'
#' For example scale a coverage table of a all human CDS to width 100
#'
#' Nice for making metaplots, the score will be mean of merged positions.
#' @param scaleTo an integer (100), if windows have different size,
#'  a meta window can not directly be created, since a meta window must
#'  have equal size for all windows. Rescale all windows to scaleTo.
#'  i.e c(1,2,3) -> size 2 -> c(1, mean(2,3)) etc. Can also be a vector,
#'  1 number per grl group.
#' @inheritParams coverageScorings
#' @inheritParams coveragePerTiling
#' @return A data.table with scored counts (counts) of
#' reads mapped to positions (position) specified in windows along with
#' frame (frame).
#' @family coverage
#' @export
#' @examples
#' library(GenomicRanges)
#' windows <- GRangesList(GRanges("chr1", IRanges(1, 200), "-"))
#' x <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges =  IRanges::IRanges(c(1, 100, 199), c(2, 101, 200)),
#'   strand = "-")
#' scaledWindowPositions(windows, x, scaleTo = 100)
#'
scaledWindowPositions <- function(grl, reads, scaleTo = 100,
                                  scoring = "meanPos", weight = "score",
                                  is.sorted = FALSE,
                                  drop.zero.dt = FALSE) {
  if ((length(scaleTo) != 1) & length(scaleTo) != length(grl))
    stop("length of scaleTo must either be 1 or length(grl)")

  cov <- coveragePerTiling(grl, reads, is.sorted = is.sorted, keep.names = FALSE,
                           as.data.table = TRUE, withFrames = FALSE,
                           weight = weight, drop.zero.dt = drop.zero.dt)
  cov[, scalingFactor := (scaleTo/widthPerGroup(grl, FALSE))[genes]]
  cov[, position := ceiling(scalingFactor * position)]

  if (length(scaleTo) == 1) {
    cov[position > scaleTo, position := scaleTo]
  } else {
    if (any(cov[, .(max = max(position)), by = genes]$max > scaleTo)) {
      index <- cov[, .(ind = which.max(position)), by = genes]$ind
      cov[index, position := scaleTo]
    }
  }

  # mean counts per position per group
  return(coverageScorings(cov, scoring, copy.dt = FALSE))
}

#' Add a coverage scoring scheme
#'
#' Different scorings and groupings of a coverage representation.
#'
#' Usually output of metaWindow or scaledWindowPositions is input in this
#' function.
#'
#' Content of coverage data.table:
#' It must contain the count and position columns.
#'
#' genes column: If you have multiple windows, the genes column must define
#' which gene/transcript grouping the different counts belong to. If there is
#' only a meta window or only 1 gene/transcript, then this column is
#' not needed.
#'
#' fraction column: If you have coverage of i.e RNA-seq and Ribo-seq, or TCP
#' -seq of large and small subunite, divide into fractions.
#' Like factor(RNA, RFP)
#'
#' feature column: If gene group is subdivided into parts, like gene is
#' transcripts, and feature column can be c(leader, cds, trailer) etc.
#'
#' Given a data.table coverage of counts, add a scoring scheme.
#' per: the grouping given, if genes is defined,
#' group by per gene in default scoring.\cr
#' Scorings:\cr
#'
#' *  zscore (count-mean(window))/sd(window) per) # Outlier detection
#' *  modzscore (count-median(window))/mad(window) per) # Outlier detection for signal with extreme values
#' *  transcriptNormalized (sum(count / sum of counts per)) this is scaled per group
#' *  mean (mean(count per))
#' *  median (median(count per))
#' *  sum (count per)
#' *  log2sum (count per)
#' *  log10sum (count per)
#' *  sumLength (count per) / number of windows
#' *  meanPos (mean per position per gene) used in scaledWindowPositions
#' *  sumPos (sum per position per gene) used in scaledWindowPositions
#' *  frameSum (sum per frame per gene) used in ORFScore
#' *  frameSumPerL (sum per frame per read length)
#' *  frameSumPerLG (sum per frame per read length per gene)
#' *  fracPos (fraction of counts per position per gene) non scaled version of transcriptNormalized
#' *  periodic (Fourier transform periodicity of meta coverage per fraction)
#' *  NULL (no grouping, return input directly)
#' @md
#' @param coverage a data.table containing at least columns (count, position),
#' it is possible to have additionals: (genes, fraction, feature)
#' @param scoring a character, one of (zscore, transcriptNormalized,
#' mean, median, sum, log2sum, log10sum, sumLength, meanPos and frameSum,
#' periodic, NULL). More info in details
#' @param copy.dt logical TRUE, copy object, to avoid overwriting original object.
#' Set to false to run function using reference to object,
#' a speed up if original object is not needed.
#' @return a data.table with new scores (size dependent on score used)
#' @family coverage
#' @export
#' @examples
#' dt <- data.table::data.table(count = c(4, 1, 1, 4, 2, 3),
#'                              position = c(1, 2, 3, 4, 5, 6))
#' coverageScorings(dt, scoring = "zscore")
#'
#' # with grouping gene
#' dt$genes <- c(rep("tx1", 3), rep("tx2", 3))
#' coverageScorings(dt, scoring = "zscore")
#'
coverageScorings <- function(coverage, scoring = "zscore",
                             copy.dt = TRUE) {
  if (is.null(scoring)) return(coverage)
  cov <- if (copy.dt) {
    setDT(copy(coverage))
  } else coverage

  # if (!is.null(coverage$genes)) {
  #   if (is(coverage$genes, "factor")) {
  #     coverage[, genes := as.integer(genes)]
  #   }
  # }
  # find groupings
  groupGF <- coverageGroupings(c(is.null(cov$fraction),
                               is.null(cov$genes)))
  groupFPF <- coverageGroupings(c(is.null(cov$fraction),
                                is.null(cov$feature)), "FPF")
  if (is.null(cov$count)) cov[, count := score]
  if (scoring == "meanPos") { # rare scoring schemes
    groupFPF <- quote(list(genes, position))
    scoring <- "mean"
  } else if (scoring == "sumPos") { # rare scoring schemes
    groupFPF <- quote(list(genes, position))
    scoring <- "sum"
  } else if (scoring == "frameSum") {
    if (is.null(cov$frame))
      stop("Can not use frameSum scoring when no frames are given!")
    groupFPF <- quote(list(genes, frame))
    scoring <- "sum"
  } else if (scoring == "frameSumPerL") {
    if (is.null(cov$frame))
      stop("Can not use frameSum scoring when no frames are given!")
    groupFPF <- quote(list(fraction, frame))
    scoring <- "sum"
  } else if (scoring == "frameSumPerLG") {
    if (is.null(cov$frame))
      stop("Can not use frameSum scoring when no frames are given!")
    groupFPF <- quote(list(fraction, genes, frame))
    scoring <- "sum"
  } else if (scoring %in% c("periodic", "periodicv")) {
    groupFPF <- ifelse(!is.null(cov$fraction),
                                (quote(fraction)), quote(""))
  }

  # create score
  if (scoring == "zscore") {
    # z score over transcript per fraction
    cov[, `:=` (windowSD = sd(count), windowMean = mean(count)),
        by = eval(groupGF)]
    cov[, zscore := (count-windowMean)/windowSD]
    # create mean and sum scores per position, per feature
    res <- cov[, .(score = mean(zscore, na.rm = TRUE)),
                        by = eval(groupFPF)]
  } else if (scoring == "modzscore") {
    # z score over transcript per fraction
    cov[, `:=` (windowMAD = mad(count), windowMedian = median(count)),
        by = eval(groupGF)]
    cov[, modzscore := (count-windowMedian)/windowMAD]
    # create mean and sum scores per position, per feature
    res <- cov[, .(score = mean(modzscore, na.rm = TRUE)),
               by = eval(groupFPF)]
  } else if (scoring == "transcriptNormalized") {
    cov[, `:=` (gene_sum = sum(count)), by = eval(groupGF)]
    res <- cov[, .(score = sum(count / gene_sum, na.rm = TRUE)),
                        by = eval(groupFPF)]
  } else if (scoring == "fracPos") {
    cov[, `:=` (gene_sum = sum(count)), by = eval(groupGF)]
    res <- cov[, `:=` (score = count / gene_sum)]

  }  else if (scoring == "mean") {
    res <- cov[, .(score = mean(count, na.rm = TRUE)),
                        by = eval(groupFPF)]
  } else if (scoring == "median") {
    res <- cov[, .(score = median(count, na.rm = TRUE)),
                        by = eval(groupFPF)]
  } else if (scoring == "sum") {
    res <- cov[, .(score = sum(count, na.rm = TRUE)),
                        by = eval(groupFPF)]
  } else if (scoring == "log2sum") {
    res <- cov[, .(score = log2(sum(count, na.rm = TRUE))),
               by = eval(groupFPF)]
  } else if (scoring == "log10sum") {
    res <- cov[, .(score = log10(sum(count, na.rm = TRUE))),
               by = eval(groupFPF)]
  } else if (scoring == "sumLength") {
    uniques <- ifelse(!is.null(cov$genes),
                      length(unique(cov$genes)), 1)
    res <- cov[, .(score = sum(count, na.rm = TRUE)),
                        by = eval(groupFPF)]
    res[, `:=`  (score = score / uniques)]
  } else if (scoring %in% c("periodic", "periodicv")) {
    cov <- coverageScorings(cov, "transcriptNormalized")
    verbose <- scoring == "periodicv"
    if (!is.null(cov$fraction)) {
      res <- cov[, .(score = isPeriodic(score, unique(fraction), verbose = verbose)),
                 by = eval(groupFPF)]
    } else {
      res <- cov[, .(score = isPeriodic(score, verbose = verbose)),
                 by = eval(groupFPF)]
    }

  } else stop(paste("Invalid scoring: ", scoring))
  res[] # for print
  return(res)
}

#' Get coverage per group
#'
#' It tiles each GRangesList group to width 1, and finds hits per position.\cr
#' A range from 1:5 will split into c(1,2,3,4,5) and count hits on each.
#' This is a safer speedup of coverageByTranscript from GenomicFeatures.
#' It also gives the possibility to return as data.table, for faster
#' computations.
#'
#' NOTE: If reads contains a $score column, it will presume that this is
#' the number of replicates per reads, weights for the
#' coverage() function.
#' So delete the score column or set weight to something else if this
#' is not wanted.
#' @param grl a \code{\link{GRangesList}} of 5' utrs, CDS, transcripts, etc.
#' @param reads a \code{\link{GAlignments}}, \code{\link{GRanges}}, or
#' precomputed coverage as \code{\link{covRle}} (one for each strand) of
#' RiboSeq, RnaSeq etc.\cr Weigths for scoring is default the 'score'
#' column in 'reads'. Can also be random access paths to bigWig or fstwig file.
#' Do not use random access for more than a few genes, then loading the entire files
#' is usually better. File streaming is still in beta, so use with care!
#' @param is.sorted logical (FALSE), is grl sorted. That is + strand groups in
#' increasing ranges (1,2,3), and - strand groups in decreasing ranges (3,2,1)
#' @param keep.names logical (TRUE), keep names or not. If as.data.table is TRUE,
#' names (genes column) will be a factor column, if FALSE it will be an
#' integer column (index of gene), so first input grl element is 1. Dropping names
#' gives ~ 20 \% speedup. If drop.zero.dt is FALSE, data.table will not
#' return names, will use index (to avoid memory explosion).
#' @param as.data.table a logical (FALSE), return as data.table with 2 columns,
#' position and count.
#' @param withFrames a logical (FALSE), only available if as.data.table is
#' TRUE, return the ORF frame, 1,2,3, where position 1 is 1, 2 is 2 and
#' 4 is 1 etc.
#' @param weight (default: 'score'), if defined a character name
#' of valid meta column in subject. GRanges("chr1", 1, "+", score = 5),
#' would mean score column tells that this alignment region was found 5 times.
#' Formats which loads a score column like this:
#' Bigwig, wig, ORFik ofst, collapsed bam, bedoc and .bedo.
#' As do CAGEr CAGE files and many other package formats.
#' You can also assign a score column manually.
#' @param drop.zero.dt logical FALSE, if TRUE and as.data.table is TRUE,
#' remove all 0 count positions.
#' This greatly speeds up and most importantly, greatly reduces memory usage.
#' Will not change any plots, unless 0 positions are used in some sense.
#' (mean, median, zscore coverage will only scale differently)
#' @param fraction integer or character, a description column. Useful for grouping
#' multiple outputs together. If returned as Rle, this
#' is added as: metadata(coverage) <- list(fraction = fraction). If as.data.table
#' it will be added as an additional column.
#' @return a numeric RleList, one numeric-Rle per group with # of hits per position.
#' Or data.table if as.data.table is TRUE,
#' with column names c("count" [numeric or integer], "genes" [integer], "position" [integer])
#' @importFrom GenomicFeatures coverageByTranscript
#' @export
#' @family ExtendGenomicRanges
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20),
#'                                 end = c(5, 15, 25)),
#'                strand = "+")
#' grl <- GRangesList(tx1_1 = ORF)
#' RFP <- GRanges("1", IRanges(25, 25), "+")
#' coveragePerTiling(grl, RFP, is.sorted = TRUE)
#' # now as data.table with frames
#' coveragePerTiling(grl, RFP, is.sorted = TRUE, as.data.table = TRUE,
#'                   withFrames = TRUE)
#' # With score column (usually replicated reads on that position)
#' RFP <- GRanges("1", IRanges(25, 25), "+", score = 5)
#' dt <- coveragePerTiling(grl, RFP, is.sorted = TRUE,
#'                         as.data.table = TRUE, withFrames = TRUE)
#' class(dt$count) # numeric
#' # With integer score column (faster and less space usage)
#' RFP <- GRanges("1", IRanges(25, 25), "+", score = 5L)
#' dt <- coveragePerTiling(grl, RFP, is.sorted = TRUE,
#'                         as.data.table = TRUE, withFrames = TRUE)
#' class(dt$count) # integer
#'
coveragePerTiling <- function(grl, reads, is.sorted = FALSE,
                               keep.names = TRUE, as.data.table = FALSE,
                               withFrames = FALSE, weight = "score",
                               drop.zero.dt = FALSE, fraction = NULL) {
  if (!is.null(fraction)) stopifnot(length(fraction) == 1)
  if (!is.sorted) grl <- sortPerGroup(grl)

  if (is(reads, "character")) {
    return(coverage_random_access_file(reads, grl, withFrames, fraction))
  }

  if (is(reads, "covRle") | is(reads, "RleList")) {
    coverage <- coverageByTranscriptC(reads, grl)
  } else {
    score.defined <- is.numeric(weight) | (weight[1] %in% colnames(mcols(reads)))
    if (score.defined) {
      seqinfo.is.correct <- identical(seqlengths(reads), seqlengths(grl)) & !anyNA(seqlengths(grl))
      coverage <- coverageByTranscriptW(reads, grl, weight = weight,
                                        seqinfo.x.is.correct = seqinfo.is.correct)
    } else coverage <- coverageByTranscript(reads, grl)
  }
  if (!keep.names) names(coverage) <- NULL

  if (as.data.table) {
    return(coverage_to_dt(coverage, keep.names = keep.names,
                          withFrames = withFrames, weight = weight,
                          drop.zero.dt = drop.zero.dt, fraction = fraction))
  }
  if (!is.null(fraction)) metadata(coverage) <- list(fraction = fraction)

  return(coverage)
}

coverage_random_access_file <- function(reads, grl, withFrames, fraction = NULL) {
  if (all(is.null(fraction))) fraction <- "all"
  file_ext <- tools::file_ext(reads)
  if (all(file_ext %in% c("bigWig", "bw", "bigWig"))) {
    if (length(grl) > 1) stop("bigwig random access only supported for single gene for now!")
    rl <- ranges(grl)
    names(rl) <- seqnamesPerGroup(grl, FALSE)
    strands <- strandPerGroup(grl, FALSE)
    coverage <- NumericList()
    for (strand in unique(strands)) {
      strand_now <- ifelse(strand == "+", 1,2)
      if (strand_now == 2) {
        if (length(reads) == 1) strand_now <- 1
        temp_cov <- import.bw(reads[strand_now], as = "NumericList",
                              which = rl[strands == strand])
        temp_cov2 <- rev(temp_cov)
        temp_cov2 <- relist(rev(unlist(temp_cov2)), skeleton = temp_cov)
      } else {
        temp_cov2 <- import.bw(reads[strand_now], as = "NumericList",
                              which = rl[strands == strand])
      }
      coverage <- c(coverage, temp_cov2)
    }

    # To data.table
    coverage <- data.table(count = unlist(coverage, use.names = FALSE))
    is_integers <- all(round(coverage$count, 0) == coverage$count)
    if (is_integers) coverage[, count := as.integer(count)]
    coverage[, genes := rep.int(1L, nrow(coverage))]

  } else if (all(file_ext == "fstwig") | TRUE) {
    chrs <- seqnamesPerGroup(grl, FALSE)
    chr_groups <- unique(chrs)
    for (chr in chr_groups) {
      coverage <- import.fstwig(unlist(grl[chrs == chr_groups], use.names = FALSE),
                                dir = reads,
                                readlengths = fraction)
    }
    coverage[, genes := rep(seq(length(grl)), each = widthPerGroup(grl, FALSE))]
  }

  coverage[, position := seq_len(.N), by = genes]
  if (withFrames) coverage[, frame := (position - 1) %% 3]
  coverage[]
  return(coverage)
}

#' Convert coverage RleList to data.table
#' @param coverage RleList with names
#' @inheritParams coveragePerTiling
#' @return a data.table with column names c("count" [numeric or integer],
#'  "genes" [integer], "position" [integer])
coverage_to_dt <- function(coverage, keep.names = TRUE,
                           withFrames = FALSE, weight = "score",
                           drop.zero.dt = FALSE, fraction = NULL) {
  if (drop.zero.dt) { # Drop 0 count rows
    if (!keep.names) names(coverage) <- seq.int(length(coverage))
    runV <- unlist(runValue(coverage), use.names = FALSE)
    ir <- unlist(ranges(coverage), use.names = TRUE)[runV > 0]
    if (length(ir) == 0) {
      count <- data.table(count = integer(), genes = integer(), position = integer())
      if (withFrames) count[, frame := integer()]
      if (!is.null(fraction)) {
        count[, fraction := new(class(fraction))]
        warning("No coverage found for fraction: ", fraction, ". Returning empty data.table!")
      } else warning("No coverage found, Returning empty data.table!")
      return(count)
    }
    runV <- runV[runV > 0]
    runV <- rep.int(runV, width(ir))
    ir <- unlist(tile(ir, width = 1))
    count <- data.table(count = runV) #
    if (keep.names) {
      count[, `:=`(genes = factor(names(ir)), position = start(ir))]
    } else {
      count[, `:=`(genes = as.integer(names(ir)),
                   position = start(ir))]
    }
  } else { # Keep all 0 count rows
    window_size <- unique(lengths(coverage))
    count.is.integer <- is(coverage@unlistData@values, "integer")

    count <- if (count.is.integer) {
      data.table(count = unlist(IntegerList(coverage), use.names = FALSE))
    } else data.table(count = unlist(NumericList(coverage),
                                     use.names = FALSE))

    count[, genes := ORFik::groupings(coverage)]
    if (length(window_size) != 1) { # different size windows
      count[, position := seq_len(.N), by = genes]
    } else { # all same size
      count[, position := rep.int(seq.int(window_size), length(coverage))]
    }
  }
  if (withFrames) {
    count[, frame := (position - 1) %% 3]
  }
  if (!is.null(fraction)) count[, fraction := fraction]

  count[]# for print
  return(count)
}

#' Append zero values to data.table
#'
#' For every position in width max.pos - min.pos + 1, append
#' 0 values in data.table. Needed when coveragePerTiling was run
#' on coverage window with drop.zero.dt as TRUE and you need to plot
#' 0 positions after a transformation by coverageScorings.
#' @param dt a data.table from coverageByTiling that is
#' normalized by coverageScorings.
#' @param max.pos integer, max position of dt
#' @param min.pos integer, default 1L. Minimum position of dt
#' @param fractions default unique(dt$fraction), will repeat
#' each fraction max.pos - min.pos + 1 times.
#' @return a data.table with appended 0 values
#' @keywords internal
appendZeroes <- function(dt, max.pos, min.pos = 1L,
                          fractions = unique(dt$fraction)) {
  width <- max.pos - min.pos + 1L
  if (!is.null(fractions)) {
    by.columns <- c("fraction", "position")
    temp <- data.table(position = seq.int(min.pos, max.pos), fraction = rep(fractions, each = width))
  } else {
    by.columns <- c("position")
    temp <- data.table(position = seq.int(min.pos, max.pos))
  }
  if (!is.null(dt$frame)) {
    temp[, frame := (position %% 3) - 1]
  }

  temp <- data.table::merge.data.table(temp, dt, by = by.columns,
                                     all.x = TRUE, all.y = FALSE)
  if (!is.null(temp$score)) {
    temp[is.na(score), score := 0]
  }
  if (!is.null(temp$count)) {
    temp[is.na(count), count := 0]
  }

  data.table::setcolorder(temp, c(colnames(dt)))
  return(temp)
}

#' Find proportion of reads per position per read length in region
#'
#' This is defined as:
#' Given some transcript region (like CDS), get coverage per position.
#' By default only returns positions that have hits, set drop.zero.dt
#' to FALSE to get all 0 positions.
#'
#' @inheritParams windowPerReadLength
#' @param withFrames logical TRUE, add ORF frame (frame 0, 1, 2), starting
#' on first position of every grl.
#' @param exclude.zero.cov.grl logical, default TRUE. Do not include
#' ranges that does not have any coverage (0 reads on them),
#' this makes it faster to run.
#' @param drop.zero.dt logical, default TRUE.
#' If TRUE and as.data.table is TRUE, remove all 0 count positions.
#' This greatly speeds up and most importantly, greatly reduces memory usage.
#' Will not change any plots, unless 0 count positions are used in some sense.
#' @param BPPARAM how many cores/threads to use? default: bpparam()
#' @return a data.table with lengths by coverage.
#' @family coverage
#' @importFrom data.table rbindlist
#' @export
#' @examples
#' # Raw counts per gene per position
#' cds <- GRangesList(tx1 = GRanges("1", 100:129, "+"))
#' reads <- GRanges("1", seq(79,129, 3), "+")
#' reads$size <- 28 # <- Set read length of reads
#' regionPerReadLength(cds, reads, scoring = NULL)
#' ## Sum up reads in each frame per read length per gene
#' regionPerReadLength(cds, reads, scoring = "frameSumPerLG")
regionPerReadLength <- function(grl, reads, acceptedLengths = NULL,
                                withFrames = TRUE,
                                scoring = "transcriptNormalized",
                                weight = "score",
                                exclude.zero.cov.grl = TRUE,
                                drop.zero.dt = TRUE,
                                BPPARAM = bpparam()) {
  rWidth <- readWidths(reads)
  all_lengths <- sort(unique(rWidth))
  if (!is.null(acceptedLengths))
    all_lengths <- all_lengths[all_lengths %in% acceptedLengths]
  if (exclude.zero.cov.grl) {
    if (!is(reads, "covRleList")) {
      grl <- grl[hasHits(grl, reads)]
    } else message("NOTE: exclude.zero.cov.grl not implemented for covRleList yet, ignoring it.")
  }

  dt <- bplapply(all_lengths, function(l, grl, reads, weight, rWidth,
                                       scoring, withFrames, drop.zero.dt) {
    d <- coveragePerTiling(grl, reads = if (!is(reads, "covRleList")) {reads[rWidth == l]}
                                           else  reads@list[[as.character(l)]],
                           as.data.table = TRUE, withFrames = withFrames,
                           weight = weight, is.sorted = TRUE,
                           drop.zero.dt = drop.zero.dt)
    d[, fraction := l]
    return(coverageScorings(d, scoring, copy.dt = FALSE))
  }, grl = grl, reads = reads, weight = weight, rWidth = rWidth,
     scoring = scoring, withFrames = withFrames, drop.zero.dt = drop.zero.dt,
     BPPARAM = BPPARAM)

  return(rbindlist(dt))
}

#' Find proportion of reads per position per read length in window
#'
#' This is defined as:
#' Fraction of reads  per read length, per position in whole window (defined
#' by upstream and downstream)
#' If tx is not NULL, it gives a metaWindow, centered around startSite of
#' grl from upstream and downstream. If tx is NULL, it will use only downstream
#' , since it has no reference on how to find upstream region.
#' The exception is when upstream is negative, that is,
#' going into downstream region of the object.
#'
#' Careful when you create windows where not all transcripts are long enough,
#' this function usually is used first with filterTranscripts to make sure they
#' are of all of valid length!
#' @inheritParams startRegion
#' @inheritParams coveragePerTiling
#' @param reads a \code{\link{GAlignments}}, \code{\link{GRanges}}, or
#' precomputed coverage as \code{\link{covRleList}}
#' (where names of covRle objects are readlengths) of
#' RiboSeq, RnaSeq etc. \cr Weigths for scoring is default the 'score'
#' column in 'reads'. Can also be random access paths to bigWig or fstwig file.
#' Do not use random access for more than a few genes, then loading the entire files
#' is usually better.
#' @param pShifted a logical (TRUE), are Ribo-seq reads p-shifted to size
#'  1 width reads? If upstream and downstream is set, this argument is
#'  irrelevant. So set to FALSE if this is not p-shifted Ribo-seq.
#' @param upstream an integer (5), relative region to get upstream from. Default:
#' \code{ifelse(!is.null(tx), ifelse(pShifted, 5, 20), min(ifelse(pShifted, 5, 20), 0))}
#' @param downstream an integer (20), relative region to get downstream from. Default:
#' \code{ifelse(pShifted, 20, 5)}
#' @param acceptedLengths an integer vector (NULL), the read lengths accepted.
#'  Default NULL, means all lengths accepted.
#' @param zeroPosition an integer DEFAULT (upstream), what is the center point?
#' Like leaders and cds combination, then 0 is the TIS and -1 is last base in leader.
#' NOTE!: if windows have different widths, this will be ignored.
#' @param scoring a character (transcriptNormalized), which meta coverage scoring ?
#' one of (zscore, transcriptNormalized, mean, median, sum, sumLength, fracPos),
#' see ?coverageScorings for more info. Use to decide a scoring of hits
#' per position for metacoverage etc. Set to NULL if you do not want meta coverage,
#' but instead want per gene per position raw counts.
#' @param append.zeroes logical, default FALSE. If TRUE and drop.zero.dt
#' is TRUE and all windows have equal length, it will add back 0 values after transformation.
#' Sometimes needed for correct plots, if TRUE, will call abort if not all
#' windows are equal length!
#' @param windows the GRangesList windows to actually check, default:
#' \code{startRegion(grl, tx, TRUE, upstream, downstream)}.
#' @return a data.table with 4 columns: position (in window), score,
#' fraction (read length). If score is NULL, will also return genes (index of grl).
#' A note is that if no coverage is found, it returns an empty data.table.
#' @family coverage
#' @importFrom data.table rbindlist
#' @export
#' @examples
#' cds <- GRangesList(tx1 = GRanges("1", 100:129, "+"))
#' tx <- GRangesList(tx1 = GRanges("1", 80:129, "+"))
#' reads <- GRanges("1", seq(79,129, 3), "+")
#' windowPerReadLength(cds, tx, reads, scoring = "sum")
#' windowPerReadLength(cds, tx, reads, scoring = "transcriptNormalized")
windowPerReadLength <- function(grl, tx = NULL, reads, pShifted = TRUE,
                                upstream = ifelse(!is.null(tx), ifelse(pShifted, 5, 20),
                                                  min(ifelse(pShifted, 5, 20), 0)),
                                downstream = ifelse(pShifted, 20, 5),
                                acceptedLengths = NULL,
                                zeroPosition = upstream,
                                scoring = "transcriptNormalized",
                                weight = "score", drop.zero.dt = FALSE,
                                append.zeroes = FALSE,
                                windows = startRegion(grl, tx, TRUE, upstream, downstream)) {
  if(length(reads) == 0 | length(grl) == 0) {
    return(data.table())
  }
  if (is.null(tx)) upstream <- min(upstream, 0)
  if(!grl_has_any_valid_lengths(windows, upstream + downstream + 1))
    return(data.table())

  rWidth <- readWidths(reads)
  all_lengths <- sort(unique(rWidth))
  if (!is.null(acceptedLengths))
    all_lengths <- all_lengths[all_lengths %in% acceptedLengths]

  dt <- data.table()
  for(l in all_lengths) {
    dt <- rbindlist(list(dt, metaWindow(
      x = if (!is(reads, "covRleList")) {reads[rWidth == l]} else  reads@list[[as.character(l)]],
      windows = windows, scoring = scoring,
      zeroPosition =  zeroPosition, forceUniqueEven = FALSE, fraction = l,
      weight = weight,
      drop.zero.dt = drop.zero.dt, append.zeroes = append.zeroes)))
  }

  dt[] # for print
  return(dt)
}

grl_has_any_valid_lengths <- function(windows, windowSize) {
  noHits <- widthPerGroup(windows) < windowSize
  if (all(noHits)) {
    warning("No object in GRangesList had valid width!")
    message(paste("First window should be:", windowSize, "it is:",
                  widthPerGroup(windows[1])))
    print(grl[noHits][1])
    message("Maybe you forgot extendLeaders()?")
    return(FALSE)
  }
  return(TRUE)
}

#' Get grouping for a coverage table in ORFik
#'
#' Either of two groupings:
#' GF: Gene, fraction
#' FGF: Fraction, position, feature
#' It finds which of these exists, and auto groups
#'
#' Normally not used directly!
#' @param logicals size 2 logical vector, the is.null checks for each column,
#' @param grouping which grouping to perform, default "GF"
#' Gene & Fraction grouping. Alternative "FGF", Fraction & position & feature.
#' @return a quote of the grouping to pass to data.table
#' @keywords internal
coverageGroupings <- function(logicals, grouping = "GF") {
  one <- !logicals[1]
  two <- !logicals[2]
  if (grouping == "GF") { # Gene/Fraction
    groupGF <-  quote("")
    if (!(one & two)) {
      if (two) {
        groupGF <- quote(list(genes))
      } else if (one)
        groupGF <- quote(list(fraction))
    } else groupGF <- quote(list(genes, fraction))
    return(groupGF)
  } else if (grouping == "FPF") { # Fraction/Position/Feature
    groupFPF <- quote(list(position))
    if (!(one & two)) {
      if (two) {
        groupFPF <- quote(list(position, feature))
      } else if (one)
        groupFPF <- quote(list(fraction, position))
    } else groupFPF <- quote(list(fraction, position, feature))
  } else stop("undefined grouping of coverage")
  return(groupFPF)
}

#' Get number of genes per coverage table
#'
#' Used to count genes in ORFik meta plots
#' @param coverage a data.table with coverage
#' @return number of genes in coverage
#' @keywords internal
getNGenesCoverage <- function(coverage) {
  if (is.null(coverage$genes)) return(0)
  if (is.null(coverage$fraction)) {
    n <- coverage[, .(nGenes = max(genes))]
  } else {
    n <- coverage[, .(nGenes = max(genes)), by = fraction]
  }

  if (nrow(n) == 0) return(0)
  return(n$nGenes)
}



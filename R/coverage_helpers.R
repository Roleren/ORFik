#' Calculate meta-coverage of reads around input GRanges/List object.
#'
#' Sums up coverage over set of GRanges objects as a meta representation.
#' @param x GRangesList/GRanges object of your reads.
#' Remember to resize them beforehand to width of 1 to focus on
#' 5' ends of footprints, if that is wanted.
#' @param windows GRangesList or GRanges of your ranges
#' @param withFrames a logical (TRUE), return positions with the 3 frames,
#' relative to zeroPosition. zeroPosition is frame 1.
#' @param zeroPosition an integer DEFAULT (NULL), the point if all windows
#' are equal size, that should be set to position 0. Like leaders and
#' cds combination, then 0 is the TIS and -1 is last base in leader. NOTE!:
#' if windows have different widths, this will be ignored.
#' @param scaleTo an integer (100), if windows have different size,
#'  a meta window can not directly be created, since a meta window must
#'  have equal size for all windows. Rescale all windows to scaleTo.
#'  i.e c(1,2,3) -> size 2 -> c(1, sum(2,3)) etc.
#' @param returnAs a character (data.frame), do data.table for speed.
#' @param fraction a character/integer (NULL), the fraction i.e
#' (27) for read length 27, or ("LSU") for large sub-unit TCP-seq.
#' @param forceUniqueEven, a logical (TRUE), require that all windows
#' are of same width and even. To avoid bugs.
#' @return A data.frame or data.table with average counts (avg_counts) of
#' reads mapped to positions (position) specified in windows along with
#' frame (frame).
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
metaWindow <- function(x, windows, withFrames = TRUE, zeroPosition = NULL,
                       scaleTo = 100, returnAs = "data.frame",
                       fraction = NULL, forceUniqueEven = TRUE) {
  window_size <- unique(sum(width(windows)))
  if (!is.null(zeroPosition) & !is.numeric(zeroPosition))
    stop("zeroPosition must be numeric if defined")
  if (length(window_size) != 1 & forceUniqueEven)
    stop("All input 'windows' should have the same sum(width(windows)),
          when forceUniqueEven is TRUE")
  if ((window_size %% 2 != 0) & forceUniqueEven)
    stop("Width of the window has to be even number, when forceUniqueEven
          is TRUE")

  if (length(window_size) != 1) {
    hitMap <- scaledWindowPositions(windows, x, scaleTo)
  } else {
    coverage <- coverageByTranscript(x, windows)
    hitMap <- data.table(Count = unlist(IntegerList(coverage),
                                       use.names = FALSE))
    zeroPosition <- ifelse(is.null(zeroPosition), (window_size)/2,
                           zeroPosition)

    hitMap[, position := rep.int(seq.int(window_size) - (zeroPosition + 1),
                                length(coverage))]
    hitMap <- hitMap[, .(counts = sum(Count)), by = position]
  }

  if (withFrames & length(window_size) == 1) {
    hitMap[, frame := c(rev(rep_len(3L:1L, zeroPosition)),
                        rep_len(seq.int(3), window_size - zeroPosition))]
  }
  if (!is.null(fraction)) {
    hitMap[, fraction := rep(fraction, nrow(hitMap))]
  }
  if(returnAs == "data.frame")
    return(setDF(hitMap))
  return(hitMap)
}

#' Scale windows to a meta window of size
#'
#' For example scale a coverage plot of a all human CDS to width 100
#'
#' Nice for making metaplots
#' @param grl GRangesList or GRanges of your ranges
#' @param reads GRanges object of your reads.
#' @param scaleTo an integer (100), if windows have different size,
#'  a meta window can not directly be created, since a meta window must
#'  have equal size for all windows. Rescale all windows to scaleTo.
#'  i.e c(1,2,3) -> size 2 -> c(1, mean(2,3)) etc.
#' @param scoring a character, one of (mean, median, sum)
#' @return A data.frame or data.table with average counts (avg_counts) of
#' reads mapped to positions (position) specified in windows along with
#' frame (frame).
#'
scaledWindowPositions <- function(grl, reads, scaleTo = 100,
                                  scoring = "mean") {
  coverage <- coverageByTranscript(reads, grl)

  grlGroup <- groupings(coverage)
  widths <- lengths(coverage) # BIOC her <-

  count <- data.table(Count = unlist(IntegerList(coverage), use.names = FALSE))
  count[, ones := rep.int(1L, length(grlGroup))]
  count[, grlGroup := grlGroup]
  count[, cumSum := cumsum(ones), by= grlGroup]
  count[, scalingFactor := (scaleTo/widths)[grlGroup]]
  count[, scalingGroup := ceiling(scalingFactor * cumSum)]
  count[scalingGroup > scaleTo]$scalingGroup <- scaleTo

  # mean counts per position per group
  if (scoring == "mean") {
    count <- count[, .(count = mean(Count)), by = list(grlGroup, scalingGroup)]
  } else if (scoring == "median") {
    count <- count[, .(count = median(Count)), by = list(grlGroup,
                                                         scalingGroup)]
  } else if (scoring == "sum") {
    count <- count[, .(count = sum(Count)), by = list(grlGroup, scalingGroup)]
  } else stop(paste("Invalid scoring: ", scoring))

  return(count)
}

#' Add a coverage scoring scheme
#'
#' Different scorings and groupings of a coverage representation.
#'
#' Usually output of metaWindow or scaledWindowCoverage is input in this
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
#' per: the grouping given, if genes is defined, group by per gene in scoring.
#' Scorings:
#' 1. zscore (count-windowMean)/windowSD per)
#' 2. transcriptNormalized (sum(count / sum of counts per))
#' 3. mean (mean(count per))
#' 4. median (median(count per))
#' 5. sum (count per)
#' 6. sumLength (count per) / number of windows
#' @param coverage a data.table containing at least columns (count, position),
#' it is possible to have additionals: (genes, fraction, feature)
#' @param scoring a character, one of (zscore, transcriptNormalized,
#' mean, median, sum, sumLength)
#' @return a data.table with new scores
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
coverageScorings <- function(coverage, scoring = "zscore") {

  coverage_sub <- coverage
  # find groupings
  groupGF <-  quote("")
  groupFPF <- quote(list(position))
  if (!(!is.null(coverage$fraction) & !is.null(coverage$genes))) {
    if (!is.null(coverage$genes)) {
      groupGF <- quote(list(genes))
    } else if (!is.null(coverage$fraction))
        groupGF <- quote(list(fraction))
  } else groupGF <- quote(list(genes, fraction))

  if (!(!is.null(coverage$fraction) & !is.null(coverage$feature))) {
    if (!is.null(coverage$feature)) {
      groupFPF <- quote(list(position, feature))
    } else if (!is.null(coverage$fraction))
      groupFPF <- quote(list(fraction, position))
  } else groupFPF <- quote(list(fraction, position, feature))

  # create score
  if (scoring == "zscore") {
    # z score over transcript per fraction
    coverage_sub[, `:=` (windowSD = sd(count), windowMean = mean(count))
                        , by = eval(groupGF)]
    coverage_sub[, zscore := (count-windowMean)/windowSD]
    # create mean and sum scores per position, per feature
    res <- coverage_sub[, .(score = mean(zscore, na.rm = TRUE)),
                        by = eval(groupFPF)]
  } else if (scoring == "transcriptNormalized") {
    coverage_sub[, `:=` (gene_sum = sum(count)), by = eval(groupGF)]
    res <- coverage_sub[, .(score = sum(count / gene_sum, na.rm = TRUE)),
                        by = eval(groupFPF)]

  } else if (scoring == "mean") {
    res <- coverage_sub[, .(score = mean(count, na.rm = TRUE)),
                        by = eval(groupFPF)]
  } else if (scoring == "median") {
    res <- coverage_sub[, .(score = median(count, na.rm = TRUE)),
                        by = eval(groupFPF)]
  } else if (scoring == "sum") {
    res <- coverage_sub[, .(score = sum(count, na.rm = TRUE)),
                        by = eval(groupFPF)]
  } else if (scoring == "sumLength") {
    uniques <- ifelse(!is.null(coverage_sub$genes),
                      length(unique(coverage_sub$genes)), 1)
    res <- coverage_sub[, .(score = sum(count, na.rm = TRUE)),
                        by = eval(groupFPF)]
    res[, `:=`  (score = score / uniques)]
  } else stop(paste("Invalid scoring: ", scoring))

  return(res)
}

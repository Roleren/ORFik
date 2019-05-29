#' Get coverage window per transcript
#'
#' NOTE: All ranges with smaller width than windowSize, will of course be
#' removed. What is the 100 position on a 1 width size ?
#' @param txdb a TxDb object or a path to gtf/gff/db file.
#' @param reads GRanges or GAlignment of reads
#' @param splitIn3 a logical(TRUE), split window in 3 (leader, cds, trailer)
#' @param windowSize an integer (100), size of windows
#' @param fraction a character (1), info on reads (which read length,
#' or which type (RNA seq))
#' @return a data.table with columns position, score
windowPerTranscript <- function(txdb, reads, splitIn3 = TRUE,
                                windowSize = 100, fraction = "1") {
  # Load data
  txdb <- loadTxdb(txdb)

  if (splitIn3) {
    # Lets take all valid transcripts, with size restrictions:
    # leader > 100 bases, cds > 100 bases, trailer > 100 bases
    txNames <- filterTranscripts(txdb, windowSize, windowSize,
                                 windowSize) # valid transcripts
    leaders = fiveUTRsByTranscript(txdb, use.names = TRUE)[txNames]
    cds <- cdsBy(txdb, "tx", use.names = TRUE)[txNames]
    trailers = threeUTRsByTranscript(txdb, use.names = TRUE)[txNames]
    txCov <- splitIn3Tx(leaders, cds, trailers, reads, fraction)

  } else {
    tx <- exonsBy(txdb, by = "tx", use.names = TRUE)
    txNames <- widthPerGroup(tx) >= windowSize
    if (!any(txNames)) stop(paste0("no valid transcripts with length",
                                   windowSize))
    tx <- tx[txNames]
    txCov <- scaledWindowPositions(tx, reads, windowSize)
    txCov[, `:=` (fraction = fraction, feature = "transcript")]
  }
  txCov[] # for print
  return(txCov)
}

#' Create coverage of transcripts, split into the 3 parts.
#'
#' The 3 parts  of transcripts are the leaders, the cds' and trailers.
#' @param leaders a \code{\link{GRangesList}} of leaders (5' UTRs)
#' @param cds a \code{\link{GRangesList}} of coding sequences
#' @param trailers a \code{\link{GRangesList}} of trailers (3' UTRs)
#' @inheritParams windowPerTranscript
#' @return a data.table with columns position, score
splitIn3Tx <- function(leaders, cds, trailers, reads, windowSize = 100,
                       fraction = "1") {
  leaderCov <- scaledWindowPositions(leaders, reads, windowSize)
  leaderCov[, `:=` (feature = "leaders")]
  cdsCov <- scaledWindowPositions(cds, reads, windowSize)
  cdsCov[, `:=` (feature = "cds")]
  trailerCov <- scaledWindowPositions(trailers, reads, windowSize)
  trailerCov[, `:=` (feature = "trailers")]

  txCov <- rbindlist(list(leaderCov, cdsCov, trailerCov))
  txCov[, `:=` (fraction = fraction)]
  return(txCov)
}

#' Calculate meta-coverage of reads around input GRanges/List object.
#'
#' Sums up coverage over set of GRanges objects as a meta representation.
#' @param x GRangesList/GRanges object of your reads.
#' Remember to resize them beforehand to width of 1 to focus on
#' 5' ends of footprints, if that is wanted.
#' @param windows GRangesList or GRanges of your ranges
#' @param scoring a character, one of (zscore, transcriptNormalized,
#' mean, median, sum, sumLength, NULL), see ?coverageScorings
#' @param withFrames a logical (TRUE), return positions with the 3 frames,
#' relative to zeroPosition. zeroPosition is frame 0.
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
#' @param feature a character string, info on region. Usually either
#' gene name, transcript part like cds, leader, or CpG motifs etc.
#' @param forceUniqueEven, a logical (TRUE), if TRUE; require that all windows
#' are of same width and even. To avoid bugs. FALSE if score is NULL.
#' @return A data.frame or data.table with scored counts (score) of
#' reads mapped to positions (position) specified in windows along with
#' frame (frame).
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
                       returnAs = "data.frame", fraction = NULL,
                       feature = NULL, forceUniqueEven = !is.null(scoring)) {
  window_size <- unique(widthPerGroup(windows))
  if (!is.null(zeroPosition) & !is.numeric(zeroPosition))
    stop("zeroPosition must be numeric if defined")
  if (length(window_size) != 1 & forceUniqueEven)
    stop("All input 'windows' should have the same sum(width(windows)),
          when forceUniqueEven is TRUE")
  if (forceUniqueEven & any(window_size %% 2 != 0))
    stop("Width of the window has to be even number, when forceUniqueEven
          is TRUE")


  if (length(window_size) != 1) {
    hitMap <- scaledWindowPositions(windows, x, scaleTo)
  } else {
    hitMap <- coveragePerTiling(windows, x, is.sorted = TRUE,
                                keep.names = TRUE, as.data.table = TRUE,
                                withFrames = FALSE)
    zeroPosition <- ifelse(is.null(zeroPosition), (window_size)/2,
                           zeroPosition)
    hitMap[, position := position - (zeroPosition + 1) ]
  }

  hitMap <- coverageScorings(hitMap, scoring)

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

  if (returnAs == "data.frame") {
    hitMap <- setDF(hitMap)
    return(hitMap)
  }
  hitMap[] # for print
  return(hitMap)
}

#' Scale windows to a meta window of size
#'
#' For example scale a coverage plot of a all human CDS to width 100
#'
#' Nice for making metaplots, the score will be mean of merged positions.
#' @param grl GRangesList or GRanges of your ranges
#' @param reads GRanges object of your reads.
#' @param scaleTo an integer (100), if windows have different size,
#'  a meta window can not directly be created, since a meta window must
#'  have equal size for all windows. Rescale all windows to scaleTo.
#'  i.e c(1,2,3) -> size 2 -> c(1, mean(2,3)) etc. Can also be a vector,
#'  1 number per grl group.
#' @param scoring a character, one of (meanPos, sumPos)
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
                                  scoring = "meanPos") {
  if ((length(scaleTo) != 1) & length(scaleTo) != length(grl))
    stop("length of scaleTo must either be 1 or length(grl)")

  count <- coveragePerTiling(grl, reads, is.sorted = FALSE, keep.names = TRUE,
                              as.data.table = TRUE, withFrames = FALSE)
  count[, scalingFactor := (scaleTo/widthPerGroup(grl, FALSE))[genes]]
  count[, position := ceiling(scalingFactor * position)]

  if (length(scaleTo) == 1) {
    count[position > scaleTo]$position <- scaleTo
  } else {
    if (any(count[, .(max = max(position)), by = genes]$max > scaleTo)) {
      index <- count[, .(ind = which.max(position)), by = genes]$ind
      count[index, position := scaleTo]
    }
  }

  # mean counts per position per group
  count <- coverageScorings(count, scoring)
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
#' 7. meanPos (mean per position per gene) used in scaledWindowPositions
#' 8. sumPos (sum per position per gene) used in scaledWindowPositions
#' 9. frameSum (sum per frame per gene) used in ORFScore
#' 10. fracPos (fraction of counts per position per gene)
#' 11. periodic (Fourier transform periodicity of meta coverage per fraction)
#' 12. NULL (return input directly)
#' @param coverage a data.table containing at least columns (count, position),
#' it is possible to have additionals: (genes, fraction, feature)
#' @param scoring a character, one of (zscore, transcriptNormalized,
#' mean, median, sum, sumLength, meanPos and frameSum, periodic, NULL)
#' @return a data.table with new scores
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
coverageScorings <- function(coverage, scoring = "zscore") {
  if (is.null(scoring)) return(coverage)
  cov <- setDT(copy(coverage))
  # find groupings
  groupGF <- coverageGroupings(c(is.null(cov$fraction),
                               is.null(cov$genes)))
  groupFPF <- coverageGroupings(c(is.null(cov$fraction),
                                is.null(cov$feature)), "FPF")
  if (is.null(cov$count)) cov$count <- cov$score
  if (scoring == "meanPos") { # rare scoring schemes
    groupFPF <- quote(list(genes, position))
    scoring <- "mean"
  } else if (scoring == "sumPos") { # rare scoring schemes
    groupFPF <- quote(list(genes, position))
    scoring <- "sum"
  } else if (scoring == "frameSum") {
    groupFPF <- quote(list(genes, frame))
    scoring <- "sum"
  } else if (scoring == "periodic") {
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
  } else if (scoring == "sumLength") {
    uniques <- ifelse(!is.null(cov$genes),
                      length(unique(cov$genes)), 1)
    res <- cov[, .(score = sum(count, na.rm = TRUE)),
                        by = eval(groupFPF)]
    res[, `:=`  (score = score / uniques)]
  }  else if (scoring == "periodic") {
    cov <- coverageScorings(cov, "transcriptNormalized")
    res <- cov[, .(score = isPeriodic(score)), by = eval(groupFPF)]
  } else stop(paste("Invalid scoring: ", scoring))
  res[] # for print
  return(res)
}

#' Get overlaps and convert to coverage list
#'
#' @param gr a \code{\link{GRanges}} object, to get coverage of.
#' @param reads a GAlignment or GRanges object of RiboSeq, RnaSeq etc.
#' @param keep.names logical (T), keep names or not.
#' @param type a string (any), argument for countOverlaps.
#' @return a Rle, one list per group with # of hits per position.
#' @export
#' @family ExtendGenomicRanges
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20),
#'                                 end = c(5, 15, 25)),
#'                strand = "+")
#' names(ORF) <- "tx1"
#' reads <- GRanges("1", IRanges(25, 25), "+")
#' overlapsToCoverage(ORF, reads)
#'
overlapsToCoverage <- function(gr, reads, keep.names = TRUE, type = "any") {
  counts <- countOverlaps(gr, reads, type = type)

  names <- names(counts)
  names(counts) <- NULL
  countList <- split(counts, names)
  if (!keep.names) {
    countList <- IRanges::RleList(countList)
    names(countList) <- NULL
    return(countList)
  }
  return(IRanges::RleList(countList))
}


#' Get coverage per group
#'
#' It tiles each GRangesList group, and finds hits per position
#'
#' This is a safer speedup of coverageByTranscript from GenomicFeatures.
#' It also gives the possibility to return as data.table, for faster
#' computations.
#' @param grl a \code{\link{GRangesList}}
#'  of 5' utrs or transcripts.
#' @param is.sorted logical (F), is grl sorted.
#' @inheritParams overlapsToCoverage
#' @param as.data.table a logical (FALSE), return as data.table with 2 columns,
#' position and count.
#' @param withFrames a logical (FALSE), only available if as.data.table is
#' TRUE, return the ORF frame, 1,2,3, where position 1 is 1, 2 is 2 and
#' 4 is 1 etc.
#' @return a RleList, one integer-Rle per group with # of hits per position.
#' Or data.table if as.data.table is TRUE.
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
#'
coveragePerTiling <- function(grl, reads, is.sorted = FALSE,
                              keep.names = TRUE, as.data.table = FALSE,
                              withFrames = FALSE) {
  if (!is.sorted) grl <- sortPerGroup(grl)
  if (length(grl) > 10000) { # faster version for big grl
    coverage <- coverageByTranscript(reads, grl)
    if (!keep.names) names(coverage) <- NULL
  } else {
    if (length(names(grl)) != length(unique(names(grl)))) {
      stop("grl have duplicated names, be sure to make them unique!")
    }
    unlTile <- unlistGrl(tile1(grl, matchNaming = FALSE))
    coverage <- overlapsToCoverage(unlTile, reads, keep.names = keep.names)
  }
  if (as.data.table) {
    window_size <- unique(widthPerGroup(grl))
    count <- data.table(count = unlist(IntegerList(coverage),
                                       use.names = FALSE))
    count[, genes := groupings(coverage)]
    if (length(window_size) != 1) { # different size windows
      count[, ones := rep.int(1L, length(genes))]
      count[, position := cumsum(ones), by = genes]
      count$ones <- NULL
    } else { # all same size
      count[, position := rep.int(seq.int(window_size), length(coverage))]
    }

    if (withFrames) {
      count[, frame := (position - 1) %% 3]
    }
    count[]# for print
    return(count)
  }
  return(coverage)
}

#' Find proportion of reads per position in window
#'
#' This is like a more detailed floss score, where floss score takes fraction
#' of reads per read length over whole window, this is defined as:
#' Fraction of reads  per read length, per position in whole window (by
#' upstream and downstream)
#'
#' If tx is not NULL, it gives a metaWindow, centered around startSite of
#' grl from upstream and downstream. If tx is NULL, it will use only downstream
#' , since it has no reference from to find upstream from. Unless upstream is
#' negative, that is, going downstream.
#'
#' @inheritParams startRegion
#' @param reads any type of reads, usualy ribo seq. As GAlignment, GRanges
#'  or GRangesList object.
#' @param pShifted a logical (TRUE), are riboseq reads p-shifted to size
#'  1 width reads? If upstream and downstream is set, this argument is
#'  irrelevant.
#' @param upstream an integer (5), relative region to get upstream from.
#' @param downstream an integer (20), relative region to get downstream from
#' @param acceptedLengths an integer vector (NULL), the read lengths accepted.
#'  Default NULL, means all lengths accepted.
#' @param zeroPosition an integer DEFAULT (upstream), the point if all windows
#' are equal size, that should be set to position 0. Like leaders and
#' cds combination, then 0 is the TIS and -1 is last base in leader. NOTE!:
#' if windows have different widths, this will be ignored.
#' @param scoring a character (transcriptNormalized), one of
#' (zscore, transcriptNormalized, mean, median, sum, sumLength, fracPos),
#' see ?coverageScorings. Use to choose meta coverage or per transcript.
#' @return a data.frame with lengths by coverage / vector of proportions
#' @family coverage
#' @importFrom data.table rbindlist
#' @export
#'
windowPerReadLength <- function(grl, tx = NULL, reads, pShifted = TRUE,
                                upstream = if (pShifted) 5 else 20,
                                downstream = if (pShifted) 20 else 5,
                                acceptedLengths = NULL,
                                zeroPosition = upstream,
                                scoring = "transcriptNormalized") {

  if (is.null(tx)) upstream <- min(upstream, 0)

  if(length(reads) == 0 | length(grl) == 0) {
    return(data.table())
  }
  windowSize <- upstream + downstream + 1
  windows <- startRegion(grl, tx, TRUE, upstream, downstream)
  noHits <- widthPerGroup(windows) < windowSize
  if (all(noHits)) {
    warning("no grl ranges had valid window size!")
    return(data.table())
  }

  rWidth <- readWidths(reads)
  all_lengths <- sort(unique(rWidth))
  if (!is.null(acceptedLengths))
    all_lengths <- all_lengths[all_lengths %in% acceptedLengths]
  dt <- data.table()

  for(l in all_lengths){
    dt <- rbindlist(list(dt, metaWindow(
      x = reads[rWidth == l], windows = windows, scoring = scoring,
      returnAs = "data.table", zeroPosition =  zeroPosition,
      forceUniqueEven = FALSE, fraction = l)))
  }

  dt[] # for print
  return(dt)
}

#' Get grouping for a coverage table in ORFik
#'
#' Either of two groupings:
#' GF: Gene, fraction
#' FGF: Fraction, position, feature
#' It finds which of these exists, and auto groups
#'
#' Normally not used directly
#' @param logicals size 2 logical vector, the is.null checks for each column,
#' @param grouping which grouping to perform
#' @return a quote of the grouping to pass to data.table
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

#' Get number of genes per grouping
#'
#' @param coverage a data.table with coverage
#' @return number of genes in coverage
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



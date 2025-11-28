#' Restrict GRangesList
#'
#' Will restrict GRangesList to `N` bp downstream from the first base.
#' @param grl (GRangesList)
#' @param firstN (integer) Allow only this many bp downstream, maximum.
#' @return a GRangesList of reads restricted to firstN and tiled by 1
#' @keywords internal
downstreamN <- function(grl, firstN = 150L) {
  return(heads(tile1(grl, matchNaming = FALSE, mergeEqualNamed = FALSE), firstN))
}

#' Get total widths per GRangesList group
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return an integer vector (named/unnamed) of widths
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' widthPerGroup(grl)
widthPerGroup <- function(grl, keep.names = TRUE) {
  if (length(grl) == 0) return(integer(0))
  validRL(class(grl))

  widths_raw <- if (is.grl(grl)) {grl@unlistData@ranges@width
  } else grl@unlistData@width

  res <- data.table(widths = widths_raw,
                    grouping = rep.int(seq_along(grl),
                      times = width(grl@partitioning)))[, .(widths = sum(widths)), by = grouping]$widths
  empty_groups <- length(grl) != max(nrow(res), 0)
  if (empty_groups) {
    res_temp <- res
    res <- rep(0L, length(grl))
    res[width(grl@partitioning) > 0] <- res_temp
  }
  if (keep.names) {
    names(res) <- names(grl)
  }
  return(res)
}
#' Get first seqname per GRangesList group
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @importFrom IRanges heads
#' @return a character vector or Rle of seqnames(if seqnames == T)
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' seqnamesPerGroup(grl)
seqnamesPerGroup <- function(grl, keep.names = TRUE) {
  validGRL(class(grl))
  if (keep.names) {
    return(heads(seqnames(grl), 1L))
  } else {
    return(as.character(heads(seqnames(grl), 1L)))
  }
}

#' Get first strand per GRangesList group
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a vector named/unnamed of characters
#' @importFrom IRanges heads
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' strandPerGroup(grl)
strandPerGroup <- function(grl, keep.names = TRUE) {
  validGRL(class(grl))
  if (!keep.names) {
    return(as.character(grl@unlistData@strand)[grl@partitioning@end])
  }
  return(grl@unlistData@strand[grl@partitioning@end])
}


#' Sort a GRangesList, helper.
#'
#' A helper for [sortPerGroup()].
#' A faster, more versatile reimplementation of GenomicRanges::sort()
#' Normally not used directly.
#' Groups first each group, then either decreasing or increasing
#' (on starts if byStarts == T, on ends if byStarts == F)
#' @param grl a \code{\link{GRangesList}}
#' @param decreasing should the first in each group have max(start(group))
#'   ->T or min-> default(F) ?
#' @param byStarts a logical T, should it order by starts or ends F.
#' @importFrom data.table as.data.table :=
#' @return an equally named GRangesList, where each group is sorted within
#' group.
#' @keywords internal
gSort <- function(grl, decreasing = FALSE, byStarts = TRUE) {
  if (length(grl) == 0) return(GRangesList())

  DT <- as.data.table(grl)
  DT$group_name <- NULL
  group <- NULL # for not getting warning
  if (decreasing) {
    if (byStarts) {
      DT <- DT[order(group, -start)]
    } else {
      DT <- DT[order(group, -end)]
    }
  } else {
    if (byStarts) {
      DT <- DT[order(group, start)]
    } else {
      DT <- DT[order(group, end)]
    }
  }
  # TODO: test naming, this is still not perfect
  testName <- names(unlist(grl[1], use.names = FALSE)[1])
  if (!is.null(testName)) {
    DT[, grnames := names(unlist(grl, use.names = FALSE))]
  }

  asgrl <- makeGRangesListFromDataFrame(
    DT, split.field = "group",
    names.field = if(is.null(testName)) NULL else "grnames",
    keep.extra.columns = TRUE)

  names(asgrl) <- names(grl)

  return(asgrl)
}


#' Sort a GRangesList
#'
#' A faster, more versatile reimplementation of
#' \code{\link{sort.GenomicRanges}} for GRangesList,
#' needed since the original works poorly for more than 10k groups.
#' This function sorts each group, where "+" strands are
#' increasing by starts and "-" strands are decreasing by ends.
#'
#' Note: will not work if groups have equal names.
#' @param grl a \code{\link{GRangesList}}
#' @param ignore.strand a boolean, (default FALSE): should minus strands be
#' sorted from highest to lowest ends. If TRUE: from lowest to highest ends.
#' @param quick.rev default: FALSE, if TRUE, given that you know all ranges are
#' sorted from min to max for both strands, it will only reverse coordinates for
#' minus strand groups, and only if they are in increasing order. Much quicker
#' @return an equally named GRangesList, where each group is
#'  sorted within group.
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(14, 7), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(1, 4), c(3, 9)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' sortPerGroup(grl)
#'
sortPerGroup <- function(grl, ignore.strand = FALSE, quick.rev = FALSE){
  if (quick.rev) {
    return(reverseMinusStrandPerGroup(grl))
  }
  if (!ignore.strand) {
    indicesPos <- strandBool(grl)

    grl[indicesPos] <- gSort(grl[indicesPos])
    grl[!indicesPos] <- gSort(grl[!indicesPos], decreasing = TRUE,
                              byStarts = FALSE)
    return(grl)
  }
  return(gSort(grl))
}

#' Reverse minus strand
#'
#' Reverse minus strand per group in a GRangesList
#' Only reverse if minus strand is in increasing order
#' @param grl a \code{\link{GRangesList}}
#' @param onlyIfIncreasing logical, default (TRUE), only reverse if decreasing
#' @return a \code{\link{GRangesList}}
#' @keywords internal
reverseMinusStrandPerGroup <- function(grl, onlyIfIncreasing = TRUE) {
  minus <- !strandBool(grl)
  if (onlyIfIncreasing) {
    minGrl <- grl[minus & numExonsPerGroup(grl, FALSE) > 1]
    if (length(minGrl) == 0) return(grl)
    decreasing <- start(minGrl[[1]])[1] > start(minGrl[[1]])[2]
    if (decreasing) return(grl)
  }
  oldGrl <- rev(grl)
  oldGrl[rev(minus)]@unlistData@ranges <- rev(grl[minus]@unlistData@ranges)
  return(rev(oldGrl))
}

#' Get first exon per GRangesList group
#'
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a \code{\link{GRangesList}}
#' @return a GRangesList of the first exon per group
#' @importFrom IRanges heads
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' firstExonPerGroup(grl)
firstExonPerGroup <- function(grl) {
  validGRL(class(grl))
  return(heads(grl, 1L))
}


#' Get last exon per GRangesList group
#'
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a \code{\link{GRangesList}}
#' @return a GRangesList of the last exon per group
#' @importFrom IRanges tails
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' lastExonPerGroup(grl)
lastExonPerGroup <- function(grl) {
  validGRL(class(grl))
  return(tails(grl, 1L))
}


#' Get first start per granges group
#'
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a Rle(keep.names = TRUE), or integer vector(FALSE)
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' firstStartPerGroup(grl)
firstStartPerGroup <- function(grl, keep.names = TRUE) {
  validGRL(class(grl))
  if (keep.names) {
    return(start(firstExonPerGroup(grl)))
  }else {
    return(as.integer(start(firstExonPerGroup(grl))))
  }
}


#' Get first end per granges group
#'
#' grl must be sorted, call ORFik:::sortPerGroup if needed
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a Rle(keep.names = T), or integer vector(F)
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' firstEndPerGroup(grl)
firstEndPerGroup <- function(grl, keep.names = TRUE) {
  validGRL(class(grl))
  if (keep.names) {
    return(end(firstExonPerGroup(grl)))
  } else {
    return(as.integer(end(firstExonPerGroup(grl))))
  }
}


#' Get last end per granges group
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a Rle(keep.names = T), or integer vector(F)
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' lastExonEndPerGroup(grl)
lastExonEndPerGroup <- function(grl,keep.names = TRUE) {
  validGRL(class(grl))
  if (keep.names) {
    return(end(lastExonPerGroup(grl)))
  } else {
    return(as.integer(end(lastExonPerGroup(grl))))
  }
}


#' Get last start per granges group
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not, default: (TRUE)
#' @return a Rle(keep.names = T), or integer vector(F)
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' lastExonStartPerGroup(grl)
lastExonStartPerGroup <- function(grl, keep.names = TRUE) {
  validGRL(class(grl))
  if (keep.names) {
    return(start(lastExonPerGroup(grl)))
  } else {
    return(as.integer(start(lastExonPerGroup(grl))))
  }
}


#' Get list of the number of exons per group
#'
#' Can also be used generaly to get number of GRanges object
#'  per GRangesList group
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a logical, keep names or not, default: (TRUE)
#' @return an integer vector of counts
#' @export
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' numExonsPerGroup(grl)
numExonsPerGroup <- function(grl, keep.names = TRUE) {
  validGRL(class(grl))
  return(lengths(grl, keep.names))
}

#' Safe unlist
#'
#' Same as [AnnotationDbi::unlist2()], keeps names correctly.
#' Two differences is that if grl have no names, it will not
#' make integer names, but keep them as null. Also if the GRangesList has names
#' , and also the GRanges groups, then the GRanges group names will be kept.
#' @param grl a GRangesList
#' @return a GRanges object
#' @export
#' @examples
#' ORF <- GRanges(seqnames = "1",
#'                ranges = IRanges(start = c(1, 10, 20),
#'                                 end = c(5, 15, 25)),
#'                strand = "+")
#' # GRL named, GR not named
#' grl <- GRangesList(tx1_1 = ORF)
#' res <- unlistGrl(grl)
#' res_original <- unlist(grl)
#' identical(res, res_original) # TRUE
#'
#' # GRL not named, GR not named
#' grl_no_names <- grl
#' names(grl_no_names) <- NULL
#' unlistGrl(grl_no_names)
#' res <- unlistGrl(grl_no_names)
#' res_original <- unlist(grl_no_names)
#' identical(res, res_original) # TRUE
#'
#' # GRL named, GR named
#' grl_names_gr_names <- unlistGrl(grl)
#' grl_names_gr_names <-
#'  split(grl_names_gr_names, names(grl_names_gr_names))
#' res <- unlistGrl(grl_names_gr_names)
#' res_original <- unlist(grl_names_gr_names)
#' identical(res, res_original) # FALSE
#'
#' # GRL not named, GR named
#' grl_not_names_gr_names <- grl_names_gr_names
#' names(grl_not_names_gr_names) <- NULL
#' res <- unlistGrl(grl_not_names_gr_names)
#' res_original <- unlist(grl_not_names_gr_names)
#' identical(res, res_original) # TRUE
#'
#'
unlistGrl <- function(grl) {
  validGRL(class(grl))
  if (length(grl) == 0) return(.unlistGrl(grl))

  grl_is_named <- !is.null(names(grl))
  res <- .unlistGrl(grl)
  if (grl_is_named) {
    gr_is_not_named <- is.null(names(res[1]))
    if (gr_is_not_named) {
      names(res) <- rep(names(grl), width(grl@partitioning))
      return(res)
    }
  }
  return(res)
}

.unlistGrl <- function(grl) grl@unlistData


#' Removes meta columns
#'
#' @param grl a \code{\link{GRangesList}} or GRanges object
#' @return same type and structure as input without meta columns
#' @keywords internal
removeMetaCols <- function(grl) {
  wasGRL <- FALSE
  if (!is.gr_or_grl(class(grl))) {
    stop("Can only remove meta columns from GRangesList or GRanges objects")
  }
  if (is.grl(class(grl))) {
    g <- groupings(grl)
    names <- names(grl)
    grl <- unlist(grl, use.names = FALSE)
    wasGRL <- TRUE
  }
  elementMetadata(grl)  <- S4Vectors::DataFrame(
    matrix(nrow = length(grl), ncol = 0))

  if (wasGRL) {
    grl <- groupGRangesBy(grl, g)
    names(grl) <- names
  }
  return(grl)
}

#' Convert GRangesList to character vector
#'
#' Single exon format:\cr
#' "1:14598834-14598914:+"\cr
#' Multi-exon format (exon separator: ';'):\cr
#' "1:15210514-15210562:+;1:15214895-15215025:+"
#' @param x A \code{\link{GRangesList}}
#' @param ... Not used for now, to preserve generic requirement
#' @return a character vector, 1 element per element in GRangesList
#' @export
setMethod("as.character", "GRangesList", function(x, ...) {
  if (length(x) == 0) return(character())
  u <- unlist(x, use.names = FALSE)
  per <- paste0(seqnames(u), ":", start(u), "-", end(u), ":", as.character(strand(u)))
  cl <- relist(per, x)
  res <- unstrsplit(cl, sep = ";")
  # names(res) <- names(x)
  return(res)
})

#' Convert a character vector to GRangesList
#' @param x a character vector
#' @return a GRangesList
#' @export
#' @examples
#' vec <- c("1:14598834-14598914:+", "1:15210514-15210562:+;1:15214895-15215025:+")
#' makeGRangesListFromCharacter(vec)
makeGRangesListFromCharacter <- function(x) {
  if (length(x) == 0) return(GRangesList())

  str_split <- strsplit(x, ";")
  gr <- as(unlist(str_split), "GRanges")
  res <- split(gr, groupings(str_split))
  names(res) <- names(x)
  return(res)
}

#' Get number of ranges per group as an iteration
#'
#' @param grl GRangesList
#' @return an integer vector
#' @export
#' @examples
#' grl <- GRangesList(GRanges("1", c(1, 3, 5), "+"),
#'                    GRanges("1", c(19, 21, 23), "+"))
#' ORFik::groupings(grl)
#'
groupings <- function(grl){
  if (length(grl) == 0) return(integer())
  l <- lengths(grl, use.names = FALSE)
  return(rep.int(seq.int(length(l)), l))
}

#' Reverse elements within list
#'
#' A faster version of S4Vectors::revElements
#' @param x RleList
#' @return a RleList (reversed inside list elements)
#' @keywords internal
revElementsF <- function(x) {
  b <- rev(x)
  b@unlistData <- rev(x@unlistData)
  return(rev(b))
}

#' coverageByTranscript with weights
#'
#' Extends the function with weights,
#' see \code{\link{coverageByTranscript}} for original function.
#' @param x reads (\code{\link{GRanges}}, \code{\link{GAlignments}})
#' @param transcripts \code{\link{GRangesList}}
#' @param ignore.strand a logical (default: FALSE)
#' @param weight a vector (default: 1L), if single number applies for all,
#' else it must be the string name of a defined meta column in "x",
#' that gives number of times a read was found.
#' GRanges("chr1", 1, "+", score = 5), would mean score column tells
#' that this alignment was found 5 times.
#' @param seqinfo.x.is.correct logical, default FALSE. If you know x, has
#' correct seqinfo, then you can save some computation time by setting this to
#' TRUE.
#' @importFrom S4Vectors wmsg isTRUEorFALSE
#' @importFrom GenomicFeatures exonsBy
#' @return Integer Rle of coverage, 1 per transcript
coverageByTranscriptW <- function (x, transcripts, ignore.strand = FALSE,
                                    weight = 1L, seqinfo.x.is.correct = FALSE) {
  if (!is(transcripts, "GRangesList")) {
    transcripts <- try(exonsBy(transcripts, by = "tx", use.names = TRUE),
                       silent = TRUE)
    if (is(transcripts, "try-error"))
      stop(wmsg("failed to extract the exon ranges ",
                "from 'transcripts' with ", "exonsBy(transcripts, by=\"tx\", use.names=TRUE)"))
  }
  if (!isTRUEorFALSE(ignore.strand))
    stop(wmsg("'ignore.strand' must be TRUE or FALSE"))
  if (!seqinfo.x.is.correct) {
    seqinfo(x) <- GenomicFeatures:::.merge_seqinfo_and_infer_missing_seqlengths(x, transcripts)
    # We create pseudo 0 lengths
    seqlengths(x)[is.na(seqlengths(x))] <- 0
  }
  return(coverageByTranscriptC(x = covRleFromGR(x, weight = weight,
                                                ignore.strand = ignore.strand),
                               transcripts = transcripts))
}

#' coverageByTranscript with coverage input
#'
#' Extends the function with direct genome coverage input,
#' see \code{\link{coverageByTranscript}} for original function.
#' @param x a covRle (one RleList for each strand in object),
#'  must have defined and correct
#' seqlengths in its SeqInfo object.
#' @param transcripts \code{\link{GRangesList}}
#' @param ignore.strand a logical (default: length(x) == 1)
#' @importFrom S4Vectors wmsg isTRUEorFALSE List
#' @importFrom GenomicFeatures exonsBy
#' @return Integer Rle of coverage, 1 per transcript
coverageByTranscriptC <- function (x, transcripts, ignore.strand = !strandMode(x)) {
  stopifnot(is(x, "covRle"))
  stopifnot(all(!anyNA(seqlengths(f(x)))))
  if (!is(transcripts, "GRangesList")) {
    transcripts <- try(exonsBy(transcripts, by = "tx", use.names = TRUE),
                       silent = TRUE)
    if (is(transcripts, "try-error"))
      stop(wmsg("failed to extract the exon ranges ",
                "from 'transcripts' with ", "exonsBy(transcripts, by=\"tx\", use.names=TRUE)"))
  }
  if (!isTRUEorFALSE(ignore.strand))
    stop(wmsg("'ignore.strand' must be TRUE or FALSE"))

  # Create hash table of unique exons as uex
  ex <- unlist(transcripts, use.names = FALSE)
  sm <- selfmatch(ex)
  is_unique <- sm == seq_along(sm)
  uex2ex <- which(is_unique)
  uex <- ex[uex2ex]

  cvg1 <- f(x)
  if (strandMode(x)) cvg2 <- r(x)

  # Get coverage
  if (ignore.strand) {
    uex_cvg <- cvg1[uex]
  }
  else {
    is_plus_ex <- strand(uex) == "+"
    is_minus_ex <- strand(uex) == "-"
    if (!identical(is_plus_ex, !is_minus_ex))
      stop(wmsg("'transcripts' has exons on the * strand. ",
                "This is not supported at the moment."))
    uex_cvg <- RleList(rep(IntegerList(1), length(uex)))
    if (length(transcripts) > 5e3) {
      names(uex) <- seq_along(uex)
      # + strand
      v <- Views(cvg1, uex[is_plus_ex])
      names(v) <- NULL
      r <- do.call(c, lapply(v[lengths(v) > 0], RleList))
      uex_cvg[as.integer(names(r))] <- r
      # - strand
      v <- Views(cvg2, uex[is_minus_ex])
      names(v) <- NULL
      r <- do.call(c, lapply(v[lengths(v) > 0], RleList))
      uex_cvg[as.integer(names(r))] <- r
    } else {
      uex_cvg[is_plus_ex] <- cvg1[uex[is_plus_ex]]
      uex_cvg[is_minus_ex] <- cvg2[uex[is_minus_ex]]
      names(uex_cvg) <- as.character(seqnames(uex))
    }
  }
  uex_cvg[strand(uex) == "-"] <- revElementsF(uex_cvg)[strand(uex) == "-"]
  ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]
  ex_cvg <- uex_cvg[ex2uex]
  ans <- IRanges:::regroupBySupergroup(ex_cvg, transcripts)
  mcols(ans) <- mcols(transcripts)
  return(ans)
}

coverageByTranscriptSum <- function(x, transcripts, ignore.strand = !strandMode(x),
                                    return_per_exon = FALSE) {
  if (!is(transcripts, "GRangesList")) {
    transcripts <- try(exonsBy(transcripts, by = "tx", use.names = TRUE), silent = TRUE)
    if (is(transcripts, "try-error"))
      stop("Failed to extract exons with exonsBy(transcripts, by='tx', use.names=TRUE)")
  }
  if (length(transcripts) == 0) return(integer())
  stopifnot(is(x, "covRle"))
  stopifnot(all(!anyNA(seqlengths(f(x)))))

  if (!isTRUEorFALSE(ignore.strand))
    stop("'ignore.strand' must be TRUE or FALSE")

  # Unique exons to avoid recomputing same ranges
  ex <- unlist(transcripts, use.names = FALSE)
  sm <- selfmatch(ex)
  is_unique <- sm == seq_along(sm)
  uex2ex <- which(is_unique)      # indices of unique exons in 'ex'
  uex <- ex[uex2ex]

  cvg1 <- f(x)                    # + strand coverage (RleList)
  cvg2 <- if (strandMode(x)) r(x) else NULL  # - strand coverage if available

  # Pre-allocate sums for unique exons
  uex_sum <- integer(length(uex))
  names(uex) <- seq_along(uex)

  if (ignore.strand) {
    stop("Not implemented yet")
    strand(uex) <- "*"
    v <- Views(cvg1, uex)
    names(v) <- NULL
    r <- unlist(lapply(v[lengths(v) > 0], sum))
    uex_sum[as.integer(names(r))] <- r
  }
  else {
    strands <- as.character(strand(uex))
    v <- Views(cvg1, uex[strands == "+"])
    names(v) <- NULL
    r <- unlist(lapply(v[lengths(v) > 0], sum))
    uex_sum[as.integer(names(r))] <- r
    if (!is.null(cvg2)) {
      v <- Views(cvg2, uex[strands == "-"])
      names(v) <- NULL
      r <- unlist(lapply(v[lengths(v) > 0], sum))
      uex_sum[as.integer(names(r))] <- r
    }
  }

  # Map sums back from unique exons to all exons, then relist by transcript
  ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]
  ex_sum <- uex_sum[ex2uex]                       # one value per exon in 'ex'
  exon_sums_by_tx <- relist(ex_sum, transcripts)  # IntegerList parallel to 'transcripts'
  mcols(exon_sums_by_tx) <- mcols(transcripts)

  if (return_per_exon)
    return(exon_sums_by_tx)

  # Per-transcript totals (sum across each element of the IntegerList)
  tx_totals <- sum(exon_sums_by_tx)
  names(tx_totals) <- names(transcripts)
  return(tx_totals)
}

#' Get coverage from fst large coverage format
#'
#' @param grl a GRangesList
#' @param fst_index a path to an existing fst index file
#' @param columns NULL or character, default NULL. Else must be a subset of
#' names in the fst files. Run ids etc.
#' @return a list, each element is a data.table of coverage
#' @export
#' @examples
#' library(data.table)
#' library(ORFik)
#' grl <- GRangesList("1:1-5:+")
#' tempdir <- tempdir()
#' fst_index <- file.path(tempdir, "coverage_index.fst")
#' mock_run_names <- c("SRR1010101", "SRR1010102", "SRR1010103")
#' coverage_file <- file.path(tempdir, paste0("coverage_1_part1_",
#'  c("forward", "reverse"), ".fst"))
#' mock_coverage <- setnames(setDT(lapply(mock_run_names, function(x) {
#'  sample(seq(0, 100), 100, replace = TRUE, prob = c(0.95, rep(0.01, 100)))})),
#'  mock_run_names)
#' mock_index <- data.table(chr = "1", start = 1, end = nrow(mock_coverage),
#'  file_forward = coverage_file[1], file_reverse = coverage_file[2])
#'
#' fst::write_fst(mock_index, fst_index)
#' fst::write_fst(mock_coverage, coverage_file[1])
#'
#' coverageByTranscriptFST(grl, fst_index)
#' coverageByTranscriptFST(grl, fst_index, c("SRR1010101", "SRR1010102"))
coverageByTranscriptFST <- function(grl, fst_index, columns = NULL) {
  if (!file.exists(fst_index)) stop("No valid fst index file at location: ", fst_index)
  index <- read_fst(fst_index, as.data.table = TRUE)
  if (dirname(index$file_forward[1]) != dirname(fst_index)) {
    index[, file_forward := file.path(dirname(fst_index), basename(file_forward))]
    index[, file_reverse := file.path(dirname(fst_index), basename(file_reverse))]
  }
  if (!file.exists(index$file_forward[1])) {
    stop("Fst index found, but no coverage fst page files found, first file missing: ",
         index$file_forward[1],
         "\n  Index located at: ", fst_index)
  }

  valid_chromosomes <- unique(index$chr)
  input_chromosomes <- as.character(unique(unlist(seqnames(grl), use.names = FALSE)))
  stopifnot(all(input_chromosomes %in% valid_chromosomes))

  # Find page files needed per grl
  tiles_all <- lapply(grl, function(gr) {
    dt <-  index[chr %in% unique(as.character(seqnames(gr))),]
    setkey(dt, start, end)
    query <- as.data.table(ranges(gr))[, c(1,2), with = FALSE]
    setnames(query, c("start", "end"))
    query[, query_id := .I]  # optional: track which query range matched
    setkey(query, start, end)
    res <- foverlaps(dt, query, nomatch = 0)

    res[, start_segment := pmax(start - i.start + 1, 1)]
    res[, end_segment := pmin(end - i.start + 1, i.end - i.start + 1)]
    res[, strand := as.character(strand(gr)[query_id])]
    res[, file := ifelse(strand == "+", file_forward, file_reverse)]
    res[order(query_id),][]
  })

  # Read data from appropriate files using offset
  results <- lapply(tiles_all, function(tiles) {
    fst::threads_fst(5, reset_after_fork = FALSE)
    return(rbindlist(lapply(data.table::transpose(tiles[, .(start_segment, end_segment, file, strand)]), function(x) {
      d <- read_fst(x[3], columns,
                    from = as.numeric(format(x[1], scientific=FALSE)),
                    to = as.numeric(format(x[2], scientific=FALSE)),
                    as.data.table = TRUE)
      if (x[4] == "-") d <- d[rev(seq.int(nrow(d))),]
      return(d)
    })))
  })
  return(results)
}

#' Get names of GRangesList
#'
#' Faster version than S4Vector generic caller
#' @param x a GRangesList
#' @return a character vector
#' @export
setMethod(
  "names",
  signature(x = "GRangesList"),
  function(x) {
    return(x@partitioning@NAMES)
  }
)

#' Get names of GRangesList
#'
#' Faster version than S4Vector generic caller
#' @param x a GRangesList
#' @param value character vector of names
#' @return a GRangesList with updated names
#' @export
setReplaceMethod(
  "names",
  signature(x = "GRangesList"),
  function(x, value) {
    ## basic sanity check
    if (!is.null(value) && length(value) != length(x)) {
      stop("length(names) must equal length of GRangesList")
    }
    x@partitioning@NAMES <- if (!is.null(value)) {
       as.character(value)
    } else value

    return(x)
  }
)

#' Get length of GRangesList
#'
#' Faster version than S4Vector generic caller
#' @param x a GRangesList
#' @return an integer (length 1)
#' @export
setMethod(
  "length",
  signature(x = "GRangesList"),
  function(x) {
    return(length(x@partitioning@end))
  }
)

#' Get widths of GRanges
#'
#' Faster version than S4Vector generic caller
#' @param x a GRanges
#' @return an integer (length equal to x)
#' @export
setMethod(
  "width",
  "GRanges",
  function(x) {
    return(x@ranges@width)
  }
)

#' Get starts of GRanges
#'
#' Faster version than S4Vector generic caller
#' @param x a GRanges
#' @return an integer (length equal to x)
#' @export
setMethod(
  "start",
  "GRanges",
  function(x) {
    return(x@ranges@start)
  }
)

#' Get ends of GRanges
#'
#' Faster version than S4Vector generic caller
#' @param x a GRanges
#' @return an integer (length equal to x)
#' @export
setMethod(
  "end",
  "GRanges",
  function(x) {
    return(x@ranges@width - x@ranges@start + 1)
  }
)

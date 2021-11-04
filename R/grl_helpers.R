#' Restrict GRangesList
#'
#' Will restrict GRangesList to `N` bp downstream from the first base.
#' @param grl (GRangesList)
#' @param firstN (integer) Allow only this many bp downstream, maximum.
#' @return a GRangesList of reads restricted to firstN and tiled by 1
#' @keywords internal
downstreamN <- function(grl, firstN = 150L) {
  return(heads(tile1(grl, matchNaming = FALSE), firstN))
}

#' Get list of widths per granges group
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
  validGRL(class(grl))
  if (keep.names) {
    return(sum(width(grl)))
  } else {
    return(as.integer(sum(width(grl))))
  }
}


#' Get list of seqnames per granges group
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

#' Get list of strands per granges group
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
  if (keep.names) {
    return(heads(strand(grl), 1L))
  } else {
    return(as.character(heads(strand(grl), 1L)))
  }
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

#' Get flanks per group
#'
#' For a GRangesList, get start and end site, return back as GRL.
#' @param grl a \code{\link{GRangesList}}
#' @return a GRangesList, 1 GRanges per group with:
#'  start as minimum start of group and end as maximum per group.
#' @export
#' @examples
#' grl_to_flank <- GRangesList(tx1 = GRanges("1", IRanges(c(1,5), width = 2), "+"),
#'                      tx2 = GRanges("2", IRanges(c(10,15), width = 2), "+"))
#' flankPerGroup(grl)
flankPerGroup <- function(grl) {
  validGRL(class(grl))
  dt <- as.data.table(grl)
  dt$group <- NULL
  dt <- dt[, .(start = min(start), end = max(end), seqnames = seqnames[1],
               strand = strand[1]), by = group_name]
  rownames(dt) <- dt$group_name
  dt$group_name <- NULL
  gr <- GRanges(dt, seqinfo = seqinfo(grl))
  grl <- groupGRangesBy(gr)
  return(grl)
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
#' grl <- GRangesList(tx1_1 = ORF)
#' unlistGrl(grl)
#'
unlistGrl <- function(grl) {
  validGRL(class(grl))
  if (length(grl) == 0) return(unlist(grl))
  unl <- unlist(grl[1], use.names = FALSE)
  return(unlist(grl, use.names = is.null(names(unl))))
}


#' Removes meta columns
#'
#' @param grl a GRangesList or GRanges object
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
#' @importFrom S4Vectors wmsg isTRUEorFALSE
#' @return Integer Rle of coverage, 1 per transcript
coverageByTranscriptW <- function (x, transcripts, ignore.strand = FALSE,
                                    weight = 1L) {
  if (!is(transcripts, "GRangesList")) {
    transcripts <- try(exonsBy(transcripts, by = "tx", use.names = TRUE),
                       silent = TRUE)
    if (is(transcripts, "try-error"))
      stop(wmsg("failed to extract the exon ranges ",
                "from 'transcripts' with ", "exonsBy(transcripts, by=\"tx\", use.names=TRUE)"))
  }
  if (!isTRUEorFALSE(ignore.strand))
    stop(wmsg("'ignore.strand' must be TRUE or FALSE"))
  seqinfo(x) <- GenomicFeatures:::.merge_seqinfo_and_infer_missing_seqlengths(x,
                                                                              transcripts)
  ex <- unlist(transcripts, use.names = FALSE)
  sm <- selfmatch(ex)
  is_unique <- sm == seq_along(sm)
  uex2ex <- which(is_unique)
  uex <- ex[uex2ex]
  # Fix GAlignments not allowing mcol weight, remove when they fix it
  # in GAlignments definition of coverage.
  if ((is(x, "GAlignments") | is(x, "GAlignmentPairs"))
      & is.character(weight)) {
    if (!(weight %in% colnames(mcols(x))))
      stop("weight is character and not mcol of x,",
           " check spelling of weight.")
    weight <- mcols(x)[, weight]
    x <- grglist(x) # convert to grl
    weight = weight[groupings(x)] # repeat weight per group
  }

  if (ignore.strand) {
    cvg <- coverage(x, weight = weight)
    uex_cvg <- cvg[uex]
  }
  else {
    pluss <- BiocGenerics::`%in%`(strand(x), c("+", "*"))
    minus <- BiocGenerics::`%in%`(strand(x), c("-", "*"))
    x1 <- x[pluss]
    x2 <- x[minus]
    if (length(weight) > 1) {
      # Add unlist in case of GAlignments
      cvg1 <- coverage(x1, weight = weight[as.logical(unlist(pluss))])
      cvg2 <- coverage(x2, weight = weight[as.logical(unlist(minus))])
    } else {
      cvg1 <- coverage(x1, weight = weight)
      cvg2 <- coverage(x2, weight = weight)
    }

    is_plus_ex <- strand(uex) == "+"
    is_minus_ex <- strand(uex) == "-"
    if (!identical(is_plus_ex, !is_minus_ex))
      stop(wmsg("'transcripts' has exons on the * strand. ",
                "This is not supported at the moment."))
    uex_cvg <- cvg1[uex]
    uex_cvg[is_minus_ex] <- cvg2[uex[is_minus_ex]]
  }
  uex_cvg[strand(uex) == "-"] <- revElementsF(uex_cvg)[strand(uex) == "-"]
  ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]
  ex_cvg <- uex_cvg[ex2uex]
  ans <- IRanges:::regroupBySupergroup(ex_cvg, transcripts)
  mcols(ans) <- mcols(transcripts)
  return(ans)
}

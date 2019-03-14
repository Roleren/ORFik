#' Restrict GRangesList
#'
#' Will restrict GRangesList to `N` bp downstream from the first base.
#' @param grl (GRangesList)
#' @param firstN (integer) Allow only this many bp downstream, maximum.
#' @return a GRangesList of reads restricted to firstN and tiled by 1
#'
downstreamN <- function(grl, firstN = 150L) {
  return(heads(tile1(grl, matchNaming = FALSE), firstN))
}

#' Get list of widths per granges group
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not
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
#' @param keep.names a boolean, keep names or not
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
    return(seqnames(heads(grl, 1L)))
  } else {
    return(as.character(seqnames(heads(grl, 1L))))
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
#'
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
  # test naming, this is still not perfect
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
#' @param ignore.strand a boolean, if FALSE: should minus strands be
#' sorted from highest to lowest ends. If TRUE: from lowest to highest ends.
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
sortPerGroup <- function(grl, ignore.strand = FALSE){
  if (!ignore.strand) {
    indicesPos <- strandBool(grl)

    grl[indicesPos] <- gSort(grl[indicesPos])
    grl[!indicesPos] <- gSort(grl[!indicesPos], decreasing = TRUE,
                              byStarts = FALSE)
    return(grl)
  }
  return(gSort(grl))
}


#' Get list of strands per granges group
#' @param grl a \code{\link{GRangesList}}
#' @param keep.names a boolean, keep names or not
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
    return(strand(heads(grl, 1L)))
  } else {
    return(as.character(strand(heads(grl, 1L))))
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
#' @param keep.names a boolean, keep names or not
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
#' @param keep.names a boolean, keep names or not
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
#' @param keep.names a boolean, keep names or not
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
#' @param keep.names a boolean, keep names or not
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
#' @param grl a GRangesList
#' @param keep.names a boolean, keep names or not
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
#' grl <- GRangesList(tx1_1 = ORF)
#' unlistGrl(grl)
#'
unlistGrl <- function(grl) {
  validGRL(class(grl))
  unl <- unlist(grl[1], use.names = FALSE)
  return(unlist(grl, use.names = is.null(names(unl))))
}


#' Removes meta columns
#'
#' @param grl a GRangesList or GRanges object
#' @return same type and structure as input without meta columns
#'
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

#' Get number of ranges per group as an iterator
#'
#' @param grl GRangesList
#' @return an integer vector
#' @examples
#' grl <- GRangesList(GRanges("1", c(1, 3, 5), "+"),
#'                    GRanges("1", c(19, 21, 23), "+"))
#' ORFik:::groupings(grl)
#'
groupings <- function(grl){
  l <- lengths(grl, use.names = FALSE)
  return(rep.int(seq.int(length(l)), l))
}

#' Restrict GRangesList
#'
#' Will restrict GRangesList to \code{N} bp downstream from the first base.
#' @param grl (GRangesList)
#' @param firstN (integer) Allow only this many bp downstream
#' @return a GRangesList of reads restricted to firstN and tiled by 1
#'
downstreamN <- function(grl, firstN = 150L) {
  grl <- tile1(grl)
  return(heads(grl, firstN))
}

#' Get list of widths per granges group
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
#' A helper for \code{\link{sortPerGroup}}.
#' A faster, more versatile reimplementation of GenomicRanges::sort()
#' Normally not used directly.
#' Groups first each group, then either decreasing or increasing
#' (on starts if byStarts == T, on ends if byStarts == F)
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
  # test naming
  testName <- names(unlist(grl[1], use.names = FALSE)[1])
  if (is.null(testName)) {
    DT[, group := NULL]
    asgrl <- makeGRangesListFromDataFrame(
      DT, split.field = "group_name",
      keep.extra.columns = TRUE)
  } else {
    if (!any(grep(pattern = "_", testName))) {
      DT[, group := gsub("_[0-9]*", "", DT$group_name)]
    } else {
      DT[, group := DT$group_name]
    }
    asgrl <- makeGRangesListFromDataFrame(
      DT, split.field = "group_name",
      names.field = "group", keep.extra.columns = TRUE)
  }

  asgrl <- asgrl[names(grl)]
  return(asgrl)
}


#' Sort a GRangesList
#'
#' A faster, more versatile reimplementation of
#' \code{\link[GenomicRanges]{sort.GenomicRanges}} for GRangesList,
#' which works poorly for more than 10k groups.
#' This function sorts each group, where "+" strands are
#' increasing by starts and "-" strands are decreasing by ends.
#'
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
lastExonEndPerGroup = function(grl,keep.names = TRUE) {
  validGRL(class(grl))
  if (keep.names) {
    return(end(lastExonPerGroup(grl)))
  } else {
    return(as.integer(end(lastExonPerGroup(grl))))
  }
}


#' Get last start per granges group
#' @param grl a \code{\link[GenomicRanges]{GRangesList}}
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
lastExonStartPerGroup = function(grl, keep.names = TRUE) {
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

  # Get Rle -> to logcal -> sum, that is groups
  exonsPerGroup <- sum(IRanges::LogicalList(strand(grl) == strand(grl)))
  if(!keep.names) {
    names(exonsPerGroup) <- NULL
  }
  return(exonsPerGroup)
}


#' Safe unlist
#'
#' Same as \code{\link[AnnotationDbi]{unlist2}}, keeps names correctly.
#' One difference is that if grl have no names, it will not
#' make integer names, but keep them as null.
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
  g <- character()
  if (!is.gr_or_grl(class(grl))) {
    stop("Can only remove meta columns from GRangesList or GRanges objects")
  }

  if (is.grl(class(grl))) {
    grouping <- numExonsPerGroup(grl)

    for (i in seq_along(grouping)) {
      g <- c(g, rep(names(grouping[i]), grouping[i]))
    }
    grl <- unlistGrl(grl)
    wasGRL <- TRUE
  }
  elementMetadata(grl)  <- S4Vectors::DataFrame(
    matrix(nrow = length(grl), ncol = 0))

  if (wasGRL) {
    return(groupGRangesBy(grl, g))
  } else {
    return(grl)
  }
}

#' Regroup rle from GRangesList
#'
#' Almost direct copy of IRanges regroupBySupergroup.
#' But only works on rle and GRangesList.
#' This function will be removed if
#' IRanges regroupBySupergroup is exported.
#' @param rle A RleList to reduce groups on.
#' @param supergroups A GRangesList to group by
#' @return A regrouped RleList
regroupRleList<- function(rle, supergroups) {
  supergroups <- PartitioningByEnd(supergroups)
  ans_breakpoints <- end(PartitioningByEnd(rle))[end(supergroups)]

  ans_partitioning <- PartitioningByEnd(ans_breakpoints,
                                        names=names(supergroups))

  return(relist(unlist(rle, use.names=FALSE), ans_partitioning))
}

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
    uex_cvg[is_plus_ex] <- cvg1[uex[is_plus_ex]]
    uex_cvg[is_minus_ex] <- cvg2[uex[is_minus_ex]]
    names(uex_cvg) <- as.character(seqnames(uex))
  }
  uex_cvg[strand(uex) == "-"] <- revElementsF(uex_cvg)[strand(uex) == "-"]
  ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]
  ex_cvg <- uex_cvg[ex2uex]
  ans <- IRanges:::regroupBySupergroup(ex_cvg, transcripts)
  mcols(ans) <- mcols(transcripts)
  return(ans)
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
    res[, file := ifelse(strand == "+", file_forward, file_reverse)][]
  })

  # Read data from appropriate files using offset
  results <- lapply(tiles_all, function(tiles) {
    fst::threads_fst(5, reset_after_fork = FALSE)
    return(rbindlist(lapply(data.table::transpose(tiles[, .(start_segment, end_segment, file, strand)]), function(x) {
      d <- read_fst(x[3], columns,
                    from = as.numeric(format(x[1], scientific=FALSE)),
                    to = as.numeric(format(x[2], scientific=FALSE)),
                    as.data.table = TRUE)
      if (x[4] == "-") d <- d[rev(seq.int(nrow(dt))),]
      return(d)
    })))
  })
  return(results)
}

# Testing new version
# coverageByTranscriptW2 <- function (x, transcripts, ignore.strand = FALSE,
#                                    weight = 1L) {
#   if (!is(transcripts, "GRangesList")) {
#     transcripts <- try(exonsBy(transcripts, by = "tx", use.names = TRUE),
#                        silent = TRUE)
#     if (is(transcripts, "try-error"))
#       stop(wmsg("failed to extract the exon ranges ",
#                 "from 'transcripts' with ", "exonsBy(transcripts, by=\"tx\", use.names=TRUE)"))
#   }
#   if (!isTRUEorFALSE(ignore.strand))
#     stop(wmsg("'ignore.strand' must be TRUE or FALSE"))
#   seqinfo(x) <- GenomicFeatures:::.merge_seqinfo_and_infer_missing_seqlengths(x,
#                                                                               transcripts)
#   ex <- unlist(transcripts, use.names = FALSE)
#   sm <- selfmatch(ex)
#   is_unique <- sm == seq_along(sm)
#   uex2ex <- which(is_unique)
#   uex <- ex[uex2ex]
#   # Fix GAlignments not allowing mcol weight, remove when they fix it
#   # in GAlignments definition of coverage.
#   if ((is(x, "GAlignments") | is(x, "GAlignmentPairs"))
#       & is.character(weight)) {
#     if (!(weight %in% colnames(mcols(x))))
#       stop("weight is character and not mcol of x,",
#            " check spelling of weight.")
#     weight <- mcols(x)[, weight]
#     x <- grglist(x) # convert to grl
#     weight = weight[groupings(x)] # repeat weight per group
#   }
#
#   if (ignore.strand) {
#     cvg <- coverage(x, weight = weight)
#     uex_cvg <- cvg[uex]
#   }
#   else {
#     pluss <- BiocGenerics::`%in%`(strand(x), c("+", "*"))
#     minus <- BiocGenerics::`%in%`(strand(x), c("-", "*"))
#     x1 <- x[pluss]
#     x2 <- x[minus]
#     if (length(weight) > 1) {
#       # Add unlist in case of GAlignments
#       cvg1 <- coverage(x1, weight = weight[as.logical(unlist(pluss))])
#       cvg2 <- coverage(x2, weight = weight[as.logical(unlist(minus))])
#     } else {
#       cvg1 <- coverage(x1, weight = weight)
#       cvg2 <- coverage(x2, weight = weight)
#     }
#
#     is_plus_ex <- strand(uex) == "+"
#     is_minus_ex <- strand(uex) == "-"
#     if (!identical(is_plus_ex, !is_minus_ex))
#       stop(wmsg("'transcripts' has exons on the * strand. ",
#                 "This is not supported at the moment."))
#     # uex_cvg <- cvg1[uex]
#     # uex_cvg[is_minus_ex] <- cvg2[uex[is_minus_ex]]
#     uex_cvg <- RleList(rep(IntegerList(1), length(uex)))
#     uex_cvg[is_plus_ex] <- cvg1[uex[is_plus_ex]]
#     uex_cvg[is_minus_ex] <- cvg2[uex[is_minus_ex]]
#     names(uex_cvg) <- as.character(seqnames(uex))
#
#     # uex_irl_pos <- split(ranges(uex[is_plus_ex]), as.character(seqnames(uex[is_plus_ex])))
#     # uex_irl_neg <- split(ranges(uex[is_minus_ex]), as.character(seqnames(uex[is_minus_ex])))
#     # uex_irl_pos <- uex_irl_pos[names(cvg1)[names(cvg1) %in% unique(names(uex_irl_pos))]]
#     # uex_irl_neg <- uex_irl_neg[names(cvg2)[names(cvg2) %in% unique(names(uex_irl_neg))]]
#     # uex_cvg3_pos <- RleViewsList(rleList = cvg1[names(uex_irl_pos)], rangesList = uex_irl_pos)
#     # uex_cvg3_neg <- RleViewsList(rleList = cvg1[names(uex_irl_neg)], rangesList = uex_irl_neg)
#   }
#   uex_cvg[strand(uex) == "-"] <- ORFik:::revElementsF(uex_cvg)[strand(uex) == "-"]
#   ex2uex <- (seq_along(sm) - cumsum(!is_unique))[sm]
#   ex_cvg <- uex_cvg[ex2uex]
#   ans <- IRanges:::regroupBySupergroup(ex_cvg, transcripts)
#   mcols(ans) <- mcols(transcripts)
#   return(ans)
# }

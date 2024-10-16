#' Defines trailers for ORF.
#'
#' Creates GRanges object as a trailer for ORFranges representing ORF,
#' maintaining restrictions of transcriptRanges. Assumes that ORFranges
#' is on the transcriptRanges, strands and seqlevels are in agreement.
#' When lengthOFtrailer is smaller than space left on the transcript than
#' all available space is returned as trailer.
#'
#' It assumes that ORFranges and transcriptRanges are not
#' sorted when on minus strand. Should be like:
#' (200, 600)
#' (50, 100)
#'
#' @param ORFranges GRanges object of your Open Reading Frame.
#' @param transcriptRanges GRanges object of transtript.
#' @param lengthOftrailer Numeric. Default is 10.
#' @return A GRanges object of trailer.
#' @export
#' @import IRanges
#' @import GenomicRanges
#' @importFrom S4Vectors runValue
#' @family ORFHelpers
#' @examples
#' ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
#'                      ranges = IRanges(start = c(1, 10, 20),
#'                                       end = c(5, 15, 25)),
#'                      strand = "+")
#' transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
#'                             ranges = IRanges(start = c(1, 10, 20, 30, 40),
#'                                              end = c(5, 15, 25, 35, 45)),
#'                      strand = "+")
#' defineTrailer(ORFranges, transcriptRanges)
#'
defineTrailer <- function(ORFranges, transcriptRanges, lengthOftrailer = 200) {

  strands <- runValue(strand(ORFranges))
  leftSpace <- setdiff(transcriptRanges, ORFranges)
  if (strands == "-") {
    leftSpace <- leftSpace[end(leftSpace) <= start(ORFranges)[1]]
  } else {
    leftSpace <- leftSpace[start(leftSpace) >= end(ORFranges)[1]]
  }

  if (sum(width(leftSpace)) <= lengthOftrailer) {
    return(leftSpace)
  } else {
    widths <- if (strands == "-") { rev(cumsum(rev(width(leftSpace)))
                          - lengthOftrailer - 1) } else {
                              cumsum(width(leftSpace)) - lengthOftrailer - 1
                            }
    whichExon <- which(widths >= 0)
    whichExon <- whichExon[if (strands == "-") {length(whichExon)} else 1]

    return(c(resize(leftSpace[whichExon],
                    width(leftSpace[whichExon]) - widths[whichExon] - 1,
                    fix = "start"),
             leftSpace[which(widths < 0)]))
  }
}


#' Map orfs to genomic coordinates
#'
#' Creates GRangesList from the results of ORFs_as_List and
#'  the GRangesList used to find the ORFs
#'
#' There is no check on invalid matches, so be carefull if you use this
#' function directly.
#' @param grl A \code{\link{GRangesList}} of the original
#'  sequences that gave the orfs in Genomic coordinates. If grl_is_sorted = TRUE (default),
#'  negative exon ranges per grl object must be sorted in descending orders.
#' @param result IRangesList A list of the results of finding uorfs
#' list syntax is: Per list group in IRangesList is per grl index. In
#' transcript coordinates. The names are grl index as character.
#' @param groupByTx logical (T), should output GRangesList be grouped by
#' transcripts (T) or by ORFs (F)?
#' @param grl_is_sorted logical, default FALSE If FALSE will sort negative transcript
#' in descending order for you. If you loaded ranges with default methods this is
#' already the case, so you can set to TRUE to save some time.
#' @return A \code{\link{GRangesList}} of ORFs.
#' @family ORFHelpers
#' @keywords internal
mapToGRanges <- function(grl, result, groupByTx = TRUE, grl_is_sorted = FALSE) {
  if(length(result) == 0) return(GRangesList())
  validGRL(class(grl))
  stopifnot(is(grl_is_sorted, "logical"))
  if (is.null(names(grl))) stop("'grl' contains no names.")
  if (!is(result, "IRangesList")) stop("Invalid type of result, must be IRL.")
  if (is.null(names(result))) stop("result IRL has no names")
  if (!grl_is_sorted) grl <- sortPerGroup(grl)

  # map transcripts to genomic coordinates, reduce away false hits
  genomicCoordinates <- pmapFromTranscriptF(result, grl, TRUE)

  return(makeORFNames(genomicCoordinates, groupByTx))
}


#' Get transcript names from orf names
#'
#' Using the ORFik definition of orf name, which is:
#' example ENSEMBL:\cr
#' tx name: ENST0909090909090\cr
#' orf id: _1 (the first of on that tx)\cr
#' orf_name: ENST0909090909090_1\cr
#' So therefor txNames("ENST0909090909090_1") = ENST0909090909090\cr
#'
#' The names must be extracted from a column called names, or the names of the
#' grl object. If it is already tx names, it returns the input
#'
#' NOTE! Do not use _123 etc in end of transcript names if it is not ORFs.
#' Else you will get errors. Just _ will work, but if transcripts are called
#' ENST_123124124000 etc, it will crash, so substitute "_" with "."
#' gsub("_", ".", names)
#' @param grl a \code{\link{GRangesList}} grouped by ORF
#' , GRanges object or IRanges object.
#' @param ref a reference \code{\link{GRangesList}}. The object
#' you want grl to subset by names. Add to make sure naming is
#' valid.
#' @param unique a boolean, if true unique the names,
#'  used if several orfs map to same transcript and you only
#'  want the unique groups
#' @export
#' @return a character vector of transcript names,
#'  without _* naming
#' @family ORFHelpers
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1_1 = gr_plus, tx2_1 = gr_minus)
#' # there are 2 orfs, both the first on each transcript
#' txNames(grl)
#'
txNames <- function(grl, ref = NULL, unique = FALSE) {
  reference_can_be_checked <- !is.null(ref) && length(ref) > 0
  if (reference_can_be_checked) {
    trailing_numbers <- length(grep("_\\d+$", names(ref[1]))) > 0
    if (trailing_numbers && (names(grl[1]) %in% names(ref))) {
      # If reference is also ORFik ORF
      if (unique) {
        return(unique(names(grl)))
      } else return(names(grl))
    }
  }

  if (!is.range(grl)) {
    stop("grl must be GRangesList, GRanges, IRangesList or IRanges Object")
  }

  if (is.null(names(grl))) {
    # should the ! be removed here ?
    if (!is.grl(class(grl))) {
      otherPossibility <- unlist(grl, use.names = FALSE)$names
    } else {
      otherPossibility <- grl$names
    }

    if (is.null(otherPossibility)) {
      stop("grl have no valid orf names to convert")
    } else {
      if (unique) {
        return(sub("_\\d+$", "", unique(otherPossibility), perl = TRUE))
      }
      return(sub("_\\d+$", "", otherPossibility, perl = TRUE))
    }
  }
  if (unique) {
    return(sub("_\\d+$", "", unique(names(grl)), perl = TRUE))
  }
  return(sub("_\\d+$", "", names(grl), perl = TRUE))
}


#' Get the start sites from a GRangesList of orfs grouped by orfs
#'
#' In ATGTTTTGG, get the position of the A.
#' @param grl a \code{\link{GRangesList}} object
#' @param asGR a boolean, return as GRanges object
#' @param keep.names a logical (FALSE), keep names of input.
#' @param is.sorted a speedup, if you know the ranges are sorted
#' @return if asGR is False, a vector, if True a GRanges object
#' @export
#' @family ORFHelpers
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' startSites(grl, is.sorted = FALSE)
#'
startSites <- function(grl, asGR = FALSE, keep.names = FALSE,
                       is.sorted = FALSE) {
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
  }
  posIds <- strandBool(grl)

  startSites <- rep(NA, length(grl))
  startSites[posIds] <- firstStartPerGroup(grl[posIds], FALSE)
  startSites[!posIds] <- firstEndPerGroup(grl[!posIds], FALSE)

  if (asGR) {
    startSites <- GRanges(seqnames = seqnamesPerGroup(grl, FALSE),
                          ranges = IRanges(startSites, width = 1),
                          strand = strandPerGroup(grl, FALSE),
                          seqinfo = seqinfo(grl))
  }
  if (keep.names) {
    names(startSites) <- names(grl)
  }
  return(startSites)
}


#' Get the stop sites from a GRangesList of orfs grouped by orfs
#'
#' In ATGTTTTGC, get the position of the C.
#' @inheritParams startSites
#' @return if asGR is False, a vector, if True a GRanges object
#' @export
#' @family ORFHelpers
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' stopSites(grl, is.sorted = FALSE)
#'
stopSites <- function(grl, asGR = FALSE, keep.names = FALSE,
                      is.sorted = FALSE) {
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
  }
  posIds <- strandBool(grl)

  stopSites <- rep(NA, length(grl))
  stopSites[posIds] <- lastExonEndPerGroup(grl[posIds], FALSE)
  stopSites[!posIds] <- lastExonStartPerGroup(grl[!posIds], FALSE)

  if (asGR) {
    stopSites <- GRanges(seqnames = seqnamesPerGroup(grl, FALSE),
                         ranges = IRanges(stopSites, width = 1),
                         strand = strandPerGroup(grl, FALSE),
                         seqinfo = seqinfo(grl))
  }
  if (keep.names) {
    names(stopSites) <- names(grl)
  }
  return(stopSites)
}


#' Get the Start codons(3 bases) from a GRangesList of orfs grouped by orfs
#'
#' In ATGTTTTGA, get the positions ATG.
#' It takes care of exons boundaries, with exons < 3 length.
#' @param grl a \code{\link{GRangesList}} object
#' @param is.sorted a boolean, a speedup if you know the ranges are sorted
#' @return a GRangesList of start codons, since they might be split on exons
#' @export
#' @family ORFHelpers
#' @examples
#' gr_plus <- GRanges(seqnames = "chr1",
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = "+")
#' gr_minus <- GRanges(seqnames = "chr2",
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = "-")
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' startCodons(grl, is.sorted = FALSE)
#'
startCodons <- function(grl, is.sorted = FALSE){
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
  }
  firstExons <- firstExonPerGroup(grl)
  widths <- widthPerGroup(firstExons)
  validWidths <- widths >= 3L
  if (!all(validWidths)) { # fix short exons by tiling
    needToFix <- grl[!validWidths]
    fixedStarts <- reduceKeepAttr(downstreamN(needToFix, 3L))
    grl[!validWidths] <- fixedStarts
  }
  # fix the others the easy way
  firstExons <- firstExons[validWidths]
  posIds <- strandBool(firstExons)

  end(firstExons[posIds]) <- start(firstExons[posIds]) + 2L
  start(firstExons[!posIds]) <- end(firstExons[!posIds]) - 2L

  grl[validWidths] <- firstExons

  return(grl)
}


#' Get the Stop codons (3 bases) from a GRangesList of orfs grouped by orfs
#'
#' In ATGTTTTGA, get the positions TGA.
#' It takes care of exons boundaries, with exons < 3 length.
#' @param grl a \code{\link{GRangesList}} object
#' @param is.sorted a boolean, a speedup if you know the ranges are sorted
#' @return a GRangesList of stop codons, since they might be split on exons
#' @export
#' @family ORFHelpers
#' @examples
#' gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
#'                    ranges = IRanges(c(7, 14), width = 3),
#'                    strand = c("+", "+"))
#' gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
#'                     ranges = IRanges(c(4, 1), c(9, 3)),
#'                     strand = c("-", "-"))
#' grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
#' stopCodons(grl, is.sorted = FALSE)
#'
stopCodons <- function(grl, is.sorted = FALSE) {
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
  }
  lastExons <- lastExonPerGroup(grl)
  widths <- widthPerGroup(lastExons)
  validWidths <- widths >= 3L
  if (!all(validWidths)) { # fix short exons by tiling
    needToFix <- grl[!validWidths]
    tileBy1 <- tile1(needToFix, matchNaming = FALSE)
    fixedStops <- reduceKeepAttr(tails(tileBy1, 3L))
    grl[!validWidths] <- fixedStops
  }
  # fix the others the easy way
  lastExons <- lastExons[validWidths]
  posIds <- strandBool(lastExons)

  start(lastExons[posIds]) <- end(lastExons[posIds]) - 2L
  end(lastExons[!posIds]) <- start(lastExons[!posIds]) + 2L

  grl[validWidths] <- lastExons

  return(grl)
}

#' Start region as GRangesList
#'
#' Get the start region of each ORF. If you want the start codon only,
#' set upstream = 0 or just use \code{\link{startCodons}}.
#' Standard is 2 upstream and 2 downstream, a width 5 window centered at
#' start site. since p-shifting is not 100% accurate, this window is
#' usually the reads from the start site.
#'
#' If tx is null, then upstream will be forced to 0 and downstream to
#' a maximum of grl width (3' UTR end for mRNAs).
#' Since there is no reference for splicing.
#' @param grl a \code{\link{GRangesList}} object
#'  with usually either leaders, cds', 3' utrs or ORFs
#' @param tx default NULL, a GRangesList of transcripts or (container region),
#' names of tx must contain all grl names. The names of grl can also be the
#' ORFik orf names. that is "txName_id"
#' @param is.sorted logical (TRUE), is grl sorted.
#' @param upstream an integer (2), relative region to get upstream from.
#' @param downstream an integer (2), relative region to get downstream from
#' @family features
#' @return a GRanges, or GRangesList object if any group had > 1 exon.
#' @export
#' @examples
#' ## ORF start region
#' orf <- GRangesList(tx1 = GRanges("1", 200:300, "+"))
#' tx <- GRangesList(tx1 = GRanges("1",
#'                    IRanges(c(100, 200), c(195, 400)), "+"))
#' startRegion(orf, tx, upstream = 6, downstream = 6)
#' ## 2nd codon of ORF
#' startRegion(orf, tx, upstream = -3, downstream = 6)
startRegion <- function(grl, tx = NULL, is.sorted = TRUE,
                        upstream = 2L, downstream = 2L) {
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
    is.sorted <- TRUE
  }
  if (is.null(tx)) {
    tx <- grl
  }
  region <- windowPerGroup(startSites(grl, TRUE, TRUE, is.sorted), tx,
                           upstream, downstream)
  return(region)
}

#' Stop region as GRangesList
#'
#' Get the stop region of each ORF / region. If you want the stop codon only,
#' set downstream = 0 or just use \code{\link{stopCodons}}.
#' Standard is 2 upstream and 2 downstream, a width 5 window centered at
#' stop site.
#'
#' If tx is null, then downstream will be forced to 0 and
#' upstream to a minimum of -grl width (to the TSS).
#' .
#' Since there is no reference for splicing.
#' @param grl a \code{\link{GRangesList}} object
#'  with usually either leaders, cds', 3' utrs or ORFs
#' @param tx default NULL, a GRangesList of transcripts or (container region),
#' names of tx must contain all grl names. The names of grl can also be the
#' ORFik orf names. that is "txName_id"
#' @param is.sorted logical (TRUE), is grl sorted.
#' @param upstream an integer (2), relative region to get upstream from.
#' @param downstream an integer (2), relative region to get downstream from
#' @family features
#' @return a GRanges, or GRangesList object if any group had > 1 exon.
#' @export
#' @examples
#' ## ORF stop region
#' orf <- GRangesList(tx1 = GRanges("1", 200:300, "+"))
#' tx <- GRangesList(tx1 = GRanges("1",
#'                    IRanges(c(100, 305), c(300, 400)), "+"))
#' stopRegion(orf, tx, upstream = 6, downstream = 6)
#' ## 2nd last codon of ORF
#' stopRegion(orf, tx, upstream = 6, downstream = -3)
stopRegion <- function(grl, tx = NULL, is.sorted = TRUE,
                        upstream = 2L, downstream = 2L) {
  if (!is.sorted) {
    grl <- sortPerGroup(grl)
    is.sorted <- TRUE
  }
  if (is.null(tx)) {
    tx <- grl
  }
  region <- windowPerGroup(stopSites(grl, TRUE, TRUE, is.sorted), tx,
                           upstream, downstream)
  return(region)
}


#' Get id's for each orf
#'
#' These id's can be uniqued by isoform etc,
#' this is not supported by GenomicRanges.
#' @param grl a \code{\link{GRangesList}}
#' @param with.tx a boolean, include transcript names,
#'  if you want unique orfs, so that they dont have duplicates
#'  from different isoforms, set it to FALSE.
#' @return a character vector of ids, 1 per orf
#' @family ORFHelpers
#' @keywords internal
orfID <- function(grl, with.tx = FALSE) {
  seqnames <- seqnamesPerGroup(grl, FALSE)
  strands <- strandPerGroup(grl, FALSE)

  starts <- start(grl); names(starts) <- NULL
  widths <- width(grl); names(widths) <- NULL
  exonInfo <- paste(starts, widths)
  exonInfo <- paste(exonInfo, sep = '', collapse = ';')

  uorfID <- paste(seqnames, strands, exonInfo, sep = ",")
  if (with.tx) {
    uorfID <- paste(uorfID, txNames(grl))
  }
  return(uorfID)
}


#' Get the unique set of groups in a GRangesList
#'
#' Sometimes \code{\link{GRangesList}} groups might be identical,
#' for example ORFs from different isoforms can have identical ranges.
#' Use this function to reduce these groups to unique elements
#' in \code{\link{GRangesList}} \code{grl}, without names and metacolumns.
#' @param grl a \code{\link{GRangesList}}
#' @return a GRangesList of unique orfs
#' @export
#' @family ORFHelpers
#' @examples
#' gr1 <- GRanges("1", IRanges(1,10), "+")
#' gr2 <- GRanges("1", IRanges(20, 30), "+")
#' # make a grl with duplicated ORFs (gr1 twice)
#' grl <- GRangesList(tx1_1 = gr1, tx2_1 = gr2, tx3_1 = gr1)
#' uniqueGroups(grl)
#'
uniqueGroups <- function(grl) {
  if (length(grl) == 0) return(grl)
  ids <- orfID(grl)
  grl <- grl[!duplicated(ids)]
  gr <- unlist(grl, use.names = FALSE)
  names(gr) <- NULL
  gr$names <- NULL
  grl <- relist(gr, grl)
  names(grl) <- seq.int(1, length(grl))
  return(grl)
}

#' Get unique ordering for GRangesList groups
#'
#' This function can be used to calculate unique numerical identifiers
#' for each of the \code{\link{GRangesList}} elements. Elements of
#' \code{\link{GRangesList}} are unique when the \code{\link{GRanges}}
#' inside are not duplicated, so ranges differences matter as well as
#' sorting of the ranges.
#' @seealso uniqueGroups
#' @param grl a \code{\link{GRangesList}}
#' @return an integer vector of indices of unique groups
#' @export
#' @family ORFHelpers
#' @examples
#' gr1 <- GRanges("1", IRanges(1,10), "+")
#' gr2 <- GRanges("1", IRanges(20, 30), "+")
#' # make a grl with duplicated ORFs (gr1 twice)
#' grl <- GRangesList(tx1_1 = gr1, tx2_1 = gr2, tx3_1 = gr1)
#' uniqueOrder(grl) # remember ordering
#'
#' # example on unique ORFs
#' uniqueORFs <- uniqueGroups(grl)
#' # now the orfs are unique, let's map back to original set:
#' reMappedGrl <- uniqueORFs[uniqueOrder(grl)]
uniqueOrder <- function(grl) {
  ids <- orfID(grl)

  sortedOrder <- data.table::chgroup(ids)
  orderedIDs <- ids[sortedOrder]
  l <- S4Vectors::Rle(orderedIDs)
  grouping <- unlist(lapply(seq.int(nrun(l)), function(x) {
    rep(x, runLength(l)[x])
  }))
  reOrdering <- grouping[order(sortedOrder)]
  return(reOrdering)
}

#' Get longest ORF per stop site
#'
#' Rule: if seqname, strand and stop site is equal, take longest one.
#' Else keep.
#' If IRangesList or IRanges, seqnames are groups, if GRanges or GRangesList
#' seqnames are the seqlevels (e.g. chromosomes/transcripts)
#'
#' @param grl a \code{\link{GRangesList}}/IRangesList, GRanges/IRanges of ORFs
#' @return a \code{\link{GRangesList}}/IRangesList, GRanges/IRanges
#' (same as input)
#' @export
#' @importFrom data.table .I
#' @family ORFHelpers
#' @examples
#' ORF1 = GRanges("1", IRanges(10,21), "+")
#' ORF2 = GRanges("1", IRanges(1,21), "+") # <- longest
#' grl <- GRangesList(ORF1 = ORF1, ORF2 = ORF2)
#' longestORFs(grl) # get only longest
longestORFs <- function(grl) {
  if(length(grl) == 0) return(grl) # if empty

  if (is.grl(class(grl))) { # only for GRangesList
    stops <- stopSites(grl, is.sorted = TRUE)
    widths <- widthPerGroup(grl, FALSE)
    seqnames <- seqnamesPerGroup(grl, FALSE)
    strands <- strandPerGroup(grl, FALSE)
  } else { # GRanges, IRanges or IRangesList
    stops <- unlist(end(grl), use.names = FALSE)
    widths <- unlist(width(grl), use.names = FALSE)

    if (is.gr_or_grl(class(grl))) { # GRanges
      seqnames <- as.character(seqnames(grl))
      strands <- as.character(strand(grl))
      stops[strands == "-"] <- start(grl)[strands == "-"]
    } else { # IRanges or IRangesList
      strands <- rep("+", length(widths))
      if (is(grl, "IRanges")) {
        seqnames <- rep.int(1, length(widths))
      } else if (is(grl, "IRangesList")) {
        seqnames <- rep.int(seq.int(length(grl)), lengths(grl))
      }
    }
  }
  dt <- data.table(seqnames, strands, stops, widths)
  longestORFs <- dt[, .I[which.max(widths)],
                    by = .(seqnames, strands, stops)]$V1
  if (is(grl, "IRangesList")) {
    ir <- unlist(grl, use.names = FALSE)
    ir <- ir[longestORFs]
    irl <- split(ir, seqnames[longestORFs])
    names(irl) <- names(grl)
    return(irl)
  }
  return(grl[longestORFs])
}

#' Get number of codons
#'
#' Length of object / 3.
#' Choose either only whole codons, or with stubs.
#' ORF stubs are not relevant, since there are no correctly defined
#' ORFs that are 17 bases long etc.
#' @param grl a \code{\link{GRangesList}} object
#' @param as.integer a logical (TRUE), remove stub codons
#' @param keep.names a logical (FALSE)
#' @return an integer vector
#' @keywords internal
numCodons <- function(grl, as.integer = TRUE,  keep.names = FALSE) {
  if (as.integer) return(ceiling(widthPerGroup(grl, keep.names) / 3))
  return(widthPerGroup(grl, keep.names) / 3)
}

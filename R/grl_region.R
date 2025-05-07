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

#' Get window region of GRanges object
#'
#' Per GRanges input (gr) of single position inputs (center point),
#' create a GRangesList window output of specified
#' upstream, downstream region relative to some transcript "tx". \cr
#' If downstream is 20, it means the window will start 20 downstream of
#' gr start site (-20 in relative transcript coordinates.)
#' If upstream is 20, it means the window will start 20 upstream of
#' gr start site (+20 in relative transcript coordinates.)
#' It will keep exon structure of tx, so if -20 is on next exon, it jumps to
#' next exon.
#'
#' If a region has a part that goes out of bounds, E.g if you try to get window
#' around the CDS start site, goes longer than the 5' leader start site,
#' it will set start to the edge boundary
#' (the TSS of the transcript in this case).
#' If region has no hit in bound, a width 0 GRanges object is returned.
#' This is useful for things like countOverlaps, since 0 hits will then always
#' be returned for the correct object index. If you don't want the 0 width
#' windows, use \code{reduce()} to remove 0-width windows.
#' @param gr a GRanges/IRanges object (startSites or others,
#'  must be single point per in genomic coordinates)
#' @param tx a \code{\link{GRangesList}} of transcripts or (container region),
#' names of tx must contain all gr names. The names of gr can also be the
#' ORFik orf names. that is "txName_id".
#' @param upstream an integer, default (0), relative region to get
#'  upstream end from. (0 means start site, +1 is one upstream, -1 is one downstream)
#' @param downstream an integer, default (0), relative region to get
#'  downstream end from (0 means start site, +1 is one downstream, -1 is one upstream)
#' @return a GRanges, or GRangesList object if any group had > 1 exon.
#' @export
#' @family ExtendGenomicRanges
#' @importFrom data.table chmatch
#' @examples
#' # find 2nd codon of an ORF on a spliced transcript
#' ORF <- GRanges("1", c(3), "+") # start site
#' names(ORF) <- "tx1_1" # ORF 1 on tx1
#' tx <- GRangesList(tx1 = GRanges("1", c(1,3,5,7,9,11,13), "+"))
#' windowPerGroup(ORF, tx, upstream = 0, downstream = 0) # <- TIS
#' windowPerGroup(ORF, tx, upstream = 0, downstream = 1) # <- first and second base
#' windowPerGroup(ORF, tx, upstream = -1, downstream = 1) # <- second base
#' # find 2nd codon of an ORF on a spliced transcript
#' windowPerGroup(ORF, tx, upstream = -3, downstream = 5) # <- 2nd codon
#'
#' # With multiple extensions downstream
#' ORF <- rep(ORF, 2)
#' names(ORF)[2] <- "tx1_2"
#' windowPerGroup(ORF, tx, upstream = 0, downstream = c(2, 5))
#' # The last one gives 2nd for first ORF and (1st and 2nd) codon for
#' # second ORF, returned as two groups of class GRanges/GRangsList
#'
windowPerGroup <- function(gr, tx, upstream = 0L, downstream = 0L) {
  g <- asTX(gr, tx, tx.is.sorted = TRUE)
  indices <- chmatch(txNames(gr, tx), names(tx))
  if (anyNA(indices)) {
    indices <- chmatch(names(gr), names(tx))
    if (anyNA(indices)) {
      stop("not all references are present, so can not map to transcripts.",
           " Does your annotation contains  '_'+number naming?",
           "ORFik uses this for a special purpose.")
    }
  }
  txEnds <- widthPerGroup(tx[indices], FALSE)
  starts <- pmin(pmax(start(g) - upstream, 1L), txEnds)

  g <- ranges(g)
  if (any(downstream != 0L)) {
    ends <- pmin(pmax(end(g) + downstream, starts - 1), txEnds)
    if (is(g, "IRanges")) # Remove this ?
      g <- IRanges(starts, ends)
  } else {
    start(g) <- starts
  }
  names(g) <- indices
  region <- pmapFromTranscriptF(g, tx, TRUE)
  names(region) <- names(gr)
  return(region)
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

#' Get flanks per group
#'
#' For a GRangesList, get start and end site, return back as GRangesList.
#' @param grl a \code{\link{GRangesList}}
#' @return a GRangesList, 1 GRanges per group with:
#'  start as minimum start of group and end as maximum per group.
#' @export
#' @examples
#' grl <- GRangesList(tx1 = GRanges("1", IRanges(c(1,5), width = 2), "+"),
#'                    tx2 = GRanges("2", IRanges(c(10,15), width = 2), "+"))
#' flankPerGroup(grl)
flankPerGroup <- function(grl) {

  gr <- unlistToExtremities(grl)
  grl <- groupGRangesBy(gr, gr$group)
  names(grl) <- grl@unlistData$names
  grl@unlistData$names <- grl@unlistData$group <- NULL
  return(grl)
}

#' Get flanks as GRanges
#'
#' For a GRangesList, get start and end site, return back as GRanges.
#' @inheritParams flankPerGroup
#' @return a GRangs object with meta column "group", which gives
unlistToExtremities <- function(grl) {
  if (length(grl) == 0) return(unlist(grl))
  validGRL(class(grl))
  old_names <- names(grl)
  grl@unlistData$group <- NULL
  dt <- as.data.table(grl)
  dt[, group_name := NULL]
  dt <- dt[, .(start = min(start), end = max(end), seqnames = seqnames[1],
               strand = strand[1]), by = group]

  gr <- GRanges(seqnames = dt$seqnames, ranges = IRanges(dt$start, dt$end), strand = dt$strand,
                group = dt$group, seqinfo = seqinfo(grl))
  if (!is.null(old_names)) gr$names <- old_names
  names(gr) <- names(grl)
  return(gr)
}

#' Get pseudo introns per Group
#'
#' If an intron is of length < 'width' * 2, it will not be split into pseudo.
#' @param grl a GrangesList of length 1
#' @param width numeric, default 100. The size of pseudo flanks.
#' @return a GRangesList
#' @export
#' @examples
#' tx <- GRangesList(GRanges("1", IRanges(c(1, 150, 1e5, 1e6)), "+"))
#' pseudoIntronsPerGroup(tx) # See intron 1 is not split
#' tx_2 <- rep(GRangesList(GRanges("1", IRanges(c(1, 150, 1e5, 1e6)), "+")), 2)
#' pseudoIntronsPerGroup(tx_2)
#' pseudoIntronsPerGroup(tx_2, 1e6)
pseudoIntronsPerGroup <- function(grl, width = 100) {
  validGRL(class(grl))
  if (length(grl) == 0) return(grl)
  original_names <- names(grl)
  flanks <- flankPerGroup(grl)
  introns <- setdiff(flanks, grl)

  names(introns) <- seq(length(introns))
  gr <- unlistGrl(introns)
  # Only consider introns longer than 1 bp
  gr_is_too_short <- width(gr) < width*2
  gr_too_short <- gr[gr_is_too_short]
  gr_too_short$type <- rep("original", length(gr_too_short))
  gr_too_short$flank <- rep("whole", length(gr_too_short))
  gr <- gr[!gr_is_too_short]
  left <- right <- GRanges()
  if (length(gr) > 0) {
    # Get left flank
    left <- resize(gr, width = width, fix = "start")
    # Get right flank
    right <- resize(gr, width = width, fix = "end")

    # Trim in case flank is longer than intron
    left <- pintersect(left, gr)
    right <- pintersect(right, gr)
    left$hit <- right$hit <- NULL
    left$type <- right$type <- "pseudo"
    left$flank <- "left"
    right$flank <- "right"
  }

  # Combine
  grl_new <- sortPerGroup(groupGRangesBy(c(gr_too_short, left, right)))

  names(grl_new) <- original_names
  return(grl_new)
}

#' Get exons with pseudo introns per Group
#'
#' If an intron is of length < 'width' * 2, it will not be split into pseudo.
#' @param grl a GrangesList of length 1
#' @param width numeric, default 100. The size of pseudo flanks.
#' @return a GRangesList
#' @export
#' @examples
#' tx <- GRangesList(GRanges("1", IRanges(c(1, 150, 1e5, 1e6)), "+"))
#' exonsWithPseudoIntronsPerGroup(tx) # See intron 1 is not split
#' tx_2 <- rep(GRangesList(GRanges("1", IRanges(c(1, 150, 1e5, 1e6)), "+")), 2)
#' exonsWithPseudoIntronsPerGroup(tx_2)
#' tx_3 <- tx_2
#' names(tx_3) <- c("tx1", "tx2")
#' exonsWithPseudoIntronsPerGroup(tx_3, 1e6)
exonsWithPseudoIntronsPerGroup <- function(grl, width = 100) {
  if (length(grl) == 0 || max(lengths(grl)) == 1) return(grl)
  original_names <- names(grl)
  names(grl) <- seq_along(grl)
  pseudo_introns <- ORFik:::pseudoIntronsPerGroup(grl, width)
  exons_introns_gr <- c(unlistGrl(grl), unlistGrl(pseudo_introns))
  mcols(exons_introns_gr) <- NULL
  grl <- reduceKeepAttr(groupGRangesBy(exons_introns_gr), keep.names = TRUE)
  names(grl) <- original_names
  return(grl)
}

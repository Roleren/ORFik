
#' Extends leaders downstream
#'
#' When finding uORFs, often you want to allow them to end inside the cds.
#'
#' This is a simple way to do that
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param cds If you want to extend 5' leaders downstream,
#'  to catch uorfs going into cds, include it.
#' @param onlyFirstExon logical (F), include whole cds or only first exons.
#' @importFrom S4Vectors pc
#' @return a GRangesList of cds exons added to ends
#' @family uorfs
#'
addCdsOnLeaderEnds <- function(fiveUTRs, cds, onlyFirstExon = FALSE) {
  if (length(cds) == 0) {
    warning("cds is empty, returning without using it.")
    return(fiveUTRs)
  }
  if (is.null(names(cds))) {
    warning("cds have no names, returning without using it.")
    return(fiveUTRs)
  }
  matchingNames <- names(fiveUTRs) %in% names(cds)
  areValidNames <- (sum(matchingNames) - length(names(fiveUTRs))) != 0
  if (areValidNames) {
    warning("not all cds names matches fiveUTRs names,
            returning without using cds.")
    return(fiveUTRs)
  }
  ## get only the ones we need
  ## select first in every, they must be sorted!
  if (onlyFirstExon) {
    firstExons <- firstExonPerGroup(cds[names(fiveUTRs)])
    gr <- unlist(firstExons, use.names = FALSE)
    mcols(gr) <- as.data.frame(mcols(unlist(
      fiveUTRs, use.names = FALSE)))[seq_along(gr),]
    grl <- relist(gr, firstExons)
  } else {
    gr <- unlist(cds[names(fiveUTRs)], use.names = FALSE)
    mcols(gr) <- as.data.frame(mcols(unlist(
      fiveUTRs, use.names = FALSE)))[seq_along(gr),]
    grl <- relist(gr, cds[names(fiveUTRs)])
  }
  ## fix mcols of cds, so that pc() will work
  fiveUTRsWithCdsExons <- pc(fiveUTRs, grl)
  ## should we use reduceKeepAttr here ?, we will lose
  ## exon_id if not.
  return(reduce(fiveUTRsWithCdsExons))
}

#' Create search space to look for uORFs
#'
#' Given a GRangesList of 5' UTRs or transcripts, reassign the start
#' sites using max peaks from CageSeq data (if CAGE is given).
#' A max peak is defined as new
#' TSS if it is within boundary of 5' leader range, specified by
#' `extension` in bp. A max peak must also be higher than minimum
#' CageSeq peak cutoff specified in `filterValue`. The new TSS will then
#' be the positioned where the cage read (with highest read count in the
#' interval). If you want to include uORFs going into the CDS, add this
#' argument too.
#'
#' @inheritParams reassignTSSbyCage
#' @param cds (GRangesList) CDS of relative fiveUTRs, applicable only if you
#' want to extend 5' leaders downstream of CDS's, to allow upstream ORFs that
#' can overlap into CDS's.
#' @return a GRangesList of newly assigned TSS for fiveUTRs,
#'  using CageSeq data.
#' @export
#' @family uorfs
#' @examples
#' # example 5' leader, notice exon_rank column
#' fiveUTRs <- GenomicRanges::GRangesList(
#'   GenomicRanges::GRanges(seqnames = "chr1",
#'                          ranges = IRanges::IRanges(1000, 2000),
#'                          strand = "+",
#'                          exon_rank = 1))
#' names(fiveUTRs) <- "tx1"
#'
#' # make fake CageSeq data from promoter of 5' leaders, notice score column
#' cage <- GenomicRanges::GRanges(
#'   seqnames = "chr1",
#'   ranges =  IRanges::IRanges(500, 510),
#'   strand = "+",
#'   score = 10)
#'
#' # finally reassign TSS for fiveUTRs
#' uORFSearchSpace(fiveUTRs, cage)
#'
uORFSearchSpace <- function(fiveUTRs, cage = NULL, extension = 1000,
                            filterValue = 1, restrictUpstreamToTx = FALSE,
                            removeUnused = FALSE, cds = NULL) {
  if (!is.null(cage)) {
    fiveUTRs <- reassignTSSbyCage(fiveUTRs, cage, extension, filterValue,
                                  restrictUpstreamToTx, removeUnused)
  }
  if(!is.null(cds))
    fiveUTRs <- addCdsOnLeaderEnds(fiveUTRs, get("cds", mode = "S4"))
  return(fiveUTRs)
}

#' Remove uORFs that are false CDS hits
#'
#' This is a strong filtering, so that even if the cds is on another transcript
#' , the uORF is filtered out, this is because there is no way of knowing by
#' current ribo-seq, rna-seq experiments.
#'
#' @param uorfs (GRangesList), the uORFs to filter
#' @param cds (GRangesList), the coding sequences (main ORFs on transcripts),
#' to filter against.
#' @return (GRangesList) of filtered uORFs
#' @family uorfs
#'
filterUORFs <- function(uorfs, cds) {
  if (is.null(cds) || length(uorfs) == 0) return(uorfs)
  validGRL(class(cds), "cds")
  message("starting to filter out ourfs...")

  #check start upstream, ending downstream
  uorfs <- removeORFsWithinCDS(uorfs, cds)
  uorfs <- sortPerGroup(uorfs)
  uorfs <- removeORFsWithSameStopAsCDS(uorfs, cds)
  uorfs <- removeORFsWithSameStartAsCDS(uorfs, cds)
  uorfs <- removeORFsWithStartInsideCDS(uorfs, cds)
  message("finished filtering ourfs")

  return(uorfs)
}

#' Remove ORFs that are within cds
#'
#' @param grl (GRangesList), the ORFs to filter
#' @param cds (GRangesList), the coding sequences (main ORFs on transcripts),
#' to filter against.
#' @return (GRangesList) of filtered uORFs
#' @family uorfs
removeORFsWithinCDS <- function(grl, cds) {
  overlaps <- findOverlaps(query = grl, cds, type = "within")
  if (length(overlaps) > 0) return(grl[-unique(from(overlaps))])
  return(grl)
}

#' Remove ORFs that have same stop site as the CDS
#'
#' @inheritParams removeORFsWithinCDS
#' @return (GRangesList) of filtered uORFs
#' @family uorfs
removeORFsWithSameStopAsCDS <- function(grl, cds) {
  overlaps <- findOverlaps(query =  stopSites(grl, asGR = TRUE,
                                              is.sorted = TRUE),
                           stopSites(cds, asGR = TRUE, is.sorted = TRUE),
                           type = "within")
  if (length(overlaps) > 0) return(grl[-unique(from(overlaps))])
  return(grl)
}

#' Remove ORFs that have same start site as the CDS
#'
#' @inheritParams removeORFsWithinCDS
#' @return (GRangesList) of filtered uORFs
#' @family uorfs
removeORFsWithSameStartAsCDS <- function(grl, cds) {
  starts <- startSites(grl, asGR = TRUE, is.sorted = TRUE)
  cdsstarts <- startSites(cds, asGR = TRUE, is.sorted = TRUE)
  overlaps <- findOverlaps(starts, cdsstarts, type = "within")
  if (length(overlaps) > 0) return(grl[-unique(from(overlaps))])
  return(grl)
}

#' Remove ORFs that have start site within the CDS
#'
#' @inheritParams removeORFsWithinCDS
#' @return (GRangesList) of filtered uORFs
#' @family uorfs
removeORFsWithStartInsideCDS <- function(grl, cds) {
  starts <- startSites(grl, asGR = TRUE, is.sorted = TRUE)
  overlaps <- findOverlaps(starts, cds, type = "within")
  if (length(overlaps) > 0) return(grl[-unique(from(overlaps))])
  return(grl)
}

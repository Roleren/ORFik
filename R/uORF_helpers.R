
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
  if(onlyFirstExon) {
    firstExons <- firstExonPerGroup(cds[names(fiveUTRs)])
    gr <- unlist(firstExons, use.names = FALSE)
  } else {
    gr <- unlist(cds[names(fiveUTRs)], use.names = FALSE)
  }

  ## fix mcols of cds, so that pc() will work
  mcols(gr) <- as.data.frame(mcols(unlist(
    fiveUTRs, use.names = FALSE)))[seq_along(gr),]
  grl <- relist(gr, cds[names(fiveUTRs)])
  fiveUTRsWithCdsExons <- pc(fiveUTRs, grl)
  ## should we use reduceKeepAttr here ?, we will lose
  ## exon_id if not.
  return(reduce(fiveUTRsWithCdsExons))
}

#' Create search space to look for uORFs
#'
#' Given a GRangesList of 5' UTRs or transcripts, reassign the start
#' sites using max peaks from CageSeq data. A max peak is defined as new
#' TSS if it is within boundary of 5' leader range, specified by
#' `extension` in bp. A max peak must also be higher than minimum
#' CageSeq peak cutoff specified in `filterValue`. The new TSS will then
#' be the positioned where the cage read (with highest read count in the
#' interval). If you want to include uORFs going into the CDS, add this
#' argument too.
#'
#' @param fiveUTRs (GRangesList) The 5' leaders or transcript sequences
#' @param cage Either a filePath for CageSeq file, or already loaded
#' CageSeq peak data as GRanges.
#' @param extension The maximum number of basses upstream of the TSS to search
#' for CageSeq peak.
#' @param filterValue The minimum number of reads on cage position,
#' for it to be counted as possible new tss.
#' (represented in score column in CageSeq data)
#' If you already filtered, set it to 0.
#' @param cds (GRangesList) CDS of relative fiveUTRs, applicable only if you
#' want to extend 5' leaders downstream of CDS's, to allow upstream ORFs that
#' can overlap into CDS's.
#' @return a GRangesList of newly assigned TSS for fiveUTRs,
#'  using CageSeq data.
#' @export
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
uORFSearchSpace <- function(fiveUTRs, cage, extension = 1000,
                            filterValue = 1, cds = NULL){
  searchSpace <- reassignTSSbyCage(fiveUTRs, cage, extension, filterValue)

  if(!is.null(cds)) searchSpace <- addCdsOnLeaderEnds(searchSpace, cds)

  return(searchSpace)
}


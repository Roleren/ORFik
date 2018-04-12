
#' Filter peak of cage-data by value
#' @param rawCage The raw cage-data
#' @param filterValue The number of counts(score) to filter on
#'  for a tss to pass as hit
#' @return the filtered Granges object
#'
filterCage <- function(rawCage, filterValue = 1) {
  if (filterValue == 0) {
    return(rawCage)
  }
  if (is.null(rawCage$score)) stop("Found no 'score' column in the CageSeq.")
  filteredCage <- rawCage[rawCage$score > filterValue, ] #filter on score
  return(filteredCage)
}

#' Match seqnames
#'
#' Check that seqlevels of fiveUTRs and cage uses
#' the same standard, i.g chr1 vs 1.
#' @param filteredCage Cage-data to check seqnames in
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @return filteredCage with matched seqnames convention
#'
matchSeqlevels <- function(filteredCage, fiveUTRs) {
  fiveSeqlevels <- seqlevels(unlist(fiveUTRs, use.names = FALSE))
  cageSeqlevels <- seqlevels(filteredCage)
  if (length(grep(pattern = "chr", fiveSeqlevels)) > 0 &&
      length(grep(pattern = "chr", cageSeqlevels)) == 0) {
    message("seqnames use different chromosome naming conventions,",
            " trying to fix them")
    # chr1, chr2, not chrX, chrY etc. ->
    regexNormalChr <- '(^[a-zA-Z])*([0-9]+)'
    normalChr <- paste0("chr", grep(
      regexNormalChr, cageSeqlevels, value = TRUE))
    normalChrInd <- grep(regexNormalChr, cageSeqlevels)

    for(i in normalChrInd) {
      seqlevels(filteredCage)[i] <- sub(
        regexNormalChr, normalChr[i], cageSeqlevels[i])
    }
  }
  if (length(grep("chrY", fiveSeqlevels)) == 0 &&
      length(grep("chrY", cageSeqlevels)) != 0)
    seqlevels(filteredCage) <- sub("chrY", "Y", seqlevels(filteredCage))
  if (length(grep("chrX", fiveSeqlevels)) == 0 &&
      length(grep("chrX", cageSeqlevels)) != 0)
    seqlevels(filteredCage) <- sub("chrX", "X", seqlevels(filteredCage))
  if (length(grep("chrM", fiveSeqlevels)) == 0 &&
      length(grep("chrM", cageSeqlevels)) != 0)
    seqlevels(filteredCage) <- sub("chrM", "MT", seqlevels(filteredCage))
  return(filteredCage)
}


#' Extends leaders downstream
#'
#' When reassigning Transcript start sites,
#' often you want to add downstream too.
#' This is a simple way to do that
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param cds If you want to extend 5' leaders downstream,
#'  to catch uorfs going into cds, include it.
#' @importFrom S4Vectors pc
#' @return a GRangesList of cds exons added to ends
#'
addFirstCdsOnLeaderEnds <- function(fiveUTRs, cds) {
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
  firstExons <- firstExonPerGroup(cds[names(fiveUTRs)])
  gr <- unlist(firstExons, use.names = FALSE)
  ## fix mcols of cds, so that pc() will work
  mcols(gr) <- as.data.frame(mcols(unlist(
    fiveUTRs, use.names = FALSE)))[seq_along(gr),]
  grl <- relist(gr, firstExons)
  fiveUTRsWithCdsExons <- pc(fiveUTRs, grl)
  ## should we use reduceKeepAttr here ?, we will lose
  ## exon_id if not.
  return(reduce(fiveUTRsWithCdsExons))
}


#' Extend first exon of each transcript with length specified
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param extension The number of basses upstream to add on transcripts
#' @return granges object of first exons
#'
extendsTSSexons <- function(fiveUTRs, extension = 1000) {
  fiveAsgr <- unlist(fiveUTRs, use.names = TRUE)
  if (is.null(fiveAsgr$exon_rank))
    stop("fiveUTRs need column called exon_rank see ?makeTranscriptDbFromGFF")
  ##TODO: I should make this optional,
  #so that we dont need the exon_rank column..
  firstExons <- fiveAsgr[fiveAsgr$exon_rank == 1]
  posIDs <- strandBool(firstExons)
  promo <- promoters(firstExons, upstream = extension)
  start(firstExons[posIDs]) <- start(promo[posIDs])
  end(firstExons[!posIDs]) <- end(promo[!posIDs])
  return(firstExons)
}


#' Find max peak for each transcript,
#' returns as data.table, without names, but with index
#' @param cageOverlaps The cageOverlaps between cage and extended 5' leaders
#' @param filteredCage The filtered raw cage-data
#'  used to reassign 5' leaders
#' @importFrom data.table as.data.table
#' @return a data.table of max peaks
#'
findMaxPeaks <- function(cageOverlaps, filteredCage) {

  dt <- as.data.table(filteredCage)
  dt <- dt[from(cageOverlaps)]
  dt$to <- to(cageOverlaps)
  if (nrow(dt) == 0) return(dt)

  maxPeaks <- dt[, max(score), by = to]

  names(maxPeaks) <- c("to", "score")
  maxPeaks <-  merge(maxPeaks, dt)
  return(maxPeaks[!duplicated(maxPeaks$to)])
}


#' Finds max peaks per trancsript from reads in the cagefile
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param cageData The location of the cage-file
#' @param extension The number of basses upstream to add on transcripts
#' @return a Hits object
#'
findNewTSS <- function(fiveUTRs, cageData, extension) {

  shiftedfiveUTRs <- extendsTSSexons(fiveUTRs, extension)
  cageOverlaps <- findOverlaps(query = cageData, subject = shiftedfiveUTRs)
  maxPeakPosition <- findMaxPeaks(cageOverlaps, cageData)
  return(maxPeakPosition)
}


#' After all transcript start sites have been updated from cage,
#' put GRangesList back together
#' @param firstExons The first exon of every transcript from 5' leaders
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @return a GRangesList
#'
assignFirstExons <- function(firstExons, fiveUTRs){

  fiveAsgr <- unlist(fiveUTRs, use.names = TRUE)
  fiveAsgr[fiveAsgr$exon_rank == 1] <- firstExons
  return(relist(fiveAsgr, fiveUTRs))
}


#' add cage max peaks as new transcript start sites for each 5' leader
#' (*) strands are not supported, since direction must be known.
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param maxPeakPosition The max peak for each 5' leader found by cage
#' @return a GRanges object of first exons
#'
addNewTSSOnLeaders <- function(fiveUTRs, maxPeakPosition){

  fiveAsgr <- unlistGrl(fiveUTRs)
  firstExons <- fiveAsgr[fiveAsgr$exon_rank == 1]

  maxPeakPosition$names <- names(firstExons[maxPeakPosition$to])
  posIDs <- maxPeakPosition$to[maxPeakPosition$strand == "+"]
  minIDs <- maxPeakPosition$to[maxPeakPosition$strand == "-"]

  firstExons[posIDs] <- resize(
    firstExons[posIDs],
    width = end(firstExons[posIDs]) - maxPeakPosition$start[
      maxPeakPosition$strand == "+"] + 1,
    fix = "end")
  firstExons[minIDs] <- resize(
    firstExons[minIDs],
    width = maxPeakPosition$end[
      maxPeakPosition$strand == "-"] - start(firstExons[minIDs]) + 1,
    fix = "end")
  # Might need an chromosome boundary here? current test show no need.
  return(firstExons)
}


#' Reassign all Transcript Start Sites (TSS)
#'
#' Given a GRangesList of 5' UTRs or transcripts, reassign the start
#' postitions using max peaks from CageSeq data. A max peak is defined as new
#' TSS if it is within boundary of 5' leader range, specified by
#' \code{extension} in bp. A max peak must also be higher than minimum
#' CageSeq peak cutoff specified in \code{filterValue}. The new TSS will then
#' be the positioned where the cage read (with highest read count in the
#' interval).
#'
#' @param fiveUTRs (GRangesList) The 5' leaders or transcript sequences
#' @param cage Either a filePath for CageSeq file, or already loaded
#' CageSeq peak data as GRanges.
#' @param extension The maximum number of basses upstream of the TSS to search
#' for CageSeq peak.
#' @param filterValue The minimum number of reads on cage position,
#' for it to be counted as possible new tss.
#' (represented in score column in CageSeq data) to
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
#' reassignTSSbyCage(fiveUTRs, cage)
#'
reassignTSSbyCage <- function(fiveUTRs, cage, extension = 1000,
                              filterValue = 1, cds = NULL) {
  validGRL(class(fiveUTRs), "fiveUTRs")
  if (is.character(cage)) {
    filteredCage <- filterCage(fread.bed(cage),
                               filterValue) # get the cage data
  } else if (class(cage) == "GRanges") {
    filteredCage <- filterCage(cage, filterValue)
  } else {
    stop("Cage-file must be either a valid character",
         " filepath or GRanges object.")
  }
  # check that seqnames match
  filteredCage <- matchSeqlevels(filteredCage, fiveUTRs)
  maxPeakPosition <- findNewTSS(fiveUTRs, filteredCage, extension)
  fiveUTRs <- assignFirstExons(
    addNewTSSOnLeaders(fiveUTRs, maxPeakPosition), fiveUTRs)
  if(!is.null(cds)) fiveUTRs <- addFirstCdsOnLeaderEnds(fiveUTRs, cds)
  return(fiveUTRs)
}

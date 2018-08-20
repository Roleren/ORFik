
#' Filter peak of cage-data by value
#' @param rawCage The raw cage-data
#' @param filterValue The number of counts(score) to filter on
#'  for a tss to pass as hit
#' @param fiveUTRs a GRangesList (NULL), if added will filter out cage reads by
#' these following rules:
#' all reads in region (-5:-1, 1:5) for each tss will be removed, removes noise.
#' @return the filtered Granges object
#'
filterCage <- function(rawCage, filterValue = 1, fiveUTRs = NULL) {
  if (filterValue == 0) {
    return(rawCage)
  }
  if (is.null(rawCage$score)) stop("Found no 'score' column in the CageSeq.")
  filteredCage <- rawCage[rawCage$score > filterValue, ] #filter on score

  if (!is.null(fiveUTRs)) {
    # check that seqnames match
    seqlevelsStyle(filteredCage) <- seqlevelsStyle(fiveUTRs)
    tss <- startSites(fiveUTRs, asGR = TRUE, is.sorted = TRUE)
    # remove all reads in region (-5:-1, 1:5) of tss
    filteredCage <- filteredCage[!((countOverlaps(filteredCage, tss,
                                                  maxgap = 4) > 0) &
                      (countOverlaps(filteredCage, tss, maxgap = 0) == 0))]
  }
  return(filteredCage)
}

#' Extend first exon of each transcript with length specified
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param extension The number of basses to extend transcripts upstream
#' @return GRangesList object of fiveUTRs
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
  start(firstExons[posIDs]) <- pmax(rep(1, length(promo[posIDs])),
                                    start(promo[posIDs]))
  end(firstExons[!posIDs]) <- end(promo[!posIDs])
  fiveAsgr[fiveAsgr$exon_rank == 1] <- firstExons

  return(relist(fiveAsgr, fiveUTRs))
}

#' Restrict extension of 5' UTRs to closest upstream TSS
#'
#' That is minimum(CAGE extension, upstream TSS)
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @return GRangesList object of restricted fiveUTRs
restrictTSSByUpstreamLeader <- function(fiveUTRs){
  #TODO
  #find overlapping ranges
  # fiveUTRs
  # tss region

  overlaps <- findOverlaps(fiveUTRs, IRanges(startSites(fiveUTRs, asGR = T, keep.names = F, is.sorted = T), width = 5))
  overlaps <- overlaps[to(overlaps) != from(overlaps)]
  # for each leader, with more than 1 overlapping tss, choose most downstream


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
#' @param cageData The CAGE as GRanges object
#' @param extension The number of basses upstream to add on transcripts
#' @return a Hits object
#'
findNewTSS <- function(fiveUTRs, cageData, extension) {

  shiftedfiveUTRs <- extendsTSSexons(fiveUTRs, extension)
  #TODO add function that remove hits upstream of closest upstream tx
  #restrictTSSByUpstreamLeader()
  cageOverlaps <- findOverlaps(query = cageData, subject = shiftedfiveUTRs)
  maxPeakPosition <- findMaxPeaks(cageOverlaps, cageData)
  return(maxPeakPosition)
}

#' add cage max peaks as new transcript start sites for each 5' leader
#' (*) strands are not supported, since direction must be known.
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param maxPeakPosition The max peak for each 5' leader found by cage
#' @return a GRanges object of first exons
#'
addNewTSSOnLeaders <- function(fiveUTRs, maxPeakPosition){

  newTSS <- startSites(fiveUTRs, asGR = FALSE, keep.names = FALSE,
                       is.sorted = TRUE)
  newTSS[maxPeakPosition$to[maxPeakPosition$strand == "+"]] <-
    maxPeakPosition$start[maxPeakPosition$to[maxPeakPosition$strand == "+"]]
  newTSS[maxPeakPosition$to[maxPeakPosition$strand == "-"]] <-
    maxPeakPosition$end[maxPeakPosition$to[maxPeakPosition$strand == "-"]]
  fiveUTRsNew <- downstreamFromPerGroup(fiveUTRs, newTSS)

  return(fiveUTRsNew)
}


#' Reassign all Transcript Start Sites (TSS)
#'
#' Given a GRangesList of 5' UTRs or transcripts, reassign the start
#' sites using max peaks from CageSeq data. A max peak is defined as new
#' TSS if it is within boundary of 5' leader range, specified by
#' `extension` in bp. A max peak must also be higher than minimum
#' CageSeq peak cutoff specified in `filterValue`. The new TSS will then
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
#' (represented in score column in CageSeq data)
#' If you already filtered, set it to 0.
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
#'   seqnames = "1",
#'   ranges =  IRanges::IRanges(500, width = 1),
#'   strand = "+",
#'   score = 10) # <- Number of tags (reads) per position
#' # notice also that seqnames use different naming, this will be fixed by ORFik
#' # finally reassign TSS for fiveUTRs
#' reassignTSSbyCage(fiveUTRs, cage)
#'
reassignTSSbyCage <- function(fiveUTRs, cage, extension = 1000,
                              filterValue = 1) {
  validGRL(class(fiveUTRs), "fiveUTRs")
  if (is.character(cage)) {
    filteredCage <- filterCage(fread.bed(cage),
                               filterValue, fiveUTRs) # get the cage data
  } else if (is(cage, "GRanges")) {
    filteredCage <- filterCage(cage, filterValue, fiveUTRs)
  } else {
    stop("Cage-file must be either a valid character",
         " filepath or GRanges object.")
  }

  maxPeakPosition <- findNewTSS(fiveUTRs, filteredCage, extension)
  fiveUTRs <- addNewTSSOnLeaders(fiveUTRs, maxPeakPosition)
  return(fiveUTRs)
}

#' Input a txdb and reassign the TSS for each transcript by CAGE
#'
#'' Given a TxDb object, reassign the start site per transcript
#' using max peaks from CageSeq data. A max peak is defined as new
#' TSS if it is within boundary of 5' leader range, specified by
#' `extension` in bp. A max peak must also be higher than minimum
#' CageSeq peak cutoff specified in `filterValue`. The new TSS will then
#' be the positioned where the cage read (with highest read count in the
#' interval).
#' @param txdb a TxDb object, normally from a gtf file.
#' @param cage Either a filePath for CageSeq file, or already loaded
#' CageSeq peak data as GRanges.
#' @param extension The maximum number of basses upstream of the TSS to search
#' for CageSeq peak.
#' @param filterValue The minimum number of reads on cage position,
#' for it to be counted as possible new tss.
#' (represented in score column in CageSeq data)
#' If you already filtered, set it to 0.
#' @importFrom data.table setkeyv
#' @export
#' @examples
#'  \dontrun{
#'  library(GenomicFeatures)
#'  # Get the gtf txdb file
#'  txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'  package = "GenomicFeatures")
#'  txdb <- loadDb(txdbFile)
#'  cagePath <- system.file("extdata", "cage-seq-heart.bed.bgz",
#'  package = "ORFik")
#'  reassignTxDbByCage(txdb, cagePath)
#'  }
#' @return a TxDb obect of reassigned transcripts
reassignTxDbByCage <- function(txdb, cage, extension = 1000,
                               filterValue = 1) {
  if (!is(txdb,"TxDb")) stop("txdb must be a TxDb object")
  fiveUTRs <- fiveUTRsByTranscript(txdb, use.names = TRUE)
  fiveUTRs <- reassignTSSbyCage(fiveUTRs, cage, extension, filterValue)
  starts <- startSites(fiveUTRs, keep.names = TRUE)

  txList <- as.list(txdb)
  # find all transcripts with 5' UTRs
  hitsTx <- txList$transcripts$tx_name %in% names(fiveUTRs)
  idTx <- which(hitsTx)
  # reassign starts for positive strand
  txList$transcripts$tx_start[hitsTx][strandBool(fiveUTRs)] <-
    starts[strandBool(fiveUTRs)]
  txList$splicings$exon_start[(txList$splicings$tx_id %in% idTx) &
                                (txList$splicings$exon_rank == 1) &
                                (txList$splicings$exon_strand == "+")] <-
    starts[strandBool(fiveUTRs)]
  # reassign stops for negative strand
  txList$transcripts$tx_end[hitsTx][!strandBool(fiveUTRs)] <-
    starts[!strandBool(fiveUTRs)]
  txList$splicings$exon_end[(txList$splicings$tx_id %in% idTx) &
                              (txList$splicings$exon_rank == 1) &
                              (txList$splicings$exon_strand == "-")] <-
    starts[!strandBool(fiveUTRs)]

  # remake exon ids
  DT <- data.table(start = txList$splicings$exon_start,
                   end = txList$splicings$exon_end,
                   chr = txList$splicings$exon_chr,
                   strand = txList$splicings$exon_strand,
                   id = seq.int(1,length(txList$splicings$exon_start)))
  setkeyv(DT, c("start", "end", "chr", "strand"))
  d <- duplicated(DT[,.(start, end, chr, strand)])
  c <- rep.int(1,length(d))
  # find unique exon ids
  for(x in seq.int(2,length(d))){
    if (d[x]) {
      c[x] <- c[x-1]
    } else {
      c[x] <- c[x-1]+1
    }
  }

  # reassign exon ids
  txList$splicings$exon_id <- as.integer(c[order(DT$id)])
  # Since exons have changed, their official exon names can not be preserved.
  txList$splicings$exon_name <- NULL

  return(do.call(makeTxDb, txList))
}

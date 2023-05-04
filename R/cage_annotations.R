
#' Filter peak of cage-data by value
#' @inheritParams reassignTSSbyCage
#' @param fiveUTRs a GRangesList (NULL), if added will filter out cage reads by
#' these following rules:
#' all reads in region (-5:-1, 1:5) for each tss will be removed, removes noise.
#' @return the filtered Granges object
#' @keywords internal
#'
filterCage <- function(cage, filterValue = 1, fiveUTRs = NULL,
                       preCleanup = TRUE) {
  cage <- fimport(cage, seqlevelsStyle(fiveUTRs))
  noScore <- is.null(score(cage))
  if (noScore | (max(width(cage)) > 1)) {
    cage <- convertToOneBasedRanges(cage, addScoreColumn = noScore)
  }

  if (filterValue == 0) {
    return(cage)
  }
  if (is.null(cage$score)) stop("Found no 'score' column in the 'cage' input!")
  filteredCage <- cage[cage$score > filterValue, ] #filter on score

  if (!is.null(fiveUTRs) & preCleanup) {
    # check that seqnames match
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
#' @keywords internal
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

#' Restrict extension of 5' UTRs to closest upstream leader end
#'
#' Basicly this function restricts all startSites, to the upstream GRangesList
#' objects end. Usually leaders, for CAGE.
#' Example: leader1: start on 10, leader2: stop on 8,
#' extend leader1 to 5 -> this function will resize leader1
#' to 9, to be outside leader2, so that CAGE reads can not wrongly overlap.
#' @param fiveUTRs The 5' leader sequences as GRangesList
#' @param shiftedfiveUTRs The 5' leader sequences as GRangesList
#'  shifted by CAGE
#' @return GRangesList object of restricted fiveUTRs
#' @keywords internal
restrictTSSByUpstreamLeader <- function(fiveUTRs, shiftedfiveUTRs) {
  if (length(fiveUTRs) != length(shiftedfiveUTRs))
    stop("fiveUTRs and shifted must be equal size!")

  startSitesGR <- startSites(fiveUTRs, asGR = TRUE, is.sorted = TRUE)
  stopSitesGR <- stopSites(fiveUTRs, asGR = TRUE, is.sorted = TRUE)

  overlaps <- findOverlaps(fiveUTRs, stopSitesGR)
  overlapsShifted <- findOverlaps(shiftedfiveUTRs, stopSitesGR)
  overlaps <- overlaps[from(overlaps) != to(overlaps)]
  overlapsShifted <- overlapsShifted[from(overlapsShifted) !=
                                       to(overlapsShifted)]

  # remove overlaps also in overlapsShifted
  if (length(overlapsShifted) > 0 & length(overlaps) > 0)
    overlapsShifted <- overlapsShifted[!(overlapsShifted %in% overlaps)]

  # for all duplicates, choose closest
  dt <- data.table(from = from(overlapsShifted), to = to(overlapsShifted))
  dt$distance <- distance(startSitesGR[from(overlapsShifted)],
                          stopSitesGR[to(overlapsShifted)])
  minDT <- dt[, .I[which.min(distance)], by=from]

  # assign closest
  shiftedfiveUTRs[minDT$from]  <- downstreamFromPerGroup(
    shiftedfiveUTRs[minDT$from],
    start(stopSitesGR[to(overlapsShifted)[minDT$V1]]))

  return(shiftedfiveUTRs)
}

#' Find max peak for each transcript,
#' returns as data.table, without names, but with index
#' @param cageOverlaps The cageOverlaps between cage and extended 5' leaders
#' @param filteredCage The filtered raw cage-data
#'  used to reassign 5' leaders
#' @importFrom data.table as.data.table
#' @return a data.table of max peaks
#' @keywords internal
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
#' @param extension The number of basses to extends transcripts upstream.
#' @param restrictUpstreamToTx a logical (FALSE), if you want to restrict
#'  leaders to not extend closer than 5 bases from closest upstream leader,
#'  set this to TRUE.
#' @return a Hits object
#' @keywords internal
findNewTSS <- function(fiveUTRs, cageData, extension, restrictUpstreamToTx) {

  shiftedfiveUTRs <- extendsTSSexons(fiveUTRs, extension)
  if (restrictUpstreamToTx){
    shiftedfiveUTRs <- restrictTSSByUpstreamLeader(fiveUTRs, shiftedfiveUTRs)
  }

  cageOverlaps <- findOverlaps(cageData, shiftedfiveUTRs, type = "within")
  maxPeakPosition <- findMaxPeaks(cageOverlaps, cageData)
  return(maxPeakPosition)
}

#' Add cage max peaks as new transcript start sites for each 5' leader
#' (*) strands are not supported, since direction must be known.
#' @param maxPeakPosition The max peak for each 5' leader found by cage
#' @inheritParams reassignTSSbyCage
#' @return a GRanges object of first exons
#' @keywords internal
addNewTSSOnLeaders <- function(fiveUTRs, maxPeakPosition, removeUnused,
                               cageMcol) {

  newTSS <- startSites(fiveUTRs, asGR = FALSE, keep.names = FALSE,
                       is.sorted = TRUE)
  newTSS[maxPeakPosition$to[maxPeakPosition$strand == "+"]] <-
    maxPeakPosition$start[maxPeakPosition$strand == "+"]
  newTSS[maxPeakPosition$to[maxPeakPosition$strand == "-"]] <-
    maxPeakPosition$end[maxPeakPosition$strand == "-"]

  if (cageMcol) {
    gr <- unlist(fiveUTRs, use.names = FALSE)
    cageM <- rep.int(0, length(gr))
    groupOld <- groupings(fiveUTRs)
    groupHits <- groupings(fiveUTRs[maxPeakPosition$to])

    cageM[groupOld %in% maxPeakPosition$to] <- maxPeakPosition$score[groupHits]
    mcols(gr) <- DataFrame(row.names = names(gr), mcols(gr),
                           cage = cageM)
    fiveUTRs <- groupGRangesBy(gr, names(fiveUTRs)[groupOld])
  }
  if (removeUnused) {
    fiveUTRsNew <- downstreamFromPerGroup(fiveUTRs[maxPeakPosition$to],
                                          newTSS[maxPeakPosition$to])
  } else {
    fiveUTRsNew <- downstreamFromPerGroup(fiveUTRs, newTSS)
  }
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
#' interval). If removeUnused is TRUE, leaders without cage hits, will be
#' removed, if FALSE the original TSS will be used.
#'
#' Note: If you used CAGEr, you will get reads of a probability region, with
#' always score of 1. Remember then to set filterValue to 0. And you should use
#' the 5' end of the read as input, use: ORFik:::convertToOneBasedRanges(cage)
#' NOTE on filtervalue: To get high quality TSS, set filtervalue to median
#' count of reads overlapping per leader. This will make you discard a lot of
#' new TSS positions though. I usually use 10 as a good standard.
#'
#' TIP: do summary(countOverlaps(fiveUTRs, cage)) so you can find a good
#' cutoff value for noise.
#' @param fiveUTRs (GRangesList) The 5' leaders or full transcript sequences
#' @param cage Either a filePath for the CageSeq file as .bed .bam or .wig,
#' with possible compressions (".gzip", ".gz", ".bgz"), or already loaded
#' CageSeq peak data as GRanges or GAlignment.
#' NOTE: If it is a .bam file, it will add a score column by running:
#' convertToOneBasedRanges(cage, method = "5prime", addScoreColumn = TRUE)
#' The score column is then number of replicates of read, if score column is
#' something else, like read length, set the score column to NULL first.
#' @param extension The maximum number of basses upstream of the TSS to search
#' for CageSeq peak.
#' @param filterValue The minimum number of reads on cage position, for it to
#' be counted as possible new tss. (represented in score column in
#' CageSeq data) If you already filtered, set it to 0.
#' @param restrictUpstreamToTx a logical (FALSE). If TRUE: restrict leaders to
#' not extend closer than 5 bases from closest upstream leader, set this
#' to TRUE.
#' @param removeUnused logical (FALSE), if False:  (standard is to set them to
#' original annotation), If TRUE: remove leaders that did not have any cage
#' support.
#' @param preCleanup logical (TRUE), if TRUE,
#' remove all reads in region (-5:-1, 1:5) of all original tss in leaders.
#' This is to keep original TSS if it is only +/- 5 bases from the original.
#' @param cageMcol a logical (FALSE), if TRUE, add a meta column to the
#' returned object with the raw CAGE counts in support for new TSS.
#' @return a GRangesList of newly assigned TSS for fiveUTRs,
#'  using CageSeq data.
#' @family CAGE
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
#' # make fake CAGE data from promoter of 5' leaders, notice score column
#' cage <- GenomicRanges::GRanges(
#'   seqnames = "1",
#'   ranges =  IRanges::IRanges(500, width = 1),
#'   strand = "+",
#'   score = 10) # <- Number of tags (reads) per position
#' # notice also that seqnames use different naming, this is fixed by ORFik
#' # finally reassign TSS for fiveUTRs
#' reassignTSSbyCage(fiveUTRs, cage)
#' # See vignette for example using gtf file and real CAGE data.
#'
reassignTSSbyCage <- function(fiveUTRs, cage, extension = 1000,
                              filterValue = 1, restrictUpstreamToTx = FALSE,
                              removeUnused = FALSE, preCleanup = TRUE,
                              cageMcol = FALSE) {
  validGRL(class(fiveUTRs), "fiveUTRs")
  filteredCage <- filterCage(cage, filterValue, fiveUTRs, preCleanup)

  maxPeakPosition <- findNewTSS(fiveUTRs, filteredCage, extension,
                                restrictUpstreamToTx)

  return(addNewTSSOnLeaders(fiveUTRs, maxPeakPosition, removeUnused, cageMcol))
}

#' Input a txdb and reassign the TSS for each transcript by CAGE
#'
#' Given a TxDb object, reassign the start site per transcript
#' using max peaks from CageSeq data. A max peak is defined as new
#' TSS if it is within boundary of 5' leader range, specified by
#' `extension` in bp. A max peak must also be higher than minimum
#' CageSeq peak cutoff specified in `filterValue`. The new TSS will then
#' be the positioned where the cage read (with highest read count in the
#' interval).
#'
#' Note: If you used CAGEr, you will get reads of a probability region, with
#' always score of 1. Remember then to set filterValue to 0. And you should use
#' the 5' end of the read as input, use: ORFik:::convertToOneBasedRanges(cage)
#' @inheritParams loadTxdb
#' @inheritParams reassignTSSbyCage
#' @importFrom data.table setkeyv
#' @family CAGE
#' @export
#' @examples
#'  \dontrun{
#'  library(GenomicFeatures)
#'  # Get the gtf txdb file
#'  txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'  package = "GenomicFeatures")
#'  cagePath <- system.file("extdata", "cage-seq-heart.bed.bgz",
#'  package = "ORFik")
#'  reassignTxDbByCage(txdbFile, cagePath)
#'  }
#' @return a TxDb obect of reassigned transcripts
reassignTxDbByCage <- function(txdb, cage, extension = 1000,
                               filterValue = 1, restrictUpstreamToTx = FALSE,
                               removeUnused = FALSE, preCleanup = TRUE) {
  txdb <- loadTxdb(txdb)
  fiveUTRs <- fiveUTRsByTranscript(txdb, use.names = TRUE)
  fiveUTRs <- reassignTSSbyCage(fiveUTRs, cage, extension, filterValue,
                                restrictUpstreamToTx, removeUnused, preCleanup)

  txList <- as.list(txdb)
  # find all transcripts with 5' UTRs
  txList <- updateTxdbStartSites(txList, fiveUTRs, removeUnused)

  # reassign exon ids
  txList$splicings$exon_id <- remakeTxdbExonIds(txList)
  # Since exons have changed, their official exon names can not be preserved.
  txList$splicings$exon_name <- NULL
  # Copy useful metadata
  original_meta <- metadata(txdb)
  original_meta <- original_meta[!(original_meta$name %in% c("Db type", "Supporting package","DBSCHEMAVERSION",
                                                             "Creation time",
                                                             "RSQLite version at creation time",
                                                             "GenomicFeatures version at creation time")),]
  txList$metadata <- original_meta

  return(do.call(makeTxDb, txList))
}

#' Input a txdb and add a 5' leader for each transcript, that does not have one.
#'
#' For all cds in txdb, that does not have a 5' leader:
#' Start at 1 base upstream of cds and use CAGE, to assign leader start.
#' All these leaders will be 1 exon based, if you really want exon splicings,
#' you can use exon prediction tools, or run sequencing experiments.
#'
#' Given a TxDb object, reassign the start site per transcript
#' using max peaks from CageSeq data. A max peak is defined as new
#' TSS if it is within boundary of 5' leader range, specified by
#' `extension` in bp. A max peak must also be higher than minimum
#' CageSeq peak cutoff specified in `filterValue`. The new TSS will then
#' be the positioned where the cage read (with highest read count in the
#' interval). If no CAGE supports a leader, the width will be set to 1 base.
#'
#' @inheritParams reassignTxDbByCage
#' @param pseudoLength a numeric, default 1. Either if no CAGE supports
#' the leader, or if CAGE is set to NULL,
#' add a pseudo length for all the UTRs. Will not extend a leader if
#' it would make it go outside the defined seqlengths of the genome.
#' So this length is not guaranteed for all!
#' @importFrom data.table setkeyv
#' @family CAGE
#' @return a TxDb obect of reassigned transcripts
#' @export
#' @examples
#' txdbFile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
#'  package = "GenomicFeatures")
#' cagePath <- system.file("extdata", "cage-seq-heart.bed.bgz",
#'  package = "ORFik")
#'
#' \dontrun{
#'   assignTSSByCage(txdbFile, cagePath)
#'   #Minimum 20 cage tags for new TSS
#'   assignTSSByCage(txdbFile, cagePath, filterValue = 20)
#'   # Create pseudo leaders for the ones without hits
#'   assignTSSByCage(txdbFile, cagePath, pseudoLength = 100)
#'   # Create only pseudo leaders (in example 2 leaders are added)
#'   assignTSSByCage(txdbFile, cage = NULL, pseudoLength = 100)
#' }
assignTSSByCage <- function(txdb, cage, extension = 1000,
                            filterValue = 1, restrictUpstreamToTx = FALSE,
                            removeUnused = FALSE, preCleanup = TRUE,
                            pseudoLength = 1) {
  if (!is.numeric(pseudoLength) & (length(pseudoLength) == 1)) {
    stop("pseudoLength must be length 1 and numeric!")
  }
  txdb <- loadTxdb(txdb)
  fiveUTRs <- loadRegion(txdb, part = "leaders")
  cds <- loadRegion(txdb, part = "cds")

  cds01 <- cds[!(names(cds) %in% names(fiveUTRs))]
  if (length(cds01) == 0) {
    message("- All mRNA have 5' leader, nothing to be done")
    return(txdb)
  }

  undefinedLeaders <- trim(promoters(startSites(cds01,asGR = T,
                                                is.sorted = TRUE,
                                                keep.names = T), 1, 0))
  cds_had_names <- !is.null(names(cds01))
  undefinedLeaders <- reduce(undefinedLeaders, TRUE,
                             with.revmap= cds_had_names)
  if (length(undefinedLeaders) == 0) {
    message("- All mRNA without leaders, had cds starting on end of contig,",
            " nothing to be done")
    return(txdb)
  } else if (length(undefinedLeaders) != length(cds01)) {
    message("Some mRNA had cds starting on end of contig, can not be added")
    message("Number of original CDSs: ", length(cds01))
    message("Number of valid 5' UTR genes: ", length(undefinedLeaders))
  }
  if (cds_had_names) {
    names(undefinedLeaders) <- names(cds01)[unlist(undefinedLeaders$revmap)]
    undefinedLeaders$revmap <- NULL
  }
  #make exon_rank col as integer
  undefinedLeaders$exon_rank <- rep.int(1L, length(undefinedLeaders))
  #make GRangesListfrom GRanges
  undefinedLeaders <- groupGRangesBy(undefinedLeaders)
  if (pseudoLength != 1) {
    undefinedLeaders <- extendLeaders(undefinedLeaders, extension = pseudoLength)
  }
  newLeaders <- if (!is.null(cage)) {
    reassignTSSbyCage(undefinedLeaders, cage, extension,
                      filterValue, restrictUpstreamToTx,
                      removeUnused, preCleanup)
  } else undefinedLeaders

  txList <- as.list(txdb)
  # find all transcripts with 5' UTRs
  txList <- updateTxdbStartSites(txList, newLeaders, removeUnused)

  # reassign exon ids
  txList$splicings$exon_id <- remakeTxdbExonIds(txList)
  # Since exons have changed, their official exon names can not be preserved.
  txList$splicings$exon_name <- NULL
  # Copy useful metadata
  original_meta <- metadata(txdb)
  original_meta <- original_meta[!(original_meta$name %in% c("Db type", "Supporting package","DBSCHEMAVERSION",
                                                             "Creation time",
                                                             "RSQLite version at creation time",
                                                             "GenomicFeatures version at creation time")),]
  txList$metadata <- original_meta
  return(do.call(makeTxDb, txList))
}

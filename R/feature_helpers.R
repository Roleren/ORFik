#' Get length of leaders ordered after oldTxNames
#'
#' Normally only a helper function for ORFik
#' @param fiveUTRs a GRangesList object of leaders
#' @param oldTxNames a character vector of names to group fiveUTRs by.
#' @return a GRangesList of reordered leaders.
#'
findCageUTRFivelen <- function(fiveUTRs, oldTxNames){
  newfiveprimeLen <- widthPerGroup(fiveUTRs)
  return(newfiveprimeLen[match(oldTxNames, names(newfiveprimeLen))])
}


#' Get transcript lengths
#'
#' A helper function for easy length retrieval
#' @param Gtf a TxDb object of a gtf file
#' @param changedFiveUTRs a GRangesList object of leaders.
#' Only add this if you used cage data or other things to change the
#' leaders, therefor we need it to update transcript lengths.
#' @return a vector of transcript lengths
#'
txLen <- function(Gtf = NULL, changedFiveUTRs = NULL){
  tx_len_temp <- transcriptLengths(Gtf)[, c("tx_name", "tx_len")]
  tx_len <- tx_len_temp[, "tx_len"]

  if (!is.null(changedFiveUTRs)) {
    if (!is.null(Gtf)) {
      new5Length <- findCageUTRFivelen(changedFiveUTRs, tx_len_temp$tx_name)
      tx_len_temp <- transcriptLengths(Gtf, TRUE, TRUE, TRUE)
      tx_len_temp$utr5_len <- new5Length
      tx_len <- tx_len_temp$utr5_len +
      tx_len_temp$cds_len + tx_len_temp$utr3_len
    }
  }
  names(tx_len) <- tx_len_temp$tx_name
  return(tx_len)
}


#' Get hits per codon
#'
#' Helper for entropy function, normally not used directly
#' Seperate each group into tuples (abstract codons)
#' Gives sum for each tuple within each group
#' Example: c(1,0,0,1), with reg_len = 2, gives
#' c(1,0) and c(0,1), these are summed and returned as list
#' @param countList a Rle of count repetitions (000,1,00,1 etc)
#' @param reg_len integer vector, size of runs
#' @param runLengths integer vector, duplications per run
#' @return a list of codon sums
#'
codonSumsPerGroup <- function(countList, reg_len,
                              runLengths ){

  len <- BiocGenerics::lengths(countList)
  if (length(len) > 1) { # if more than 1 hit total
    acums <- cumsum(as.numeric(len[seq.int(1, length(len)-1)]))
    acums <- rep.int(c(1,acums), runLengths)
  } else { # special case for 1 group only
    acums <- 1
  }

  # Need to reassign variables to equal length,
  # to be able to vectorize
  # h: The sequences we make the tuplets per orf,
  # if h[1] is: c(0,1,2,3,4,5) and reg_len[1] is: c(3,3)
  # you get int_seqs: ->  1: c(1,2,3 , 4,5,6) <- 2 triplets
  reg_len <- rep.int(reg_len, runLengths)
  h <- unlist(lapply(runLengths - 1, function(x) {
    seq.int(0, x)
  }), use.names = FALSE)

  which_reads_start <- (acums + h * reg_len)
  which_reads_end <- (h * reg_len + (reg_len + acums - 1))
  # the actual triplets ->
  int_seqs <- lapply(seq_along(which_reads_start), function(x) {
    which_reads_start[x]:which_reads_end[x]
  })

  unlintcount <- unlist(IntegerList(countList), use.names = FALSE)
  # get the assigned tuplets per orf, usually triplets
  triplets <- lapply(int_seqs, function(x) {
    unlintcount[x]
  })
  return(sum(IntegerList(triplets)))
}


#' Create normalizations of counts
#'
#' A helper for [fpkm()]
#' Normally use function [fpkm()], if you want unusual normalization
#' , you can use this.
#' Short for: Fragments per kilobase of transcript per million fragments
#' Normally used in Translations efficiency calculations
#' @references doi: 10.1038/nbt.1621
#' @param counts a list, # of read hits per group
#' @param lengthSize a list of lengths per group
#' @param librarySize a numeric of size 1, the # of reads in library
#' @family features
#' @return a numeric vector
#'
fpkm_calc <- function(counts, lengthSize, librarySize){
  return((as.numeric(counts) * (10^9)) /
           (as.numeric(lengthSize) * as.numeric(librarySize)))
}

#' Get -20,20 start region as DNA characters per gr group
#'
#' @param grl a GRangesList to find regions
#' @param tx a GRangesList of transcripts containing grl
#' @param faFile a FaFile from the fasta file, see ?FaFile.
#'  Can also be path to fastaFile with fai file in same dir.
#' @param groupBy (NULL) column to group grl by, if NULL group by names(gr)
#' @return a character vector of start regions
startRegionString <- function(grl, tx, faFile, groupBy = NULL){
  gr <- startSites(grl, TRUE, TRUE, TRUE)

  grl <- groupGRangesBy(windowPerGroup(grl, tx, 20, 20), groupBy)


  return(as.character(txSeqsFromFa(grl, faFile, is.sorted = TRUE)))
}

#' Hits from reads
#'
#' Finding GRanges groups that have overlap hits with reads
#' @param grl a GRanges or GRangesList
#' @param reads a GAlignment or GRanges object with reads
#' @param keep.names logical (F), keep names or not
#' @return a list of logicals, T == hit, F == no hit
#'
hasHits <- function(grl, reads, keep.names = FALSE) {
  overlaps <- countOverlaps(grl, reads)
  if (!keep.names) {
    names(overlaps) <- NULL
  }
  return(overlaps > 0)
}

#' Helper Function to check valid RNA input
#' @param class, the given class of RNA object
#' @return NULL, stop if unvalid object
#'
checkRNA <- function(class){
  if (is.null(class) || (class == "NULL")) {
    message("No RNA added, skipping feature te and fpkm of RNA, ",
            "also ribosomeReleaseScore will also be not normalized best ",
            "way possible.")
  } else {
    if (class != "GAlignments" & class != "GRanges") {
      stop("RNA must be either GAlignments or GRanges")
    }
  }
}


#' Helper Function to check valid RFP input
#' @param class, the given class of RFP object
#' @return NULL, stop if invalid object
#'
checkRFP <- function(class) {
  if (class != "GAlignments" & class != "GRanges") {
    stop("RFP must be either GAlignments or GRanges")
  }
}


#' Helper function to check valid combinations of extension and cageFiveUTRs
#' @param extension a numeric/integer to reassign 5' utrs.
#' @param cageFiveUTRs a GRangesList, if you used cage-data to extend 5' utrs,
#' @return NULL, stop if invalid object
validExtension <- function(extension, cageFiveUTRs) {
  if (is.null(extension)) {
    stop("please specify extension, to avoid bugs ",
         "if you did not use cage, set it to 0, ",
         "standard cage extension is 1000")
  } else if (!is.numeric(extension) && !is.integer(extension)) {
      stop("extension must be numeric or integer")
  }
  if (extension != 0 && !is.grl(cageFiveUTRs)) {
    stop("if extension is not 0, then cageFiveUTRs must be defined")
  }
}

#' Find proportion of reads per position in ORF TIS window
#'
#' Proportion defined as:
#' average count per position in -20,20 normalized by counts per gene
#'
#' This pattern can be averaged on CDS's, to find other ORFs.
#' When detecting new ORFs, this CDS average can be used as a template.
#' @param grl a \code{\link{GRangesList}} object
#'  with usually new ORFs, but can also be
#'  either leaders, cds', 3'UTRs.
#' @param tx a GrangesList of transcripts covering grl.
#' @param footprints ribo seq reads as GAlignment, GRanges
#'  or GRangesList object.
#' @param onlyProportion a logical (FALSE), return whole data.frame or only
#'  proportions
#' @param average a logical (FALSE), take average over coverage in all grl ?
#' @param pShifted a logical (TRUE), are riboseq reads p-shifted?
#' @param keep.names a logical(FALSE), only applies when onlyProportion
#'  is TRUE.
#' @return a data.frame with lengths by coverage / vector of proportions
#' @importFrom BiocGenerics Reduce
#'
riboTISCoverageProportion <- function(grl, tx, footprints,
                                      onlyProportion = FALSE, average = FALSE,
                                        pShifted = TRUE, keep.names = FALSE){
  upStart <- if (pShifted) 5 else 20
  downStop <- if (pShifted) 20 else 5

  windowSize <- upStart + downStop + 1
  window <- windowPerGroup(startSites(grl, TRUE, FALSE, TRUE), tx, upStart,
                           downStop)
  noHits <- widthPerGroup(window) < windowSize
  if (all(noHits)) {
    warning("no grl had valid window size!")
    return(RleList(Rle(values = 0, lengths = windowSize)))
  }
  window <- window[!noHits]
  rwidth <- readWidths(footprints)
  footprints <- footprints[rwidth < 31 & rwidth > 26]
  rwidth <- readWidths(footprints)
  allLengths <- sort(unique(rwidth))
  gr <- resize(granges(footprints), 1)
  lengthProportions <- c()

  unlTile <- unlistGrl(tile1(window, matchNaming = FALSE))
  if (!is.null(unlTile$names)) { # for orf case
    names(unlTile) <- unlTile$names
  }

  for (l in allLengths) {
    ends_uniq <- gr[rwidth == l]

    cvg <- overlapsToCoverage(unlTile, ends_uniq, FALSE, type = "within")

    cvg <- cvg /sum(cvg)
    cvg[is.nan(unlist(sum(runValue(cvg)), use.names = FALSE))] <-
      RleList(Rle(values = 0, lengths = windowSize))
    if (average) {
      cvg <- Reduce(`+`, cvg)
      lengthProportions <- c(lengthProportions, as.vector(cvg/sum(cvg)))
    } else {
      lengthProportions <- c(lengthProportions, cvg)
    }
  }
  if (!onlyProportion) {
    if (average) {
      lengths <- unlist(lapply(allLengths, function(x){rep.int(x,windowSize)}),
                        use.names = FALSE)
      df <- data.frame(prop = lengthProportions, length = lengths,
                       pos = rep(seq.int(-upStart,downStop),
                                 length(allLengths)))
      return(df)
    } else {
      warning("Can only return proportion when average == FALSE")
      return(lengthProportions)
    }
  } else {
    if (keep.names & !average) {
      names(lengthProportions[[1]]) <- names(window)
    }
    return(lengthProportions)
  }
}


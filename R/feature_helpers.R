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
#' @param h indices per tuple
#' @param indeces whole list of indices
#' @param L Lengths
#' @param N hit sums
#' @param reg_len size of runs
#' @param runLengths integers per run
#' @param countList a Rle of count repetitions (000,1,00,1 etc)
#' @return a list of codon sums
#'
codonSumsPerGroup <- function(h, indeces, L, N, reg_len,
                              runLengths, countList){
  reg_len <- lapply(indeces, function(x) {
    rep(reg_len[x], runLengths[x])
  })

  unlh <- unlist(h, use.names = FALSE)
  unlreg_len <- unlist(reg_len, use.names = FALSE)

  # make sequences for reads, start -> stop
  # gives a triplet reading, 1:3, 3,6

  if (length(L) > 1) { # if more than 1 hit total
    acums <- cumsum(as.numeric(L))
    unlacums <- unlist(lapply(seq(length(runLengths)-1), function(x) {
      rep(acums[x], runLengths[x+1])
    }), use.names = FALSE)
    unlacums <- unlist(c(rep(1, runLengths[1]), unlacums))
  } else { # special case for 1 group only
    unlacums <- 1
  }

  which_reads_start <- (unlacums + unlh * unlreg_len)
  which_reads_end <- (unlh * unlreg_len + (unlreg_len + unlacums -1))
  # the actual triplets ->
  int_seqs <- lapply(seq_along(which_reads_start), function(x) {
    which_reads_start[x]:which_reads_end[x]
  })

  intcountList <- IntegerList(countList)
  unlintcount <- unlist(unlist(intcountList, use.names = FALSE),
                        use.names = FALSE)
  # get the assigned tuplets per orf, usually triplets
  triplets <- lapply(int_seqs, function(x) {
    unlintcount[x]
  })
  tripletSums <- unlist(lapply(triplets, function(x) {
    sum(x)
  }), use.names = FALSE)

  return(tripletSums)
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

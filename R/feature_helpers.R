#' Get hits per codon
#'
#' Helper for entropy function, normally not used directly
#' Seperate each group into tuples (abstract codons)
#' Gives sum for each tuple within each group
#' Example: c(1,0,0,1), with reg_len = 2, gives
#' c(1,0) and c(0,1), these are summed and returned as list
#' @param countList a Rle of count repetitions (000, 1, 00, 1 etc)
#' @param reg_len integer vector, size of runs
#' @param runLengths integer vector, duplications per run
#' @return a list of codon sums
#'
codonSumsPerGroup <- function(countList, reg_len,
                              runLengths ) {
  # TODO: USE NEW SCORING TO MAKE THIS VECTOR SIMPLER
  len <- lengths(countList)
  if (length(len) > 1) { # if more than 1 hit total
    acums <- cumsum(as.numeric(len[seq.int(length(len) - 1)]))
    acums <- rep.int(c(1, acums), runLengths)
  } else { # special case for 1 group only
    acums <- 1
  }

  # Need to reassign variables to equal length, to be able to vectorize
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
fpkm_calc <- function(counts, lengthSize, librarySize) {
  return((as.numeric(counts) * (10^9)) /
           (as.numeric(lengthSize) * as.numeric(librarySize)))
}

#' Get start region as DNA-strings per GRanges group
#'
#' One window per start site, if upstream and downstream are both 0, then
#' only the startsite is returned.
#' @param grl a \code{\link{GRangesList}} of ranges to find regions in.
#' @inheritParams windowPerGroup
#' @param faFile a FaFile from the fasta file, see ?FaFile.
#'  Can also be path to fastaFile with fai file in same dir.
#' @return a character vector of start regions
startRegionString <- function(grl, tx, faFile, upstream = 20,
                              downstream = 20) {
  grl <- startRegion(grl, tx, is.sorted = TRUE, upstream, downstream)
  return(as.character(txSeqsFromFa(grl, faFile, is.sorted = TRUE)))
}

#' Hits from reads
#'
#' Finding GRanges groups that have overlap hits with reads
#' Similar to %over%.
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

#' Subset GRanges to get coverage.
#'
#' GRanges object should be beforehand
#' tiled to size of 1. This subsetting takes account for strand.
#' @param cov A coverage object from coverage()
#' @param y GRanges object for which coverage should be extracted
#' @return numeric vector of coverage of input GRanges object
#' @family features
#'
subsetCoverage <- function(cov, y) {
  cov1 <- cov[[as.vector(seqnames(y)[1])]]
  return(as.vector(cov1[ranges(y)]))
}

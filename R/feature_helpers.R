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

  len <- lengths(countList)
  if (length(len) > 1) { # if more than 1 hit total
    acums <- cumsum(as.numeric(len[seq.int(1, length(len)-1)]))
    acums <- rep.int(c(1,acums), runLengths)
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

#' Get -20,20 start region as DNA characters per GRanges group
#'
#' @param grl a GRangesList to find regions
#' @param tx a GRangesList of transcripts containing grl
#' @param faFile a FaFile from the fasta file, see ?FaFile.
#'  Can also be path to fastaFile with fai file in same dir.
#' @param groupBy (NULL) column to group grl by, if NULL group by names(gr)
#' @return a character vector of start regions
startRegionString <- function(grl, tx, faFile, groupBy = NULL) {
  gr <- startSites(grl, TRUE, TRUE, TRUE)
  grl <- groupGRangesBy(windowPerGroup(grl, tx, 20, 20), groupBy)

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
#' @param upStart upstream region boundary (5 or 20 as standard), relative
#'  (5, mean 5 upstream from TIS)
#' @param downStop downstream region boundary (5 or 20 as standard), relative
#'  (5, mean 5 downstream from TIS)
#' @return a data.frame with lengths by coverage / vector of proportions
#' @importFrom BiocGenerics Reduce
#'
riboTISCoverageProportion <- function(grl, tx, footprints,
                                      onlyProportion = FALSE, average = FALSE,
                                      pShifted = TRUE, keep.names = FALSE,
                                      upStart = if (pShifted) 5 else 20,
                                      downStop = if (pShifted) 20 else 5) {
  windowSize <- upStart + downStop + 1
  window <- windowPerGroup(startSites(grl, TRUE, TRUE, TRUE), tx, downStop,
                           upStart)
  noHits <- widthPerGroup(window) < windowSize
  if (all(noHits)) {
    warning("no grl had valid window size!")
    return(RleList(Rle(values = 0, lengths = windowSize)))
  }
  window <- window[!noHits]
  # fix names, find a better way to store this, should be a function.
  if (is.ORF(grl)) {
    names(window) <- names(grl[!noHits])
    g <- unlist(window, use.names = TRUE)
    names(g) <- sub("\\..*", "", names(g), perl = TRUE)
    mcols(g) <- DataFrame(row.names = names(g), names = names(g))
    window <- groupGRangesBy(g)
  }

  unlTile <- tile1(window, matchNaming = FALSE)
  if(length(unlTile) != length(window)) stop("Bad naming, most likely _
                                             is not used for ORF correctly")
  unlTile <- unlistGrl(unlTile)

  rwidth <- readWidths(footprints)
  footprints <- footprints[rwidth < 31 & rwidth > 26]
  rwidth <- readWidths(footprints)
  allLengths <- sort(unique(rwidth))
  if (!(is.gr_or_grl(footprints) & unique(width(footprints)) == 1)) {
    footprints <- resize(granges(footprints), 1)
  }

  lengthProportions <- c()
  for (l in allLengths) {
    ends_uniq <- footprints[rwidth == l]

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
                       pos = rep(seq.int(-upStart, downStop),
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

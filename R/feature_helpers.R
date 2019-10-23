#' Get read hits per codon
#'
#' Helper for entropy function, normally not used directly
#' Seperate each group into tuples (abstract codons)
#' Gives sum for each tuple within each group
#'
#' Example: counts c(1,0,0,1), with reg_len = 2, gives
#' c(1,0) and c(0,1), these are summed and returned as data.table
#' 10 bases, will give 3 codons, 1 base codons does not exist.
#' @param grl a GRangesList
#' @param reads a GRanges or GAlignment
#' @return a data.table with codon sums
#'
codonSumsPerGroup <- function(grl, reads) {
  dt <- scaledWindowPositions(grl, reads, numCodons(grl))
  dt[, `:=` (codonSums = score / sum(score)), by = genes]
  return(dt)
}


#' Create normalizations of read counts
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
#' @export
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

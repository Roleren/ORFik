#' CountOverlaps with weights
#'
#' Similar to countOverlaps, but takes an optional weight column.
#' This is usually the score column
#'
#' @param query IRanges, IRangesList, GRanges, GRangesList object.
#' Usually transcript a transcript region.
#' @param subject GRanges, GRangesList, GAlignment, usually reads.
#' @param weight (default: NULL), if defined either numeric or character name
#' of valid meta column in subject. If weight is single numeric, it is used
#' for all. A normall weight is the score column given as weight = "score".
#' GRanges("chr1", 1, "+", score = 5), would mean score column tells
#' that this alignment region was found 5 times.
#' @param ... additional arguments passed to countOverlaps/findOverlaps
#' @return a named vector of number of overlaps to subject weigthed
#'  by 'weight' column.
#' @family features
#' @export
#' @examples
#' gr1 <- GRanges(seqnames="chr1",
#'                ranges=IRanges(start = c(4, 9, 10, 30),
#'                               end = c(4, 15, 20, 31)),
#'                strand="+")
#' gr2 <- GRanges(seqnames="chr1",
#'                ranges=IRanges(start = c(1, 4, 15, 25),
#'                               end = c(2, 4, 20, 26)),
#'                strand=c("+"),
#'                score=c(10, 20, 15, 5))
#' countOverlaps(gr1, gr2)
#' countOverlapsW(gr1, gr2, weight = "score")
countOverlapsW <- function(query, subject, weight = NULL, ...) {
  if (is.null(weight)) return(countOverlaps(query, subject, ...))

  weight <- getWeights(subject, weight)
  hits = as(findOverlaps(query, subject, ...), "List")
  weightedCount = sum(extractList(weight, hits))

  names(weightedCount) <- names(query)
  return(weightedCount)
}




#' Get read hits per codon
#'
#' Helper for entropy function, normally not used directly
#' Seperate each group into tuples (abstract codons)
#' Gives sum for each tuple within each group
#'
#' Example: counts c(1,0,0,1), with reg_len = 2, gives
#' c(1,0) and c(0,1), these are summed and returned as data.table
#' 10 bases, will give 3 codons, 1 base codons does not exist.
#' @inheritParams scaledWindowPositions
#' @return a data.table with codon sums
#'
codonSumsPerGroup <- function(grl, reads, weight = "score", is.sorted = FALSE) {
  dt <- scaledWindowPositions(grl, reads, numCodons(grl), weight = weight,
                              is.sorted = is.sorted)
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
#' @inheritParams findFa
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

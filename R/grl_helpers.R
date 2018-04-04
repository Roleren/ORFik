#' Restrict GRangesList
#'
#' Will restrict GRangesList to \code{N} bp downstream from the first base.
#' @param grl (GRangesList)
#' @param firstN (integer) Allow only this many bp downstream
#' @return a GRangesList of reads restricted to firstN and tiled by 1
#'
downstreamN <- function(grl, firstN = 150L) {
  grl <- tile1(groupGRangesBy(unlist(grl, use.names = TRUE)))
  return(phead(grl, firstN))
}


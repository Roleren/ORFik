#' Get logical list of strands
#'
#' Helper function to get a logical list of True/False,
#'  if GRangesList group have + strand = T, if - strand = F
#' Also checks for * strands, so a good check for bugs
#' @param grl a \code{\link{GRangesList}} or GRanges object
#' @return a logical vector
#' @export
#' @examples
#' gr <- GRanges(Rle(c("chr2", "chr2", "chr1", "chr3"), c(1, 3, 2, 4)),
#'               IRanges(1:10, width = 10:1),
#'               Rle(strand(c("-", "+", "*", "+", "-")), c(1, 2, 2, 3, 2)))
#' strandBool(gr)
#'
strandBool <- function(grl) {
  if (is(grl, "GRanges")) {
    posIndices <- as.character(strand(grl)) == "+"
  } else {
    posIndices <- strandPerGroup(grl, FALSE) == "+"
  }

  sums <- sum(posIndices) + sum(!posIndices)
  if (is.na(sums)) {
    stop("could not get strands from grl object",
         " most likely NULL object was passed.")
  }
  if (sums != length(grl)) {
    stop("grl contains * strands, set them to either + or -")
  }
  return(posIndices)
}


#' Match naming of GRangesList
#'
#' Given a GRangesList and a reference, make the naming convention and
#' the number of metacolumns equal to reference
#' @param gr a \code{\link{GRangesList}}
#'  or GRanges object
#' @param reference a GRangesList of a reference
#' @return a GRangesList
#'
matchNaming <- function(gr, reference) {
  if (is.grl(gr)) gr <- unlistGrl(gr)

  ## now get a reference
  grTest <- unlist(reference[1], use.names = FALSE)

  # TODO: This can still be optimized for strange cases.
  # One case is that you can keep, even though ncol new > ncol old,
  # if all are equal within group
  mcols(gr) <- DataFrame(row.names = names(gr),
                         rep(mcols(grTest[1]), length(gr)))

  if (is.ORF(grTest) | any(grep(names(reference[1]), pattern = "_"))) {
    return(groupGRangesBy(gr, gr$names))
  } else {
    return(groupGRangesBy(gr))
  }
}

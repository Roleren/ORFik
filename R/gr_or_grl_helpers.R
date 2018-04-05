#' Get logical list of strands
#'
#' Helper function to get a logical list of True/False, if GRangesList group have
#' + strand = T, if - strand = F
#' Also checks for * strands, so a good check for bugs
#' @param grl a \code{\link[GenomicRanges]{GRangesList}} or GRanges object
#' @export
#' @return a logical vector
#'
strandBool <- function(grl){
  if (class(grl) == "GRanges") {
    posIndices <- as.character(strand(grl)) == "+"
  } else {
    posIndices <- strandPerGroup(grl, FALSE) == "+"
  }

  sums <- sum(posIndices) + sum(!posIndices)
  if (is.na(sums)) {
    stop("could not get strands from grl object,\n
          most likely NULL object was passed.")
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
#' @param gr a \code{\link[GenomicRanges]{GRangesList}}
#'  or GRanges object
#' @param reference a GRangesList of a reference
#' @return a GRangesList
matchNaming <- function(gr, reference){
  ## First check if unlist should be T or F
  if (is.grl(gr))
    gr <- unlistGrl(gr)

  ## now get a reference
  grTest <- unlist(reference[1], use.names = FALSE)

  ## TODO: add possibily to add unknown columns, i.g. exon_rank etc.
  ## Will try this now.
  if(ncol(elementMetadata(grTest)) == 1 &&
     colnames(elementMetadata(grTest)) == "names") {

  } else {
    equalMetaCols <- ncol(elementMetadata(grTest)) ==
      ncol(elementMetadata(gr))
    if (equalMetaCols) {
      equalMetaCols <- all(colnames(elementMetadata(grTest)) ==
                             colnames(elementMetadata(gr)))
    }

    if (!equalMetaCols) {
      refMeta <- elementMetadata(unlist(reference, use.names = FALSE))
      grMeta <- elementMetadata(gr)
      # if same number of elements just replace
      if (nrow(refMeta) == nrow(grMeta)) {
        elementMetadata(gr) <- refMeta
      } else {
        # else set them to NA
        df <- S4Vectors::DataFrame(matrix(NA,
                ncol = ncol(refMeta), nrow = length(gr)))

        colnames(df) <- colnames(refMeta)
        # set col classes
        for(i in 1:ncol(df)){
          class(df[,i]) <- class(refMeta[,i])
        }

        elementMetadata(gr) <- df
      }
    }
  }
  ## Name column is special case:
  ## add names column if reference have it
  if (!is.null(grTest$names)) gr$names <- names(gr)
  ## if reference have names, add them to gr
  if (!is.null(names(grTest))) {
    if (!any(grep(pattern = "_", names(grTest)[1]))) {
      names(gr) <- txNames(gr)
    }
  }

  if (!is.null(gr$names)) {
    return(groupGRangesBy(gr, gr$names))
  } else {
    return(groupGRangesBy(gr))
  }
}

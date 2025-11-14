#' Helper function to find overlaping seqlevels
#'
#' Keep only seqnames in reads that are in grl
#' Useful to avoid seqname warnings in bioC
#' @param grl a \code{\link{GRangesList}} or GRanges object
#' @param reads a GRanges, GAlignment or GAlignmentPairs object
#' @return a character vector of valid seqlevels
#' @family validity
#' @keywords internal
validSeqlevels <- function(grl, reads) {
  readNames <- unique(seqnames(reads))
  seqMatch <- readNames %in%
    unique(seqnamesPerGroup(grl, FALSE))
  return(readNames[seqMatch])
}

#' Helper function to check for GRangesList
#' @param class the class you want to check if is GRL,
#' either a character from class or the object itself.
#' @return a boolean
#' @family validity
#' @keywords internal
is.grl <- function(class) {
  if (!is.character(class)) {
    class <- class(class)
  }
  return((class == "GRangesList" || class == "CompressedGRangesList"))
}

is.rl <- function(class) {
  if (!is.character(class)) {
    class <- class(class)
  }
  return(is.grl(class) ||
           (class == "IRangesList" || class == "CompressedIRangesList"))
}


#' Helper function to check for GRangesList or GRanges class
#' @param class the class you want to check if is GRL or GR,
#'  either a character from class or the object itself.
#' @return a boolean
#' @family validity
#' @keywords internal
is.gr_or_grl <- function(class) {
  if (!is.character(class)) {
    class <- class(class)
  }
  return(is.grl(class) || class == "GRanges")
}

#' Helper function to check for ranged object
#' @param x the object to check is a ranged object.
#' Either GRangesList, GRanges, IRangesList, IRanges.
#' @return a boolean
#' @family validity
#' @keywords internal
is.range <- function(x) {
  return(is.gr_or_grl(x) | is(x, "IRanges") |
           is(x, "IRangesList"))
}

#' Check if all requirements for an ORFik ORF is accepted.
#' @param grl a GRangesList or GRanges to check
#' @return a logical (TRUE/FALSE)
#' @family validity
#' @keywords internal
is.ORF <- function(grl) {
  if (is.gr_or_grl(class(grl))){
    if (is.grl(grl)) {
      names <- unlist(grl[1], use.names = FALSE)$names
    } else {
      names <-grl[1]$names
    }
    return(any(grep("_", names)))
  }
  return(FALSE)
}

#' Helper Function to check valid GRangesList input
#' @param class as character vector the given class of
#'  supposed GRangesList object
#' @param type a character vector, is it gtf, cds, 5', 3', for messages.
#' @param checkNULL should NULL classes be checked and return indeces of these?
#' @return either NULL or indices (checkNULL == TRUE)
#' @family validity
#' @keywords internal
validGRL <- function(class, type = "grl", checkNULL = FALSE) {
  if(length(class) != length(type)) stop("not equal length of classes",
                                         " and types, see validGRL")
  if (checkNULL) {
    indeces <- "NULL" == class
    class <- class[!indeces]
    if (length(class) == 0) return(rep(TRUE, length(type)))
    type <- type[!indeces]
  }
  for (classI in seq_along(class)) {
    if (!is.grl(class[classI])) {
      messageI <- paste(type[classI], "must be given and be type GRangesList")
      stop(messageI)
    }
  }
  if (checkNULL) {
    return(indeces)
  }
}

validRL <- function(class, type = "rl", checkNULL = FALSE) {
  if(length(class) != length(type)) stop("not equal length of classes",
                                         " and types, see validsRL")
  if (checkNULL) {
    indeces <- "NULL" == class
    class <- class[!indeces]
    if (length(class) == 0) return(rep(TRUE, length(type)))
    type <- type[!indeces]
  }
  for (classI in seq_along(class)) {
    if (!is.rl(class[classI])) {
      messageI <- paste(type[classI], "must be given and be type GRangesList or IRangesList")
      stop(messageI)
    }
  }
  if (checkNULL) {
    return(indeces)
  }
}

#' Helper Function to check valid RNA input
#' @param class, the given class of RNA object
#' @return NULL, stop if unvalid object
#' @family validity
#' @keywords internal
checkRNA <- function(class){
  if (is.null(class) || (class == "NULL")) {
    message("No RNA added, skipping feature te and fpkm of RNA, ",
            "also ribosomeReleaseScore will also be not normalized best ",
            "way possible.")
  } else {
    if (class != "GAlignmentPairs" & class != "GAlignments" & class != "GRanges"
        & class != "covRle") {
      stop("RNA must be either of class GAlignmentPairs, GAlignments,",
           " covRle or GRanges")
    }
  }
}


#' Helper Function to check valid RFP input
#' @param class, the given class of RFP object
#' @return NULL, stop if invalid object
#' @family validity
#' @keywords internal
checkRFP <- function(class) {
  if (class != "GAlignments" & class != "GRanges" & class != "covRle") {
    stop("RFP must be either of class GAlignments, covRle or GRanges")
  }
}

#' A paste function for directories
#' Makes sure slashes are corrected, and not doubled.
#' @param ... any amount of arguments that are possible to convert
#'  to characters
#' @return the pasted string
#' @keywords internal
pasteDir <- function(...) {
  temp <- gsub(pattern = "//", x = paste(..., sep = "/"), replacement = "/")
  return(gsub(pattern = "/\\.", temp, replacement = "."))
}

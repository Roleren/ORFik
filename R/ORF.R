hasRanges <- function(object) {
    errors <- character()

    if (length(object@orfRanges) == 0) {
        errors <- c(errors, "orfRanges must have some IRanges in it")
    }
    if (length(errors) == 0) {
        return(TRUE)
    } else {
        return(errors)
    }
}

#' An S4 class to represent an ORF.
#'
#' @field orfRanges A GRanges of ORF
#' @field trailerRanges A GRanges of trailer
#' @field family character. Family of ORF e.g. 'lincRNA'
#' @field isoform character. What is isoform of this ORF toward the transcript
#' @field frame numeric.
#' @field ORFLength numeric. Length of ORF.
#' @field trailerLength numeric. Length of trailer.
#' @field transcript character. Name of the transcript.
#' @field start character. Defined start of the ORF. e.g. 'ATG'
#' @field identifier character. Unique identifier of this ORF.
#' @export
#' @import GenomicRanges
#' @import methods

setClass("ORF",
         representation(orfRanges = "GRanges",
                        trailerRanges = "GRanges",
                        family = "character",
                        isoform = "character",
                        frame = "numeric",
                        ORFlength = "numeric",
                        trailerLength = "numeric",
                        transcript = "character",
                        start = "character",
                        identifier = "character"),
         prototype(trailerRanges = GRanges(),
                   family = NA_character_,
                   isoform = NA_character_,
                   frame = NA_integer_,
                   transcript = NA_character_,
                   start = NA_character_),
         validity = hasRanges)

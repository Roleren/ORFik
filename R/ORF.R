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

#' An S4 class to represent a bank account.
#'
#' @field orfRanges A GRanges of ORF
#' @field trailerRanges A GRanges of trailer
#' @field family character. Family of ORF e.g. 'lincRNA'
#' @field ORFLength numeric. Length of ORF.
#' @field trailerLength numeric. Length of trailer.
#' @field transcript character. Name of the transcript.
#' @field start character. Defined start of the ORF. e.g. 'ATG'
#' @export
#' @import GenomicRanges
#' @import methods

setClass("ORF", representation(orfRanges = "GRanges", trailerRanges = "GRanges", family = "character", ORFlength = "numeric", 
    trailerLength = "numeric", transcript = "character", start = "character"), prototype(trailerRanges = GRanges(), family = NA_character_, 
    transcript = NA_character_, start = NA_character_), validity = hasRanges)

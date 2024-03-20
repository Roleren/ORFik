#' ORFik for analysis of open reading frames.
#'
#' Main goals:
#' \enumerate{
#' \item Finding Open Reading Frames (very fast) in the genome of interest or
#'  on the set of transcripts/sequences.
#' \item Utilities for metaplots of RiboSeq coverage over gene START and STOP
#' codons allowing to spot the shift.
#' \item Shifting functions for the RiboSeq data.
#' \item Finding new Transcription Start Sites with the use of CageSeq data.
#' \item Various measurements of gene identity e.g. FLOSS, coverage, ORFscore,
#' entropy that are recreated based on many scientific publications.
#' \item Utility functions to extend GenomicRanges for faster grouping,
#' splitting, tiling etc.
#' }
#'
#' @useDynLib ORFik
#' @import GenomicRanges GenomicAlignments GenomeInfoDb
#' @importFrom Rcpp sourceCpp
#'
"_PACKAGE"

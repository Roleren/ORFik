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

.onAttach <- function(libname, pkgname) {
  if (requireNamespace("qs2", quietly = TRUE)) {
    if (isFALSE(qs2:::check_TBB())) {
      packageStartupMessage(
        "qs2 detected without TBB support.\n",
        "For best performance, reinstall with:\n",
        "remotes::install_cran(\"qs2\", type = \"source\", ",
        "configure.args = \"--with-TBB --with-simd=AVX2\", ",
        "force = TRUE)"
      )
    }
  }
}

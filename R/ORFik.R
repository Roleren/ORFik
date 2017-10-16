#' ORFik for analysis of open reading frames.
#'
#' Main goals:
#' \enumerate{
#' \item Find open reading frames
#' \item Compatible with GRanges
#' \item Precise quantification of mutation rates.
#' \item Prepare automatic reports as .Rmd files that are flexible
#' and open for manipulation.
#' \item Provide specialized plots for deletions, insertions, mismatches,
#' variants, heterogeneity of the reads.
#' }
#'
#' To learn more about amplican, start with the vignettes:
#' \code{browseVignettes(package = "amplican")}
#'
#' @docType package
#' @name amplican
#'
#' @useDynLib ORFik
#' @importFrom Rcpp sourceCpp
#'
#' @import IRanges
"_PACKAGE"


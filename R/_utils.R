#' Source bioconductor
#'
#' Helper function for quick update of bioconductor packages,
#' @param packages either NULL if only source and no update/install
#' or "all" if you want to update all your bioconductor packages
#' or c(package1, package2, ...)
#' for specific packages as a character vector
#' @name sourceBioc
#' @return NULL
#'
sourceBioc <- function(packages = NULL) {
  source("https://bioconductor.org/biocLite.R")
  if (!is.null(packages)) {
    biocLite <- NULL # Get this from bioc
    if (packages == "all") {
      biocLite()
    } else {
      biocLite(packages)
    }
  }
}

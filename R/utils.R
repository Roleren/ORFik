#' Converts different type of files to Granges
#'
#' column 5 will be set to score
#' Only Accepts bed files for now, standard format from Fantom5
#' @param x An data.frame from imported bed-file,
#'  to convert to GRanges
#' @param bed6 If bed6, no meta column is added
#' @return a GRanges object from bed
bedToGR <- function(x, bed6 = TRUE){

  if (!bed6) {
    gr <- GRanges(x[, 1], IRanges(x[, 2] + 1, x[, 3]))
    return(gr)
  }
  starts <- x[, 2] + 1
  ends <- x[, 3]
  gr <- GRanges(x[, 1], IRanges(starts, ends),
                strand = x[, 6])
  score(gr) <- x[, 5]
  if (ncol(x) > 6) mcols(gr) <- x[, 7:ncol(x)]
  return(gr)
}

#' Get bed file from a file-path
#'
#' Tries to speed up rtracklayer import.bed
#' If speedup is not supported, it will use normal
#' rtracklayer::import.bed
#' Supports gzip, gz, bgz and bed formats
#' @param filePath The location of the bed file
#' @importFrom data.table fread setDF
#' @importFrom tools file_ext
#' @importFrom rtracklayer import.bed
#' @return a GRanges object
#' @export
#' @examples
#' # path to example CageSeq data from hg19 heart sample
#' cageData <- system.file("extdata", "cage_data_heart.bed.bgz",
#'                        package = "ORFik")
#'
#' fread.bed(cageData)
#'
fread.bed <- function(filePath) {

  if (.Platform$OS.type == "unix") {
    if (file.exists(filePath)) {
      if (any(file_ext(filePath) == c("gzip", "gz", "bgz"))) {
        bed <- bedToGR(setDF(
          fread(paste("gunzip -c", filePath), sep = "\t")))
      } else if (file_ext(filePath) == "bed"){
        bed <- bedToGR(setDF(fread(filePath, sep = "\t")))
      } else {
        bed <- import.bed(con =  filePath)
      }
    } else {stop("Filepath specified does not name existing file.")}
  } else {
    ## NB: Windows user will have slower loading
    bed <- import.bed(con =  filePath)
  }

  return(bed)
}

#' Source bioconductor
#'
#' Helper function for quick update of bioconductor packages,
#' @param packages either NULL if only source and no update/install
#' or "all" if you want to update all your bioconductor packages
#' or c(package1, package2, ...)
#' for specific packages as a character vector
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

#' Helper function to check for GRangesList
#' @param class the class you want to check if is GRL,
#' either a character from class or the object itself.
#' @return a boolean
#'
is.grl <- function(class) {
  if (!is.character(class)) {
    class <- class(class)
  }
  return((class == "GRangesList" || class == "CompressedGRangesList"))
}


#' Helper function to check for GRangesList or GRanges class
#' @param class the class you want to check if is GRL or GR,
#'  either a character from class or the object itself.
#' @return a boolean
#'
is.gr_or_grl <- function(class) {
  if (!is.character(class)) {
    class <- class(class)
  }
  return(is.grl(class) || class == "GRanges")
}


#' Helper Function to check valid GRangesList input
#' @param class as character vector the given class of
#'  supposed GRangesList object
#' @param type a character vector, is it gtf, cds, 5', 3', for messages.
#' @param checkNULL should NULL classes be checked and return indeces of these?
#' @return either NULL or indices (checkNULL == TRUE)
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


#' Converts different type of files to Granges
#'
#' column 5 will be set to score
#' Only Accepts bed files for now, standard format from Fantom5
#' @param x An data.frame from imported bed-file,
#'  to convert to GRanges
#' @param bed6 If bed6, no meta column is added
#' @return a GRanges object from bed
#'
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


#' Load bed file as GRanges.
#'
#' Wraps around rtracklayer::import.bed and tries to speed up loading with the
#' use of data.table. Supports gzip, gz, bgz and bed formats.
#' @param filePath The location of the bed file
#' @importFrom data.table fread setDF
#' @importFrom tools file_ext
#' @importFrom rtracklayer import.bed
#' @return a GRanges object
#' @export
#' @examples
#' # path to example CageSeq data from hg19 heart sample
#' cageData <- system.file("extdata", "cage-seq-heart.bed.bgz",
#'                         package = "ORFik")
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
#' Convenience function for quick update of bioconductor packages,
#' @param packages either NULL if only source and no update/install
#' or "all" if you want to update all your bioconductor packages
#' or c(package1, package2, ...) for specific packages as a character vector
#' @return NULL, loading the packages requested, or just source
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
#' Convenience wrapper for Rsamtools FaFile
#' @param faFile a character path or FaFile
#' @importFrom Rsamtools FaFile
#' @return a FaFile or BSgenome
findFa <- function(faFile) {
  if (is.character(faFile)) {
    if (dir.exists(faFile)) {
      return(FaFile(faFile))
    } else {
      stop("faFile does not name a valid fasta file")
    }

  } else if (class(faFile) == "FaFile" || class(faFile) == "BSgenome") {
    return(faFile)
  }
  stop("faFile must be FaFile, BSgenome or valid filePath")
}

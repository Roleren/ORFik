#' Load bed file as GRanges.
#'
#' Wraps around rtracklayer::import.bed and tries to speed up loading with the
#' use of data.table. Supports gzip, gz, bgz and bed formats.
#' Also safer chromosome naming with the argument chrStyle
#' @param filePath The location of the bed file
#' @inheritParams matchSeqStyle
#' @importFrom data.table fread setDF
#' @importFrom tools file_ext
#' @importFrom rtracklayer import.bed
#' @return a \code{\link{GRanges}} object
#' @export
#' @family utils
#' @examples
#' # path to example CageSeq data from hg19 heart sample
#' cageData <- system.file("extdata", "cage-seq-heart.bed.bgz",
#'                         package = "ORFik")
#' fread.bed(cageData)
#'
fread.bed <- function(filePath, chrStyle = NULL) {

  if (.Platform$OS.type == "unix") {
    if (file.exists(filePath)) {
      if (any(file_ext(filePath) %in% c("gzip", "gz", "bgz"))) {
        bed <- bedToGR(setDF(
          fread(cmd = paste("gunzip -c", filePath), sep = "\t")))
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

  return(matchSeqStyle(bed, chrStyle))
}

#' Custom bam reader
#'
#' Safer version that handles the most important error done.
#' In the future will use a faster .bam loader for big .bam files in R.
#' @param path a character path to .bam file. If paired end bam files,
#' input must be a data.table with two columns (forward and reverse)
#' and one row:\cr
#' if paired end reads in single bam file:\cr
#' forward contains paired end bam file path, reverse must be
#' either "paired-end" or "" (single end).\cr
#' if paired end reads split in two files:\cr
#' forward contains paired end bam file path (R1), reverse must be
#' paired end bam file path (R2 file), this is a rare case\cr
#' If all are single-end or
#' you don't need to load data as paired end, the reverse column
#' can be skipped.
#' @inheritParams matchSeqStyle
#' @return a \code{\link{GAlignments}} object of bam file
#' @export
#' @family utils
#' @examples
#' bam_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' readBam(bam_file, "UCSC")
readBam <- function(path, chrStyle = NULL) {
  if (!(length(path) %in% c(1,2))) stop("readBam must have 1 or 2 bam files!")
  if (is(path, "factor")) path <- as.character(path)
  # If data.table path
  if (is(path, "data.table")) {
    if (path$reverse == "paired-end") {
      message("ORFik reads this paired end bam as readGAlignmentPairs")
      bam <- matchSeqStyle(readGAlignmentPairs(path$forward), chrStyle)
      if (length(bam) == 0)
        stop(paste("File", path$forward,
                   "was read as one paired-end file, but had 0 paired reads!"))
      return(bam)
    } else {
      message("ORFik reads these split paired end bams as readGAlignments combination")
      return(matchSeqStyle(c(readGAlignments(path$forward),
                             readGAlignments(path$reverse)), chrStyle))
    }
  }
  # If character path
  if (is(path, "character") & length(path) == 2) {
    if (path[2] == "paired-end"){
      message("ORFik reads paired end bam in as readGAlignmentPairs")
      bam <- matchSeqStyle(readGAlignmentPairs(path[1]), chrStyle)
      if (length(bam) == 0)
        stop(paste("File", path[1],
                   "was read as paired-end file, but had 0 paired reads!"))
      return(bam)
    } else {
      message("ORFik reads these split paired end bams as readGAlignments combination")
      return(matchSeqStyle(c(readGAlignments(path[1]),
                             readGAlignments(path[2])), chrStyle))
    }
  }
  # else single end bam file
  return(matchSeqStyle(readGAlignments(path), chrStyle))
}

#' Custom wig reader
#'
#' Given 2 wig files, first is forward second is reverse.
#' Merge them and return as GRanges object.
#' If they contain name reverse and forward, first and second order
#' does not matter, it will search for forward and reverse.
#'
#' @param path a character path to two .wig files, or a data.table
#' with 2 columns, (forward, filepath) and reverse, only 1 row.
#' @inheritParams matchSeqStyle
#' @importFrom rtracklayer import.wig
#' @return a \code{\link{GRanges}} object of the file/s
#' @family utils
#'
readWig <- function(path, chrStyle = NULL) {
  if (is(path, "character")) {
    if (length(path) != 2) stop("readWig must have 2 wig files,
                              one forward strand and one reverse!")

    forwardPath <- grep("forward|fwd", path)
    reversePath <- grep("reverse|rev", path)
    if (length(forwardPath) == 1 & length(reversePath) == 1){
      forwardIndex <- forwardPath
      reverseIndex <- reversePath
    }

    forward <- import.wig(path[forwardIndex])
    reverse <- import.wig(path[reverseIndex])
  } else if (is(path, "data.table")) {
    if (!is.null(path$forward)) {
      forward <- import.wig(path$forward)
    } else forward <- import.wig(path$filepath)
    reverse <- import.wig(path$reverse)
  }
  strand(forward) <- "+"
  strand(reverse) <- "-"
  return(matchSeqStyle(c(forward, reverse), chrStyle))
}

#' Load GRanges object from .bedo
#'
#' .bedo is .bed ORFik, an optimized bed format for coverage reads with read lengths
#' .bedo is a text based format with columns (6 maximum):\cr
#' 1. chromosome\cr 2. start\cr 3. end\cr 4. strand\cr
#' 5. ref width (cigar # M's, match/mismatch total)\cr
#' 6. duplicates of that read\cr
#'
#' Positions are 1-based, not 0-based as .bed.
#' export with export.bedo
#' @param path a character, location on disc (full path)
#' @return GRanges object
#' @importFrom tools file_ext
#' @export
import.bedo <- function(path) {
  if (file_ext(path) != "bedo") stop("import.bedo can only load .bedo files!")
  dt <- fread(input = path)
  if (is.null(dt$end)) dt$end <- dt$start
  return(makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE))
}

#' Load GAlignments object from .bedoc
#'
#' A much faster way to store, load and use bam files.\cr
#' .bedoc is .bed ORFik, an optimized bed format for coverage reads with
#' cigar and replicate number.\cr
#' .bedoc is a text based format with columns (5 maximum):\cr
#' 1. chromosome\cr 2. cigar: (cigar # M's, match/mismatch total)
#' \cr 3. start (left most position) \cr 4. strand (+, -, *)\cr
#' 5. score: duplicates of that read\cr
#'
#' Positions are 1-based, not 0-based as .bed.
#' export with export.bedo
#' @param path a character, location on disc (full path)
#' @return GAlignments object
#' @importFrom tools file_ext
#' @export
import.bedoc <- function(path) {
  if (file_ext(path) != "bedoc")
    stop("import.bedoc can only load .bedoc files!")
  return(makeGAlignmentsFromDataFrame(fread(input = path)))
}

#' #' Load GRanges / GAlignments object from .ofst
#'
#' A much faster way to store, load and use bam files.\cr
#' .ofst is ORFik fast serialized object,
#' an optimized bed format for coverage reads with
#' cigar and replicate number.\cr
#' .ofst is a text based format with minimum 4 columns:\cr
#' 1. chromosome\cr  2. start (left most position) \cr 3. strand (+, -, *)\cr
#' 4. width (not added if cigar exists)\cr
#' 5. cigar (not needed if width exists):
#'  (cigar # M's, match/mismatch total) \cr
#' 5. score: duplicates of that read\cr
#' 6. size: qwidth according to reference of read
#' Other columns can be named whatever you want and added to meta columns.
#' Positions are 1-based, not 0-based as .bed.
#' Import with import.ofst
#' @param file a path to a .ofst file
#' @return a GAlignment or GRanges object, dependent of if cigar is
#' defined in .ofst file.
#' @importFrom fst read_fst
#' @export
import.ofst <- function(file) {
  df <- read_fst(file)
  if ("cigar" %in% colnames(df)) {
    getGAlignments(df)
  } else getGRanges(df)
}

#' Load any type of sequencing reads
#'
#' Wraps around rtracklayer::import and tries to speed up loading with the
#' use of data.table. Supports gzip, gz, bgz compression formats.
#' Also safer chromosome naming with the argument chrStyle
#'
#' NOTE: For wig you can send in 2 files, so that it automaticly merges
#' forward and reverse stranded objects. You can also just send 1 wig file,
#' it will then have "*" as strand.
#'
#' @param path a character path to file (1 or 2 files),
#'  or data.table with 2 colums(forward&reverse)
#'  or a GRanges/Galignment/GAlignmentPairs object etc.
#'  If it is ranged object it will presume to be
#'  already loaded, so will return the object as it is,
#'  updating the seqlevelsStyle if given.
#' @inheritParams matchSeqStyle
#' @importFrom tools file_ext
#' @importFrom tools file_path_sans_ext
#' @importFrom rtracklayer import
#' @return a \code{\link{GAlignments}}/\code{\link{GRanges}} object,
#'  depending on input.
#' @export
#' @family utils
#' @examples
#' bam_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' fimport(bam_file)
#' # Certain chromosome naming
#' fimport(bam_file, "NCBI")
#'
fimport <- function(path, chrStyle = NULL) {
  if (is(path, "data.table")) {
    if (ncol(path) == 2 & colnames(path) == c("forward", "reverse")) {
      path <- c(path$forward, path$reverse)
    } else stop("When path is data.table,",
                "it must have 2 columns (forward&reverse)")
  } else if (is(path, "list")) {
    if (!(all(lengths(path) %in% c(1,2))) | length(path) != 1) {
      stop("When path is list,",
           "it must have 1 or 2 elements (forward, reverse)")
    }
  }
  if (is(path, "list") | is(path, "data.table") | is(path, "character")) {
    path <- unlist(path)
    path <- path[path != ""]
    pairedEndBam <- FALSE
    if (length(path) == 2) {
      if (path[2] == "paired-end") {
        pairedEndBam <- TRUE
        path <- path[1]
      }
    }
  }

  if (is.character(path)) {
    if (all(file.exists(path))) {
      fext <- file_ext(path)
      compressions <- c("gzip", "gz", "bgz", "zip")
      areCompressed <- fext %in% compressions
      fext[areCompressed] <- file_ext(file_path_sans_ext(path[areCompressed],
                                                         compression = FALSE))
      if (length(path) == 2) { # Multiple file paths
        if (all(fext %in% c("wig"))) {
          return(readWig(path, chrStyle))
        } else if (all(fext %in% c("bam"))) {
          return(readBam(path, chrStyle))
        } else stop("only wig and valid bam format allowed for 2 files input!")
      } else if (length(path) == 1) { # Only 1 file path given
        if (fext == "bam") {
          if (pairedEndBam) path <- c(path, "paired-end")
          return(readBam(path, chrStyle))
        } else if (fext == "bed" |
                   file_ext(file_path_sans_ext(path,
                                               compression = TRUE)) == "bed" |
                   file_ext(file_path_sans_ext(path, compression = FALSE))
                   == "bed") {
          return(fread.bed(path, chrStyle))
        } else if (fext == "bedo") {
          return(matchSeqStyle(import.bedo(path), chrStyle))
        } else if (fext == "bedoc") {
          return(matchSeqStyle(import.bedoc(path), chrStyle))
        } else if (fext == "ofst") {
          return(matchSeqStyle(import.ofst(path), chrStyle))
        }else return(matchSeqStyle(import(path), chrStyle))
      } else stop("fimport takes either 1 or 2 files!")
    } else stop(paste(path, "does not exist as File/Files!"))
  } else if (is.gr_or_grl(path) | is(path, "GAlignments") |
             is(path, "GAlignmentPairs")) {
    return(matchSeqStyle(path, chrStyle))
  } else {
    stop("path must be either a valid character",
         " filepath or ranged object.")
  }
}

#' Convenience wrapper for Rsamtools FaFile
#'
#' Get fasta file object, to find sequences in file.\cr
#' Will load and import file if necessarry.
#' @param faFile \code{\link{FaFile}}, BSgenome, fasta/index file path or an
#' ORFik \code{\link{experiment}}. This file is usually used to find the
#' transcript sequences from some GRangesList.
#' @importFrom Rsamtools FaFile
#' @importFrom methods is
#' @return a \code{\link{FaFile}} or BSgenome
#' @family utils
#' @export
#'
findFa <- function(faFile) {
  if (is.character(faFile)) {
    if (file.exists(faFile)) {
      return(FaFile(faFile))
    } else {
      stop("faFile does not name a valid fasta/index file")
    }
  } else if (is(faFile, "FaFile") || is(faFile, "BSgenome")) {
    return(faFile)
  } else if (is(faFile, "experiment")) {
    return(FaFile(faFile@fafile))
  }
  stop("faFile must be FaFile, BSgenome, valid filePath, or ORFik experiment")
}


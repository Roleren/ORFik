
#' Converts different type of files to Granges
#'
#' column 5 will be set to score
#' Only Accepts bed files for now, standard format from Fantom5
#' @param x A \code{\link{data.frame}} from imported bed-file,
#'  to convert to GRanges
#' @param bed6 If bed6, no meta column is added
#' @return a \code{\link{GRanges}} object from bed
#' @family utils
#'
bedToGR <- function(x, bed6 = TRUE) {

  if (!bed6) {
    gr <- GRanges(x[, 1L], IRanges(x[, 2L] + 1L, x[, 3L]))
    return(gr)
  }
  starts <- x[, 2L] + 1L
  ends <- x[, 3L]
  gr <- GRanges(x[, 1L], IRanges(starts, ends),
                strand = x[, 6L])
  mcols(gr) <- S4Vectors::DataFrame(mcols(gr), score = x[, 5L])
  if (ncol(x) > 6L) mcols(gr) <- x[, 7L:ncol(x)]
  return(gr)
}

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
#' Only for single end reads
#' Safer version that handles the most important error done.
#' In the future will use a faster .bam loader for big .bam files in R.
#' @param path a character path to .bam file
#' @inheritParams matchSeqStyle
#' @return a \code{\link{GAlignments}} object of bam file
#' @export
#' @family utils
#' @examples
#' bam_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' readBam(bam_file, "UCSC")
readBam <- function(path, chrStyle = NULL) {
  if (is(path, "factor")) path <- as.character(path)
  return(matchSeqStyle(readGAlignments(path), chrStyle))
}

#' Custom wig reader
#'
#' Given 2 wig files, first is forward second is reverse.
#' Merge them and return as GRanges object.
#' If they contain name reverse and forward, first and second order
#' does not matter, it will search for forward and reverse.
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
#' Find pair of forward and reverse strand wig files
#' @param paths a character path at least one .wig file
#' @return if not all are paired, return original list,
#' if they are all paired, return a data.table with matches as 2 columns
findWigPairs <- function(paths) {
  forwardPath <- grep("forward\\.wig*|fwd\\.wig*", paths)
  reversePath <- grep("reverse\\.wig*|rev\\.wig*", paths)

  if ((length(forwardPath) != length(reversePath)) |
      length(forwardPath) == 0 | length(reversePath) == 0) return(paths)
  dt <- data.table(forward = paths[forwardPath],
                   reverse = paths[reversePath], match = FALSE)
  for (row in seq(nrow(dt))) {
    if (gsub(pattern = "forward\\.wig*|fwd\\.wig*", x = dt$forward[row], "")
        == gsub(pattern = "reverse\\.wig*|rev\\.wig*", x = dt$reverse[row], ""))
        dt$match[row] = TRUE
  }
  if (all(dt$match)) {
    return(dt)
  }
  return(paths)
}

#' Store GRanges object as .bedo
#'
#' .bedo is .bed ORFik, an optimized bed format for coverage reads with read lengths
#' .bedo is a text based format with columns:
#' chromosome, start, stop, width, strand, (cigar # M's, match/mismatch total)
#' , duplicates of that read
#' @param object a GRanges object
#' @param out a character, location on disc (full path)
#' @return NULL, object saved to disc
#' @importFrom data.table fwrite
#'
export.bedo <- function(object, out) {
  if (!is(object, "GRanges")) stop("object must be GRanges")
  dt <- setDT(as.data.frame(object))
  fwrite(dt, file = out)
}

#' Load GRanges object from .bedo
#'
#' .bedo is .bed ORFik, an optimized bed format for coverage reads with read lengths
#' .bedo is a text based format with columns:
#' chromosome, start, stop, width, strand, (cigar # M's, match/mismatch total)
#' , duplicates of that read
#'
#' export with export.bedo
#' @param path a character, location on disc (full path)
#' @return GRanges object
#' @importFrom tools file_ext
#' @export
import.bedo <- function(path) {
  if (file_ext(path) != "bedo") stop("export.bedo can only load .bedo files!")
  return(makeGRangesFromDataFrame(fread(input = path), keep.extra.columns = TRUE))
}

#' Remove file extension of path
#'
#' Allows removal of compression
#' @param path character path (allows multiple paths)
#' @param basename relative path (TRUE) or full path (FALSE)? (default: FALSE)
#' @importFrom tools file_ext
#' @return character path without file extension
remove.file_ext <- function(path, basename = FALSE) {
  out <- c()
  for (p in path) {
    fext <- file_ext(path)
    compressions <- c("gzip", "gz", "bgz", "zip")
    areCompressed <- fext %in% compressions
    if (areCompressed) {
      ext <- file_ext(file_path_sans_ext(path, compression = FALSE))
      regex <- paste0("*\\.",ext,"\\.", fext,"$")
    } else {
      regex <- paste0("*\\.",fext,"$")
    }
    new <- gsub(pattern = regex, "", path)
    out <- c(out, new)
  }
  return(ifelse(basename, basename(out), out))
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
#' @param path a character path to file or a GRanges/Galignment object etc.
#' Any Ranged object.
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
  if (is.character(path)) {
    if (all(file.exists(path))) {
      fext <- file_ext(path)
      compressions <- c("gzip", "gz", "bgz", "zip")
      areCompressed <- fext %in% compressions
      fext[areCompressed] <- file_ext(file_path_sans_ext(path[areCompressed],
                                               compression = FALSE))
      if (length(path) > 1) { # Multiple file paths
        if (all(fext %in% c("wig"))) {
          return(readWig(path, chrStyle))
        } else if (all(fext %in% c("bam"))) {
          stop("only wig format allowed for multiple files!,
               is this paired end bam?")
        } else stop("only wig format allowed for multiple files!")
      } else { # Only 1 file path given
        if (fext == "bam") {
          return(readBam(path, chrStyle))
        } else if (fext == "bed" |
                   file_ext(file_path_sans_ext(path,
                                               compression = TRUE)) == "bed" |
                   file_ext(file_path_sans_ext(path, compression = FALSE))
                                                                  == "bed") {
          return(fread.bed(path, chrStyle))
        } else if (fext == "bedo") {
          return(matchSeqStyle(import.bedo(path), chrStyle))
        } else return(matchSeqStyle(import(path), chrStyle))
      }
    } else stop(paste0(path, "does not exist as File/Files!"))
  } else if (is.gr_or_grl(path) | is(path, "GAlignments")) {
    return(matchSeqStyle(path, chrStyle))
  } else {
    stop("path must be either a valid character",
         " filepath or ranged object.")
  }
}

#' A wrapper for seqlevelsStyle
#'
#' To make sure chromosome naming is correct (chr1 vs 1 vs I etc)
#' @param range a ranged object, (GRanges, GAlignment etc)
#' @param chrStyle a GRanges object, TxDb, FaFile,
#' or a \code{\link{seqlevelsStyle}}
#' (Default: NULL) to get seqlevelsStyle from. Is chromosome 1
#' called chr1 or 1, is mitocondrial chromosome called MT or chrM etc.
#' Will use 1st seqlevel-style if more are present.
#' Like: c("NCBI", "UCSC") -> pick "NCBI"
#' @return a GAlignment/GRanges object depending on input.
matchSeqStyle <- function(range, chrStyle = NULL) {
  # if needed add this ->
  # if (tryCatch(seqlevelsStyle(cage) <- seqlevelsStyle(fiveUTRs),
  #              error = function(e) {TRUE}) == TRUE) {
  #   warning("seqlevels of CAGE/fiveUTRs are not standardized, check them.")
  # } else {
  #   seqlevelsStyle(cage) <- seqlevelsStyle(fiveUTRs)
  # }
  if (!is.null(chrStyle)) {
    if (is.character(chrStyle)) {
      seqlevelsStyle(range) <- chrStyle[1]
    } else if (is.gr_or_grl(chrStyle) | is(chrStyle, "TxDb") |
               is(chrStyle, "FaFile")) {
      seqlevelsStyle(range) <- seqlevelsStyle(chrStyle)[1]
    } else stop("chrStyle must be valid GRanges object, or a valid chr style!")
  }
  return(range)
}

#' Find optimized subset of valid reads
#'
#' Keep only the ones that overlap within the grl ranges.
#' Also sort them in the end
#' @inheritParams validSeqlevels
#' @return the reads as GRanges or GAlignment
#' @family utils
#'
optimizeReads <- function(grl, reads) {
  seqMatch <- validSeqlevels(grl, reads)
  reads <- keepSeqlevels(reads, seqMatch, pruning.mode = "coarse")

  reads <- reads[countOverlaps(reads, grl, type = "within") > 0]
  reads <- sort(reads)
  if (length(reads) == 0) warning("No reads left in 'reads' after",
                                    "optimisation!")

  return(reads)
}

#' Convert a GRanges Object to 1 width reads
#'
#' There are 5 ways of doing this
#' 1. Take 5' ends, reduce away rest (5prime)
#' 2. Take 3' ends, reduce away rest (3prime)
#' 3. Tile to 1-mers and include all (tileAll)
#' 4. Take middle point per GRanges (middle)
#' 5. Get original with metacolumns (None)
#'
#' You can also do multiple at a time, then output is GRangesList, where
#' each list group is the operation (5prime is [1], 3prime is [2] etc)
#'
#'
#' Many other ways to do this have their own functions, like startSites and
#' stopSites etc.
#' To retain information on original width, set addSizeColumn to TRUE.
#' To compress data, 1 GRanges object per unique read, set addScoreColumn to
#' TRUE. This will give you a score column with how many duplicated reads there
#' were in the specified region.
#'
#' NOTE: Does not support paired end reads for the moment!
#' @param gr GRanges, GAlignment Object to reduce
#' @param method the method to reduce, see info. (5prime defualt)
#' @param addScoreColumn logical (FALSE), if TRUE, add a score column that
#'  sums up the hits per unique range This will make each read unique, so
#'  that each read is 1 time, and score column gives the number of hits.
#'  A useful compression. If addSizeColumn is FALSE, it will not differentiate
#'  between reads with same start and stop, but different length.
#' @param addSizeColumn logical (FALSE), if TRUE, add a size column that
#'  for each read, that gives original width of read. Useful if you need
#'  original read lengths. This takes care of soft clips etc.
#' @param after.softclips logical (TRUE), TRUE: include softclips in size
#' @return  Converted GRanges object
#' @export
#' @family utils
#'
convertToOneBasedRanges <- function(gr, method = "5prime",
                                    addScoreColumn = FALSE,
                                    addSizeColumn = FALSE,
                                    after.softclips = TRUE) {
  if (is(gr, "GAlignmentPairs")) stop("Paired end reads not supported,
                                      load as GAlignments instead!")

  if (addSizeColumn & is.null(mcols(gr)$size)) {
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr),
                                      size = readWidths(gr, after.softclips))
  }
  if (addScoreColumn) {
    dt <- data.table(seqnames = as.character(seqnames(gr)),
                     start = start(gr),
                     end = end(gr),
                     strand = as.character(strand(gr)))
    if (addSizeColumn) {
      dt[, size := mcols(gr)$size]
      dt <- dt[, .(score = .N), .(seqnames, start, end, strand, size)]
    } else {
      dt <- dt[, .(score = .N), .(seqnames, start, end, strand)]
    }
    gr <- makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE)
  }

  gr <- GRanges(gr)
  if (method == "5prime") {
    gr <- resize(gr, width = 1, fix = "start")
  } else if(method == "3prime") {
    gr <- resize(gr, width = 1, fix = "end")
  } else if(method == "None") {
  } else if(method == "tileAll") {
    gr <- unlist(tile(gr, width = 1), use.names = FALSE)
  } else if (method == "middle") {
    ranges(gr) <- IRanges(start(gr) + ceiling((end(gr) - start(gr)) / 2),
                          width = 1)
  } else stop("invalid type: must be 5prime, 3prime, None, tileAll or middle")

  return(gr)
}

#' Convenience wrapper for Rsamtools FaFile
#'
#' Get fasta file object, to find sequences in file.
#' @param faFile \code{\link{FaFile}}, BSgenome, fasta/index file path or an
#' ORFik \code{\link{experiment}}. This file used to find the
#' transcript sequences
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

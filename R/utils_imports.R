#' Load bed file as GRanges
#'
#' Wraps around \code{\link{import.bed}} and
#' tries to speed up loading with the
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
#' Read in Bam file from either single end or paired end.
#' Safer combined version of \code{\link{readGAlignments}} and
#' readGAlignmentPairs that takes care of some common errors.\cr
#' If QNAMES of the aligned reads are from collapsed fasta files
#' (if the names are formated from collapsing in either
#' (ORFik, ribotoolkit or fastx)), the
#' bam file will contain a meta column called "score" with the counts
#' of duplicates per read. Only works for single end reads, as perfect duplication
#' events for paired end is more rare.\cr
#'
#' In the future will use a faster .bam loader for big .bam files in R.
#' @param path a character / data.table with path to .bam file. There
#' are 3 input file possibilities.
#' \itemize{
#'  \item{single end : }{a character path (length 1)}
#'  \item{paired end (1 file) : }{Either a character path (length of 2), where
#'  path[2] is "paired-end", or a data.table
#'   with 2 columns, forward = path & reverse = "paired-end"}
#'  \item{paired end (2 files) : }{Either a character path (length of 2), where
#'  path[2] is path to R2, or a data.table
#'   with 2 columns, forward = path to R1 & reverse = path to R2.
#'   (This one is not used often)}
#' }
#' @inheritParams matchSeqStyle
#' @inheritParams GenomicAlignments::readGAlignments
#' @param strandMode numeric, default 0. Only used for paired end bam files.
#' One of (0: strand = *, 1: first read of pair is +, 2: first read of pair is -).
#' See ?strandMode. Note: Sets default to 0 instead of 1, as readGAlignmentPairs uses 1.
#' This is to guarantee hits, but will also make mismatches of overlapping
#' transcripts in opposite directions.
#' @return a \code{\link{GAlignments}} or \code{\link{GAlignmentPairs}} object of bam file
#' @importFrom Rsamtools scanBam BamFile ScanBamParam
#' @export
#' @family utils
#' @examples
#' bam_file <- system.file("extdata", "ribo-seq.bam", package = "ORFik")
#' readBam(bam_file, "UCSC")
readBam <- function(path, chrStyle = NULL, param = NULL, strandMode = 0) {
  if (!(length(path) %in% c(1,2))) stop("readBam must have 1 or 2 bam files!")
  if (is(path, "factor")) path <- as.character(path)
  # If data.table path
  if (is(path, "data.table")) {
    if (path$reverse == "paired-end") {
      message("ORFik reads this paired end bam as readGAlignmentPairs")
      message(paste("strandMode =", strandMode, ". Update it if is wrong!"))
      bam <- matchSeqStyle(readGAlignmentPairs(path$forward, param = param,
                                               strandMode = strandMode), chrStyle)
      if (length(bam) == 0)
        stop(paste("File", path$forward,
                   "was read as paired-end file, but had 0 paired reads!"))
      return(bam)
    } else {
      message("ORFik reads these split paired end bams as readGAlignments combination")
      return(matchSeqStyle(c(readGAlignments(path$forward, param = param),
                             readGAlignments(path$reverse, param = param)), chrStyle))
    }
  }
  # If character path
  if (is(path, "character") & length(path) == 2) {
    if (path[2] == "paired-end"){
      message("ORFik reads paired end bam in as readGAlignmentPairs")
      message(paste("strandMode =", strandMode, ". Update it if is wrong!"))
      bam <- matchSeqStyle(readGAlignmentPairs(path[1], param = param,
                                               strandMode = strandMode), chrStyle)
      if (length(bam) == 0)
        stop(paste("File", path[1],
                   "was read as paired-end file, but had 0 paired reads!"))
      return(bam)
    } else {
      message("ORFik reads these split paired end bams as readGAlignments combination")
      return(matchSeqStyle(c(readGAlignments(path[1], param = param),
                             readGAlignments(path[2], param = param)), chrStyle))
    }
  }
  # else single end bam file

  # Check if it is a collapsed reads format bam file
  # Get qnames with the data (check first 2 rows)
  headers <- unlist(scanBam(BamFile(path, yieldSize=2),
                            param = ScanBamParam(what = "qname")),
                    use.names = FALSE)
  # Check with and without extra ">" start sign
  bam.is.collapsed <- all(seq_along(headers) %in%
                            grep("(^(>|>seq)\\d+(-|_x)\\d+$)|(^(seq)\\d+(-|_x)\\d+$)", headers))
  if (bam.is.collapsed) {

    headers <- unlist(scanBam(path, param = ScanBamParam(what = "qname")),
                      use.names = FALSE)
    format.header <- ifelse(all(seq_along(headers) %in%
                                  grep("^>seq\\d+_x\\d+$|^seq\\d+_x\\d+$", headers)),
                            "ribotoolkit",
                            "fastx")
    if (format.header == "ribotoolkit") {
      scores <- as.integer(gsub(".*_x", "", headers))
    } else {
      scores <- as.integer(gsub(".*-", "", headers))
    }
    bam <- matchSeqStyle(readGAlignments(path, param = param), chrStyle)
    mcols(bam) <- DataFrame(score = scores)
    return(bam)
  } else return(matchSeqStyle(readGAlignments(path, param = param), chrStyle))
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

#' Custom bigWig reader
#'
#' Given 2 bigWig files (.bw, .bigWig), first is forward second is reverse.
#' Merge them and return as GRanges object.
#' If they contain name reverse and forward, first and second order
#' does not matter, it will search for forward and reverse.
#'
#' @param path a character path to two .bigWig files, or a data.table
#' with 2 columns, (forward, filepath) and reverse, only 1 row.
#' @inheritParams matchSeqStyle
#' @importFrom rtracklayer import.bw
#' @return a \code{\link{GRanges}} object of the file/s
#' @family utils
#'
readBigWig <- function(path, chrStyle = NULL) {
  if (is(path, "character")) {
    if (length(path) != 2) stop("readWig must have 2 wig files,
                              one forward strand and one reverse!")

    forwardPath <- grep("forward|fwd", path)
    reversePath <- grep("reverse|rev", path)
    if (length(forwardPath) == 1 & length(reversePath) == 1){
      forwardIndex <- forwardPath
      reverseIndex <- reversePath
    }

    forward <- import.bw(path[forwardIndex])
    reverse <- import.bw(path[reverseIndex])
  } else if (is(path, "data.table")) {
    if (!is.null(path$forward)) {
      forward <- import.bw(path$forward)
    } else forward <- import.wig(path$filepath)
    reverse <- import.bw(path$reverse)
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
  return(getGAlignments(fread(input = path, stringsAsFactors = TRUE)))
}

#' Load GRanges / GAlignments object from .ofst
#'
#' A much faster way to store, load and use bam files.\cr
#' .ofst is ORFik fast serialized object,
#' an optimized format for coverage reads with
#' cigar and replicate number. It uses the fst format as back-end:
#' \code{\link{fst-package}}.\cr A .ofst ribo seq file can compress the
#' information in a bam file from 5GB down to a few MB. This new files has
#' super fast reading time, only a few seconds, instead of minutes. It also has
#' random index access possibility of the file. \cr
#' .ofst is represented as a data.frane format with minimum 4 columns:\cr
#' 1. chromosome\cr  2. start (left most position) \cr 3. strand (+, -, *)\cr
#' 4. width (not added if cigar exists)\cr
#' 5. cigar (not needed if width exists):
#'  (cigar # M's, match/mismatch total) \cr
#' 5. score: duplicates of that read\cr
#' 6. size: qwidth according to reference of read\cr\cr
#' If file is from \code{\link{GAlignmentPairs}},
#' it will contain a cigar1, cigar2 instead
#' of cigar and start1 and start2 instead of start
#'
#' Other columns can be named whatever you want and added to meta columns.
#' Positions are 1-based, not 0-based as .bed.
#' Import with import.ofst
#' @param file a path to a .ofst file
#' @inheritParams readBam
#' @param seqinfo Seqinfo object, defaul NULL (created from ranges).
#' Add to avoid warnings later on differences in seqinfo.
#' @return a GAlignment, GAlignmentPairs or GRanges object,
#' dependent of if cigar/cigar1 is defined in .ofst file.
#' @importFrom fst read_fst
#' @export
#' @examples
#' ## GRanges
#' gr <- GRanges("1:1-3:-")
#' tmp <- file.path(tempdir(), "path.ofst")
#' # export.ofst(gr, file = tmp)
#' # import.ofst(tmp)
#' ## GAlignment
#' # Make input data.frame
#' df <- data.frame(seqnames = "1", cigar = "3M", start = 1L, strand = "+")
#' ga <- ORFik:::getGAlignments(df)
#' # export.ofst(ga, file = tmp)
#' # import.ofst(tmp)
import.ofst <- function(file, strandMode = 0, seqinfo = NULL) {
  df <- read_fst(file)
  if ("cigar" %in% colnames(df)) {
    getGAlignments(df, seqinfo = seqinfo)
  } else if ("cigar1" %in% colnames(df)) {
    getGAlignmentsPairs(df, strandMode = strandMode, seqinfo = seqinfo)
  } else getGRanges(df, seqinfo = seqinfo)
}

#' Load any type of sequencing reads
#'
#' Wraps around ORFik file format loaders and rtracklayer::import
#' and tries to speed up loading with the
#' use of data.table. Supports gzip, gz, bgz compression formats.
#' Also safer chromosome naming with the argument chrStyle
#'
#' NOTE: For wig/bigWig files you can send in 2 files, so that it automatically
#' merges forward and reverse stranded objects. You can also just send 1 wig/bigWig file,
#' it will then have "*" as strand.
#' @param path a character path to file (1 or 2 files),
#'  or data.table with 2 colums(forward&reverse)
#'  or a GRanges/Galignment/GAlignmentPairs object etc.
#'  If it is ranged object it will presume to be
#'  already loaded, so will return the object as it is,
#'  updating the seqlevelsStyle if given.
#' @inheritParams readBam
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
#' # Paired end bam strandMode 1:
#' fimport(bam_file, strandMode = 1)
#' # (will have no effect in this case, since it is not paired end)
#'
fimport <- function(path, chrStyle = NULL, param = NULL, strandMode = 0) {
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
          return(readBam(path, chrStyle, param, strandMode))
        } else if (all(fext %in% c("bigWig", "bw"))) {
          return(readBigWig(path, chrStyle))
        } else stop("only wig and valid bam format allowed for 2 files input!")
      } else if (length(path) == 1) { # Only 1 file path given
        if (fext == "bam") {
          if (pairedEndBam) path <- c(path, "paired-end")
          return(readBam(path, chrStyle, param, strandMode))
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
          return(matchSeqStyle(import.ofst(path, strandMode), chrStyle))
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
#' @examples
#' # Some fasta genome with existing fasta index in same folder
#' path <- system.file("extdata", "genome.fasta", package = "ORFik")
#' findFa(path)
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

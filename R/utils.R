
#' Converts different type of files to Granges
#'
#' column 5 will be set to score
#' Only Accepts bed files for now, standard format from Fantom5
#' @param x An data.frame from imported bed-file,
#'  to convert to GRanges
#' @param bed6 If bed6, no meta column is added
#' @return a GRanges object from bed
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
#' @return a GRanges object
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
      if (any(file_ext(filePath) == c("gzip", "gz", "bgz"))) {
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
#' @return a GAlignment object of bam file
#' @family utils
#'
readBam <- function(path, chrStyle = NULL) {
  if (is(path, "factor")) path <- as.character(path)
  return(matchSeqStyle(readGAlignments(path), chrStyle))
}

#' Custom wig reader
#'
#' Given 2 wig files, first is forward second is reverse.
#' Merge them and return as GRanges object.
#' @param path a character path to .bam file
#' @inheritParams matchSeqStyle
#' @importFrom rtracklayer import.wig
#' @return a GAlignment object of bam file
#' @family utils
#'
readWig <- function(path, chrStyle = NULL) {
  if (length(path) != 2) stop("readWig must have 2 wig files,
                              one forward strand and one reverse!")
  forward <- import.wig(path[1])
  reverse <- import.wig(path[2])
  strand(forward) <- "+"
  strand(reverse) <- "-"
  return(matchSeqStyle(c(forward, reverse), chrStyle))
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
#' @return a GAlignment/GRanges object depending on input.
fimport <- function(path, chrStyle = NULL) {
  if (is.character(path)) {
    if (all(file.exists(path))) {
      if (length(path) > 1) { # Multiple file paths
        if (all(file_ext(path) == "wig")) {
          return(readWig(path, chrStyle))
        } else stop("only wig format allowed for multiple files!")
      } else { # Only 1 file path given
        if (file_ext(path) == "bam") {
          return(readBam(path, chrStyle))
        } else if (file_ext(path) == "bed" |
                   file_ext(file_path_sans_ext(path,
                                               compression = TRUE)) == "bed") {
          return(fread.bed(path, chrStyle))
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
#' @param chrStyle a GRanges object, or a character style (Default: NULL) to
#' get seqlevelsStyle from. Is chromosome 1 called chr1 or 1, is mitocondrial
#' chromosome called MT or chrM etc. Will use 1st seqlevel-
#' style if more are present. Like: c("NCBI", "UCSC") -> pick "NCBI"
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
    } else if (is.gr_or_grl(chrStyle)) {
      seqlevelsStyle(range) <- seqlevelsStyle(chrStyle)[1]
    } else stop("chrStyle must be valid GRanges object, or a valid chr style!")
  }
  return(range)
}

#' Find optimized subset of valid reads
#'
#' If more than a million reads, keep only the ones that overlap within the
#' grl ranges.
#' @inheritParams validSeqlevels
#' @return the reads as GRanges or GAlignment
#' @family utils
#'
optimizeReads <- function(grl, reads) {
  seqMatch <- validSeqlevels(grl, reads)
  reads <- keepSeqlevels(reads, seqMatch, pruning.mode = "coarse")
  if (length(reads) > 1e6) { # speedup on big libraries
    reads <- reads[countOverlaps(reads, grl, type = "within") > 0]
    reads <- sort(reads)
  }
  return(reads)
}

#' Convert a GRanges Object to 1 width reads
#'
#' There are 4 ways of doing this
#' 1. Take 5' ends, reduce away rest (5prime)
#' 2. Take 3' ends, reduce away rest (3prime)
#' 3. Tile to 1-mers and include all (tileAll)
#' 4. Take middle point per GRanges (middle)
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
#'  A useful compression.
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

  if (addSizeColumn) {
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr),
                                      size = readWidths(gr, after.softclips))
  }
  if (addScoreColumn) {
    if (!is.null(cigar(gr))) {
      dt <- as.data.table(gr)[, .(.N, .I),
                              .(seqnames, start, end, strand, cigar)]
    } else {
      dt <- as.data.table(gr)[, .(.N, .I), .(seqnames, start, end, strand)]
    }
    dt <- dt[!duplicated(N), ] # this is right ?
    gr <- gr[dt$I]
    gr$score <- NULL
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr), score = as.integer(dt$N))
  }
  gr <- GRanges(gr)
  if (method == "5prime") {
    gr <- resize(gr, width = 1, fix = "start")
  } else if(method == "3prime") {
    gr <- resize(gr, width = 1, fix = "end")
  } else if(method == "tileAll") {
    gr <- unlist(tile(gr, width = 1), use.names = FALSE)
  } else if (method == "middle") {
    ranges(gr) <- IRanges(start(gr) + ceiling((end(gr) - start(gr)) / 2),
                          width = 1)
  } else stop("method not defined: must be 5prime, 3prime etc.")

  return(gr)
}

#' Convenience wrapper for Rsamtools FaFile
#' @param faFile FaFile, BSgenome or fasta/index file path used to find the
#'  transcripts
#' @importFrom Rsamtools FaFile
#' @importFrom methods is
#' @return a FaFile or BSgenome
#' @family utils
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
  }
  stop("faFile must be FaFile, BSgenome or valid filePath")
}

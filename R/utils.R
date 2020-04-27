
#' Converts bed style data.frame to Granges
#'
#' For info on columns, see:
#' https://www.ensembl.org/info/website/upload/bed.html
#' @param x A \code{\link{data.frame}} from imported bed-file,
#'  to convert to GRanges
#' @param skip.name default (TRUE), skip name column (column 4)
#' @return a \code{\link{GRanges}} object from bed
#' @family utils
#'
bedToGR <- function(x, skip.name = TRUE) {
  if (ncol(x) < 3) stop("bed file must have a least 3 columns!")
  starts <- x[, 2L] + 1L
  ends <- x[, 3L]
  gr <- GRanges(x[, 1L], IRanges(starts, ends))

  bed4 <- (ncol(x) >= 4) & !skip.name #name column
  if (bed4) mcols(gr) <- S4Vectors::DataFrame(mcols(gr), name = x[, 4L])
  bed5 <- ncol(x) >= 5 # score column
  if (bed5) mcols(gr) <- S4Vectors::DataFrame(mcols(gr), score = x[, 5L])
  bed6 <- ncol(x) >= 6
  if (bed6) strand(gr) <- x[, 6L]
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

#' Export as bed12 format
#'
#' bed format for multiple exons per group, as transcripts.
#' Can be use as alternative as a sparse .gff format for ORFs.
#' Can be direct input for ucsc browser or IGV
#'
#' If grl has no names, groups will be named 1,2,3,4..
#' @param grl A GRangesList
#' @param file a character path to valid output file name
#' @return NULL (File is saved as .bed)
#' @importFrom data.table fwrite
#' @export
#' @family utils
#' @examples
#' grl <- GRangesList(GRanges("1", c(1,3,5), "+"))
#' # export.bed12(grl, "output/path/orfs.bed")
export.bed12 <- function(grl, file){
  if (!is.grl(class(grl))) stop("grl, must be of class GRangesList")
  if (!is.character(file)) stop("file must be of class character")
  if (is.null(names(grl))) names(grl) <- seq.int(length(grl))
  grl <- sortPerGroup(grl, ignore.strand = TRUE) # <- sort bed way!

  dt.grl <- data.table(seqnames = seqnamesPerGroup(grl, FALSE))
  dt.grl$start <- as.integer(firstStartPerGroup(grl,keep.names = FALSE) -1)
  dt.grl$end <- lastExonEndPerGroup(grl, keep.names = FALSE)#non inclusive end
  dt.grl$name <- names(grl)
  dt.grl$score <- widthPerGroup(grl, keep.names = FALSE)
  dt.grl$strand <- strandPerGroup(grl, FALSE)
  dt.grl$thickStart <- dt.grl$start
  dt.grl$thickEnd <- dt.grl$end
  dt.grl$rgb <- rep(0, length(grl))
  dt.grl$blockCount <- numExonsPerGroup(grl)
  blockSizes <- paste(width(grl), collapse = ",")
  names(blockSizes) <- NULL
  dt.grl$blockSizes <- blockSizes
  relativeStarts <- (start(grl) -1) - dt.grl$start
  blockStarts <- paste(relativeStarts, collapse = ",")
  names(blockStarts) <- NULL
  dt.grl$blockStarts <- blockStarts

  #chromStart + chromStarts[last] + blockSizes[last])
  #must equal chromEnd.
  data.table::fwrite(x = dt.grl, file = file,
                     sep = "\t", col.names = FALSE, row.names = FALSE,
                     quote = FALSE)
  return(invisible(NULL))
}

#' Custom bam reader
#'
#' Safer version that handles the most important error done.
#' In the future will use a faster .bam loader for big .bam files in R.
#' @param path a character path to .bam file. If paired end bam files,
#' input must be a data.table with two columns: forward and reverse,
#' or a character vector of length 2, forward will be index 1.
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
  if (is(path, "data.table")) { # Will read pairs as combination
    message("ORFik reads paired end bam in as combination: c(forward, reverse)")
    return(matchSeqStyle(c(readGAlignments(path$forward),
                           readGAlignments(path$reverse)), chrStyle))
  }
  if (length(path) == 2){ # Will read pairs as combination
    message("ORFik reads paired end bam in as combination: c(forward, reverse)")
    return(matchSeqStyle(c(readGAlignments(path[1]),
                           readGAlignments(path[2])), chrStyle))
  }

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

#' Find pair of forward and reverse strand wig / bed files and
#' paired end bam files
#'
#' @param paths a character path at least one .wig / .bed file
#' @param f Default (c("forward", "fwd")
#' a character vector for forward direction regex.
#' @param r Default (c("reverse", "rev")
#' a character vector for reverse direction regex.
#' @param format default "wig", for bed do "bed". Also searches
#' compressions of these variants.
#' @return if not all are paired, return original list,
#' if they are all paired, return a data.table with matches as 2 columns
findNGSPairs <- function(paths, f = c("forward", "fwd"),
                         r = c("reverse", "rev"), format = "wig") {
  f <- paste0(f, "\\.", format, "*", collapse = "|")
  r <- paste0(r, "\\.", format, "*", collapse = "|")
  forwardPath <- grep(f, paths)
  reversePath <- grep(r, paths)

  if ((length(forwardPath) != length(reversePath)) |
      length(forwardPath) == 0 | length(reversePath) == 0) return(paths)
  dt <- data.table(forward = paths[forwardPath],
                   reverse = paths[reversePath], match = FALSE)
  dt.sub <- dt[, .(forward = gsub(pattern = f, x = forward, ""),
                   reverse = gsub(pattern = r, x = reverse, ""))]
  matches <- chmatch(dt.sub$forward, dt.sub$reverse)
  if (!anyNA(matches)) {
    dt <- dt[, .(forward, reverse = reverse[matches], match = TRUE)]
  }
  if (all(dt$match)) {
    return(dt)
  }
  return(paths)
}

#' Store GRanges object as .bedo
#'
#' .bedo is .bed ORFik, an optimized bed format for coverage reads with
#' read lengths .bedo is a text based format with columns (6 maximum):\cr
#' 1. chromosome\cr 2. start\cr 3. end\cr 4. strand\cr
#' 5. ref width (cigar # M's, match/mismatch total)\cr
#' 6. duplicates of that read\cr
#'
#' Positions are 1-based, not 0-based as .bed.
#' End will be removed if all ends equals all starts.
#' Import with import.bedo
#' @param object a GRanges object
#' @param out a character, location on disc (full path)
#' @return NULL, object saved to disc
#' @importFrom data.table fwrite
#'
export.bedo <- function(object, out) {
  if (!is(object, "GRanges")) stop("object must be GRanges")
  dt <- setDT(as.data.frame(object))
  uwidths <- unique(dt$width)
  if (all(uwidths == 1)) dt$end <- NULL
  dt$width <- NULL
  fwrite(dt, file = out)
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
#' @param path a character path to file (1 or 2 files),
#'  or data.table with 2 colums(forward&reverse)
#'  or a GRanges/Galignment object etc. If it is
#'  ranged object it will presume to be
#'  already loaded, so will return the object as it is.
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
      path <- path[path != ""]
    } else stop("When path is data.table,",
                "it must have 2 columns (forward&reverse)")
  }
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
#' @param method the method to reduce ranges, see info. (5prime defualt)
#' @param addScoreColumn logical (FALSE), if TRUE, add a score column that
#'  sums up the hits per unique range This will make each read unique, so
#'  that each read is 1 time, and score column gives the number of hits.
#'  A useful compression. If addSizeColumn is FALSE, it will not differentiate
#'  between reads with same start and stop, but different length. If
#'  addSizeColumn is FALSE, it will remove it.
#' @param addSizeColumn logical (FALSE), if TRUE, add a size column that
#'  for each read, that gives original width of read. Useful if you need
#'  original read lengths. This takes care of soft clips etc.
#' @inheritParams readWidths
#' @return  Converted GRanges object
#' @export
#' @family utils
#'
convertToOneBasedRanges <- function(gr, method = "5prime",
                                    addScoreColumn = FALSE,
                                    addSizeColumn = FALSE,
                                    after.softclips = TRUE,
                                    along.reference = FALSE) {
  if (is(gr, "GAlignmentPairs")) stop("Paired end reads not supported,
                                      load as GAlignments instead!")

  if (addSizeColumn & is.null(mcols(gr)$size)) {
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr),
                                      size = readWidths(gr, after.softclips,
                                                        along.reference))
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

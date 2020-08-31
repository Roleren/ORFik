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

#' Internal GRanges loader from fst data.frame
#' @param df a data.frame with columns minimum 4 columns:
#' seqnames, start, strand and width.\cr
#' Additional columns will be assigned as meta columns
#' @return GRanges object
getGRanges <- function(df) {

  ranges <- new2("IRanges", start = df$start,
                 width = df$width,
                 NAMES = df$NAMES,
                 elementMetadata = NULL,
                 check = FALSE)
  seqinfo <- Seqinfo(levels(df$seqnames))
  df$NAMES <- NULL
  if (ncol(df) == 4){
    mcols <- NULL
  } else {
    mcols <- df[,5:ncol(df)]
    if (ncol(df) == 5) {
      mcols <- data.frame(mcols)
      names(mcols) <- names(df)[5]
    }
  }
  mcols <- S4Vectors:::normarg_mcols(mcols, "GRanges", nrow(df))

  new2("GRanges", seqnames = Rle(df$seqnames), ranges = ranges, strand = Rle(df$strand),
       elementMetadata = mcols, seqinfo = seqinfo, check = FALSE)
}

#' Internal GAlignments loader from fst data.frame
#' @param df a data.frame with columns minimum 4 columns:
#' seqnames, start, strand and width.\cr
#' Additional columns will be assigned as meta columns
#' @return GAlignments object
getGAlignments <- function(df) {
  seqinfo <- Seqinfo(levels(df$seqnames))
  names <- df$NAMES
  df$NAMES <- NULL
  if (ncol(df) == 4){
    mcols <- NULL
  } else {
    mcols <- df[,5:ncol(df)]
    if (ncol(df) == 5) { # Hm... Is this safe ? What if a score is there ?
      mcols <- data.frame(mcols)
      names(mcols) <- names(df)[5]
    }
  }
  mcols <- S4Vectors:::normarg_mcols(mcols, "GRanges", nrow(df))
  new2("GAlignments", NAMES = names, seqnames = Rle(df$seqnames), start = df$start,
       cigar = as.character(df$cigar), strand = Rle(df$strand), elementMetadata = mcols,
       seqinfo = seqinfo, check = FALSE)

}

#' Internal GAlignmentPairs loader from fst data.frame
#' @param df a data.frame with columns minimum 6 columns:
#' seqnames, start1/start2 (integers), cigar1/cigar2 and strand\cr
#' Additional columns will be assigned as meta columns
#' @return GAlignmentPairs object
getGAlignmentsPairs <- function(df) {
  seqinfo <- Seqinfo(levels(df$seqnames))
  names <- df$NAMES
  df$NAMES <- NULL
  if (ncol(df) == 6){
    mcols <- NULL
  } else {
    mcols <- df[,7:ncol(df)]
    if (ncol(df) == 7) { # Hm... Is this safe ? What if a score is there ?
      mcols <- data.frame(mcols)
      names(mcols) <- names(df)[7]
    }
  }
  mcols <- S4Vectors:::normarg_mcols(mcols, "GRanges", nrow(df))
  # reverse strand for last
  strand2 <- strandTemp <- factor(df$strand, levels = c("+", "-", "*"))
  strandTemp[strand2 == "+"] <- "-"
  strandTemp[strand2 == "-"] <- "+"
  strand2 <- strandTemp
  new2("GAlignmentPairs",
       first = new2("GAlignments", NAMES = names, seqnames = Rle(df$seqnames), start = df$start1,
        cigar = as.character(df$cigar1), strand = Rle(df$strand), elementMetadata = mcols,
        seqinfo = seqinfo, check = FALSE),
       last = new2("GAlignments", NAMES = names, seqnames = Rle(df$seqnames), start = df$start2,
             cigar = as.character(df$cigar2), strand = Rle(strand2), elementMetadata = mcols,
             seqinfo = seqinfo, check = FALSE), check = FALSE)

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
#' @param rgb integer vector, default (0), either single integer or
#' vector of same size as grl to specify groups. It is adviced to not
#' use more than 8 different groups
#' @return NULL (File is saved as .bed)
#' @importFrom data.table fwrite
#' @export
#' @family utils
#' @examples
#' grl <- GRangesList(GRanges("1", c(1,3,5), "+"))
#' # export.bed12(grl, "output/path/orfs.bed")
export.bed12 <- function(grl, file, rgb = 0) {
  if (!is.grl(class(grl))) stop("grl, must be of class GRangesList")
  if (!is.character(file)) stop("file must be of class character")
  if (length(rgb) != 1 & length(rgb) != length(grl))
    stop("rgb must be integer of size 1 or length(grl)")
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
  dt.grl$rgb <- if(length(rgb) > 1) rgb else rep(rgb, length(grl))
  dt.grl$blockCount <- numExonsPerGroup(grl)
  blockSizes <- paste(width(grl), collapse = ",")
  names(blockSizes) <- NULL
  dt.grl$blockSizes <- blockSizes
  relativeStarts <- (start(grl) -1) - dt.grl$start
  blockStarts <- paste(relativeStarts, collapse = ",")
  names(blockStarts) <- NULL
  dt.grl$blockStarts <- blockStarts

  # Write without colnames
  data.table::fwrite(x = dt.grl, file = file,
                     sep = "\t", col.names = FALSE, row.names = FALSE,
                     quote = FALSE)
  return(invisible(NULL))
}

#' Convert data.frame to GAlignment object
#' @param x a data.frame / data.table with at least 4 columns:
#' seqnames, pos, cigar and strand, can also have a 5th column score.
#' @return a GAlignments object
makeGAlignmentsFromDataFrame <- function(x) {
  if (!all(c("seqnames", "start", "cigar", "strand") %in% colnames(x)))
    stop("x must at minimum have 4 columns named: seqnames, pos, cigar and strand")
  if (is.factor(x$cigar)) x$cigar <- as.character(x$cigar)

  ga <- GAlignments(seqnames = as.factor(x$seqnames), cigar = x$cigar,
                    pos = as.integer(x$start),
                    strand = factor(x$strand, levels = c("+", "-", "*")))
  if (!is.null(x$score)) {
    mcols(ga) <- DataFrame(mcols(ga), x$score)
  }
  return(ga)
}

#' Find pair of forward and reverse strand wig / bed files and
#' paired end bam files split in two
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

#' Store GAlignments object as .bedoc
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
#' Import with import.bedoc
#' @param object a GAlignments object
#' @param out a character, location on disc (full path)
#' @return NULL, object saved to disc
#' @importFrom data.table fwrite
#' @export
#'
export.bedoc <- function(object, out) {
  if (!is(object, "GAlignments"))
    stop("object must be of class GAlignments!")
  dt <- data.table(seqnames = as.character(seqnames(object)),
                     start = start(ranges(object)),
                     cigar = cigar(object),
                     strand = as.character(strand(object)))
  if (!is.null(mcols(object)$score)) dt$score = mcols(object)$score
  fwrite(dt, file = out)
}

#' Store GRanges / GAlignments object as .ofst
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
#' 6. size: qwidth according to reference of read\cr\cr
#' If file is from GAlignmentPairs, it will contain a cigar1, cigar2 instead
#' of cigar and start1 and start2 instead of start
#'
#' Other columns can be named whatever you want and added to meta columns.
#' Positions are 1-based, not 0-based as .bed.
#' Import with import.ofst
#' @param x a GRanges, GAlignments or GAlignmentPairs object
#' @param ... additional arguments for write_fst
#' @return NULL, object saved to disc
#' @importFrom fst write_fst
#' @export
#'
setGeneric("export.ofst", function(x,...) standardGeneric("export.ofst"))

#' @inherit export.ofst
#' @param file a character, location on disc (full path)
setMethod("export.ofst", "GRanges",
          function(x, file, ...) {
            df <- data.frame(seqnames = x@seqnames,
                             start = x@ranges@start,
                             width = x@ranges@width,
                             strand = x@strand)
            if (!is.null(x@ranges@NAMES)) df$NAMES <- x@ranges@NAMES
            if (ncol(x@elementMetadata) > 0) df <- cbind(df,x@elementMetadata)
            write_fst(df, file,...)
          })

#' @inherit export.ofst
#' @param file a character, location on disc (full path)
setMethod("export.ofst", "GAlignments",
          function(x, file, ...) {
            df <- data.frame(seqnames = x@seqnames,
                             start = x@start,
                             cigar = factor(x@cigar),
                             strand = x@strand)
            if (!is.null(x@NAMES)) df$NAMES <- x@ranges@NAMES
            if (ncol(x@elementMetadata) > 0) df <- cbind(df,x@elementMetadata)
            write_fst(df, file,...)
          })

#' @inherit export.ofst
#' @param file a character, location on disc (full path)
setMethod("export.ofst", "GAlignmentPairs",
          function(x, file, ...) {
            # There is always equal seqname in a pair,
            # strand is always reverse of the other
            df <- data.frame(seqnames = x@seqnames,
                             start1 = x@first@start,
                             start2 = x@last@start,
                             cigar1 = factor(x@first@cigar),
                             cigar2 = factor(x@last@cigar),
                             strand = x@strand)
            if (!is.null(x@NAMES)) df$NAMES <- x@NAMES # Check that this is correct
            if (ncol(x@elementMetadata) > 0) df <- cbind(df, x@elementMetadata)
            write_fst(df, file,...)
          })

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
#' @return the reads as GRanges,  GAlignment or GAlignmentPairs
#' @importFrom GenomicAlignments first
#' @family utils
#'
optimizeReads <- function(grl, reads) {
  seqMatch <- validSeqlevels(grl, reads)
  reads <- keepSeqlevels(reads, seqMatch, pruning.mode = "coarse")

  reads <- reads[countOverlaps(reads, grl, type = "within") > 0]
  reads <- if (is(reads, "GAlignmentPairs")) {
    reads <- reads[order(GenomicAlignments::first(reads))]
    } else reads <- sort(reads)

  if (length(reads) == 0) warning("No reads left in 'reads' after",
                                    "optimisation!")

  return(reads)
}

#' Convert a GRanges Object to 1 width reads
#'
#' There are 5 ways of doing this\cr
#' 1. Take 5' ends, reduce away rest (5prime)\cr
#' 2. Take 3' ends, reduce away rest (3prime)\cr
#' 3. Tile to 1-mers and include all (tileAll)\cr
#' 4. Take middle point per GRanges (middle)\cr
#' 5. Get original with metacolumns (None)\cr
#' You can also do multiple at a time, then output is GRangesList, where
#' each list group is the operation (5prime is [1], 3prime is [2] etc)\cr
#' Many other ways to do this have their own functions, like startSites and
#' stopSites etc.
#' To retain information on original width, set addSizeColumn to TRUE.
#' To compress data, 1 GRanges object per unique read, set addScoreColumn to
#' TRUE. This will give you a score column with how many duplicated reads there
#' were in the specified region.
#'
#' NOTE: For special case of GAlignmentPairs, 5prime will only use left (first)
#' 5' end and read and 3prime will use only right (last) 3' end of read
#' in pair. tileAll and middle can possibly find poinst that are not in the
#' reads since: lets say pair is 1-5 and 10-15, middle is 7, which is not in
#' the read.
#'
#' @param gr GRanges, GAlignment or GAlignmentPairs object to reduce.
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
#' @importFrom GenomicAlignments first
#' @importFrom GenomicAlignments last
#' @return  Converted GRanges object
#' @export
#' @family utils
#' @examples
#' gr <- GRanges("chr1", 1:10,"+")
#' # 5 prime ends
#' convertToOneBasedRanges(gr)
#' # is equal to convertToOneBasedRanges(gr, method = "5prime")
#' # 3 prime ends
#' convertToOneBasedRanges(gr, method = "3prime")
#' # With lengths
#' convertToOneBasedRanges(gr, addSizeColumn = TRUE)
#' # With score (# of replicates)
#' gr <- rep(gr, 2)
#' convertToOneBasedRanges(gr, addSizeColumn = TRUE, addScoreColumn = TRUE)
#'
convertToOneBasedRanges <- function(gr, method = "5prime",
                                    addScoreColumn = FALSE,
                                    addSizeColumn = FALSE,
                                    after.softclips = TRUE,
                                    along.reference = FALSE) {

  if (addSizeColumn & is.null(mcols(gr)$size)) {
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr),
                                      size = readWidths(gr, after.softclips,
                                                        along.reference))
  }
  if (addScoreColumn) {
    dt <- data.table(seqnames = as.character(seqnames(gr)),
                     start = start(ranges(gr)),
                     end = end(ranges(gr)),
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

#' Collapse duplicated reads
#'
#' For every GAlignments read, with the same:
#' seqname, start, cigar and strand, collapse and give a new
#' meta column called "score", which contains the number of duplicates
#' of that read.
#' @param x a GAlignments object
#' @param ... alternative arguments. addScoreColumn = TRUE, if FALSE,
#' only collapse and not add score column.
#' @return a GAlignments object
#' @export
setGeneric("collapseDuplicatedReads", function(x,...) standardGeneric("collapseDuplicatedReads"))

#' @inherit collapseDuplicatedReads
#' @param addScoreColumn = TRUE, if FALSE,
#' only collapse and not add score column.
setMethod("collapseDuplicatedReads", "GAlignments",
          function(x, addScoreColumn = TRUE) {

            dt <- data.table(seqnames = as.character(seqnames(x)),
                             start = start(ranges(x)),
                             cigar = cigar(x),
                             strand = as.character(strand(x)))
            dt <- dt[, .(score = .N), .(seqnames, start, cigar, strand)]
            if (!addScoreColumn) dt$score <- NULL
            return(makeGAlignmentsFromDataFrame(dt))
          })

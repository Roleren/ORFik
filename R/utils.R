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
  if (!all(c("seqnames", "start", "width", "strand") %in% colnames(df)))
    stop("df must at minimum have 4 columns named: seqnames, start, width and strand")
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
#' seqnames, start ("pos" in final GA object), strand and width.\cr
#' Additional columns will be assigned as meta columns
#' @return GAlignments object
getGAlignments <- function(df) {
  if (!all(c("seqnames", "start", "cigar", "strand") %in% colnames(df)))
    stop("df must at minimum have 4 columns named: seqnames, start, cigar and strand")
  seqinfo <- Seqinfo(levels(df$seqnames))
  names <- df$NAMES
  if (!is.null(df$NAMES)) df$NAMES <- NULL
  if (ncol(df) == 4){
    mcols <- NULL
  } else {
    mcols <- df[,5:ncol(df)]
    if (ncol(df) == 5) { # Hm... Is this safe ? What if a score is there ?
      mcols <- data.frame(mcols)
      names(mcols) <- names(df)[5]
    }
  }
  df$strand <- factor(df$strand, levels = c("+", "-", "*"))
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
  if (!is.null(df$NAMES)) df$NAMES <- NULL
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
  strand2 <- strandTemp <- df$strand <-
    factor(df$strand, levels = c("+", "-", "*"))
  strandTemp[strand2 == "+"] <- "-"
  strandTemp[strand2 == "-"] <- "+"
  strand2 <- strandTemp
  new2("GAlignmentPairs",
       first = new2("GAlignments", NAMES = names, seqnames = Rle(df$seqnames), start = df$start1,
        cigar = as.character(df$cigar1), strand = Rle(df$strand),
        seqinfo = seqinfo, check = FALSE,
        elementMetadata = DataFrame(data.frame(matrix(nrow = nrow(df), ncol = 0)))),
       last = new2("GAlignments", NAMES = names, seqnames = Rle(df$seqnames), start = df$start2,
             cigar = as.character(df$cigar2), strand = Rle(factor(strand2, levels)),
             seqinfo = seqinfo, check = FALSE,
             elementMetadata = DataFrame(data.frame(matrix(nrow = nrow(df), ncol = 0)))),
       isProperPair = rep(TRUE, nrow(df)),
       elementMetadata = mcols, check = FALSE)
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
#'  sums up the hits per unique range. This will make each read unique, so
#'  that each read is 1 time, and score column gives the number of
#'  collapsed hits.
#'  A useful compression. If addSizeColumn is FALSE, it will not differentiate
#'  between reads with same start and stop, but different length. If
#'  addSizeColumn is FALSE, it will remove it. Collapses after conversion.
#' @param addSizeColumn logical (FALSE), if TRUE, add a size column that
#'  for each read, that gives original width of read. Useful if you need
#'  original read lengths. This takes care of soft clips etc.
#'  If collapsing reads, each unique range will be grouped also by size.
#' @param reuse.score.column logical (TRUE), if addScoreColumn is TRUE,
#'  and a score column exists, will sum up the scores to create a new score.
#'  If FALSE, will skip old score column and create new according to number
#'  of replicated reads after conversion.
#'  If addScoreColumn is FALSE, this argument is ignored.
#' @inheritParams readWidths
#' @importFrom GenomicAlignments first
#' @importFrom GenomicAlignments last
#' @return Converted GRanges object
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
                                    along.reference = FALSE,
                                    reuse.score.column = TRUE) {
  if (addSizeColumn & is.null(mcols(gr)$size)) {
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr),
                                      size = readWidths(gr, after.softclips,
                                                        along.reference))
  }
  # Convert to positions wanted
  if (!is(gr, "GRanges")) gr <- GRanges(gr)
  if (method == "5prime") {
    gr <- resize(gr, width = 1, fix = "start")
  } else if(method == "3prime") {
    gr <- resize(gr, width = 1, fix = "end")
  } else if(method %in% c("None", "none")) {
  } else if(method == "tileAll") {
    gr <- unlist(tile(gr, width = 1), use.names = FALSE)
  } else if (method == "middle") {
    ranges(gr) <- IRanges(start(gr) + ceiling((end(gr) - start(gr)) / 2),
                          width = 1)
  } else stop("invalid type: must be 5prime, 3prime, None, tileAll or middle")
  # Collapse after conversion
  if (addScoreColumn) {
    gr <- collapseDuplicatedReads(gr, addSizeColumn = addSizeColumn,
                            reuse.score.column = reuse.score.column)
  }
  return(gr)
}

#' Merge reads by sum of existing scores
#'
#' If you have multiple reads a same location but different read lengths,
#' specified in meta column "size", it will sum up the scores
#' (number of replicates) for all reads at that position
#' @param x a GRanges object
#' @return merged GRanges object
#' @examples
#' gr_s1 <- rep(GRanges("chr1", 1:10,"+"), 2)
#' gr_s2 <- GRanges("chr1", 1:12,"+")
#' gr2 <- GRanges("chr1", 21:40,"+")
#' gr <- c(gr_s1, gr_s2, gr2)
#' res <- convertToOneBasedRanges(gr,
#'    addScoreColumn = TRUE, addSizeColumn = TRUE)
#' ORFik:::collapse.by.scores(res)
#'
collapse.by.scores <- function(x) {
  dt <- data.table(seqnames = as.character(seqnames(x)),
                   start = start(ranges(x)),
                   end = end(ranges(x)),
                   strand = as.character(strand(x)),
                   score = mcols(x)$score)
  dt <- dt[, .(score = sum(score)), .(seqnames, start, end, strand)]
  # TODO change makeGRangesFromDataFrame to internal fast function
  return(makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE))
}

#' Collapse duplicated reads
#'
#' For every GRanges, GAlignments read, with the same:
#' seqname, start, (cigar) / width and strand, collapse and give a new
#' meta column called "score", which contains the number of duplicates
#' of that read. If score column already exists, will return input object!
#' @param x a GRanges, GAlignments or GAlignmentPairs object
#' @param ... alternative arguments. addScoreColumn = TRUE, if FALSE,
#' only collapse and not add score column.
#' @return a GRanges, GAlignments or GAlignmentPairs object, same as input
#' @export
#' @examples
#' gr <- rep(GRanges("chr1", 1:10,"+"), 2)
#' collapseDuplicatedReads(gr)
setGeneric("collapseDuplicatedReads", function(x,...) standardGeneric("collapseDuplicatedReads"))

#' @inherit collapseDuplicatedReads
#' @param addScoreColumn = TRUE, if FALSE,
#' only collapse and not keep score column.
#' @inheritParams convertToOneBasedRanges
setMethod("collapseDuplicatedReads", "GRanges",
          function(x, addScoreColumn = TRUE, addSizeColumn = FALSE,
                   reuse.score.column = TRUE) {
    if (addSizeColumn) {
      if (!("size" %in% colnames(mcols(x))))
        stop("addSizeColumn is TRUE, and no size column found!")
    }

    dt <- data.table(seqnames = as.character(seqnames(x)),
                     start = start(ranges(x)),
                     end = end(ranges(x)),
                     strand = as.character(strand(x)))

    if (reuse.score.column & ("score" %in% colnames(mcols(x)))) { # reuse
      dt[, score := mcols(x)$score]
      if (addSizeColumn) {
        dt[, size := mcols(x)$size]
        dt <- dt[, .(score = sum(score)), .(seqnames, start, end, strand, size)]
      } else {
        dt <- dt[, .(score = sum(score)), .(seqnames, start, end, strand)]
      }
    } else { # Do not reuse or "score" does not exist
      if (addSizeColumn) {
        dt[, size := mcols(x)$size]
        dt <- dt[, .(score = .N), .(seqnames, start, end, strand, size)]
      } else {
        dt <- dt[, .(score = .N), .(seqnames, start, end, strand)]
      }
    }
    if (!addScoreColumn) dt$score <- NULL
    # TODO change makeGRangesFromDataFrame to internal fast function
    return(makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE))
})

#' @inherit collapseDuplicatedReads
#' @param addScoreColumn = TRUE, if FALSE,
#' only collapse and not add score column.
setMethod("collapseDuplicatedReads", "GAlignments",
          function(x, addScoreColumn = TRUE) {
            if ("score" %in% colnames(mcols(x))) return(x)

            dt <- data.table(seqnames = factor(seqnames(x)),
                             start = start(ranges(x)),
                             cigar = cigar(x),
                             strand = factor(strand(x)))
            dt <- dt[, .(score = .N), .(seqnames, start, cigar, strand)]
            if (!addScoreColumn) dt$score <- NULL
            return(getGAlignments(dt))
          })

#' @inherit collapseDuplicatedReads
#' @param addScoreColumn = TRUE, if FALSE,
#' only collapse and not add score column.
setMethod("collapseDuplicatedReads", "GAlignmentPairs",
          function(x, addScoreColumn = TRUE) {
            if ("score" %in% colnames(mcols(x))) return(x)

            dt <- data.table(seqnames = factor(x@first@seqnames),
                             start1 = x@first@start,
                             start2 = x@last@start,
                             cigar1 = factor(x@first@cigar),
                             cigar2 = factor(x@last@cigar),
                             strand = factor(x@first@strand))
            dt <- dt[, .(score = .N), .(seqnames, start1, start2,
                                        cigar1, cigar2, strand)]
            if (!addScoreColumn) dt$score <- NULL
            return(getGAlignmentsPairs(dt))
          })

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

#' Export as wiggle format
#'
#' Will create 2 files, 1 for + strand (*_forward.wig)
#' and 1 for - strand (*_reverse.wig). If all
#' ranges are * stranded, will output 1 file.
#' Can be direct input for ucsc browser or IGV
#'
#' @references https://genome.ucsc.edu/goldenPath/help/wiggle.html
#' @param x A GRangesList, GAlignment GAlignmentPairs with score column.
#' Will be converted to 5' end position of original range. If score column
#' does not exist, will group ranges and give replicates as score column.
#' @param file a character path to valid output file name
#' @return invisible(NULL) (File is saved as 2 .wig files)
#' @importFrom rtracklayer export.wig
#' @export
#' @family utils
#' @examples
#' x <- c(GRanges("1", c(1,3,5), "-"), GRanges("1", c(1,3,5), "+"))
#' # export.wiggle(x, "output/path/rna.wig")
export.wiggle <- function(x, file) {
  if (!(is(x, "GRanges") | is(x, "GAlignmentPairs") | is(x, "GAlignments")))
     stop("x must be GRanges, GAlignments or GAlignmentPairs")
  if (!is(x, "GRanges")) x <- GRanges(x)

  x <- resize(x, width = 1, fix = "start")
  if (!("score" %in% colnames(mcols(x)))) {
    x <- convertToOneBasedRanges(x, method = "None",
                                 addScoreColumn = TRUE,
                                 addSizeColumn = FALSE)
  } else { # merge reads by sum of existing scores
    x <- collapse.by.scores(x)
  }
  strands <- as.character(strand(x))
  if (all(strands == "*")) {
    file <- gsub("\\.wig", "", file)
    file <- paste0(file, ".wig")
    export.wig(x, file)
  } else {
    file <- gsub("\\.wig", "", file)
    forward_file <- paste0(file, "_forward.wig")
    reverse_file <- paste0(file, "_reverse.wig")
    export.wig(x[strandBool(x)], forward_file)
    export.wig(x[!strandBool(x)], reverse_file)
  }

  return(invisible(NULL))
}

#' Export as bigWig format
#'
#' Will create 2 files, 1 for + strand (*_forward.bigWig)
#' and 1 for - strand (*_reverse.bigWig). If all
#' ranges are * stranded, will output 1 file.
#' Can be direct input for ucsc browser or IGV
#'
#' @references https://genome.ucsc.edu/goldenPath/help/bigWig.html
#' @param x A GRangesList, GAlignment GAlignmentPairs with score column.
#' Will be converted to 5' end position of original range. If score column
#' does not exist, will group ranges and give replicates as score column.
#' @param file a character path to valid output file name
#' @return invisible(NULL) (File is saved as 2 .bigWig files)
#' @importFrom rtracklayer export.bw
#' @export
#' @family utils
#' @examples
#' x <- c(GRanges("1", c(1,3,5), "-"), GRanges("1", c(1,3,5), "+"))
#' # export.bigWig(x, "output/path/rna.bigWig")
export.bigWig <- function(x, file) {
  if (!(is(x, "GRanges") | is(x, "GAlignmentPairs") | is(x, "GAlignments")))
    stop("x must be GRanges, GAlignments or GAlignmentPairs")
  if (!is(x, "GRanges")) x <- GRanges(x)

  x <- resize(x, width = 1, fix = "start")
  if (!("score" %in% colnames(mcols(x)))) {
    x <- convertToOneBasedRanges(x, method = "None",
                                 addScoreColumn = TRUE,
                                 addSizeColumn = FALSE)
  } else { # merge reads by sum of existing scores
    x <- collapse.by.scores(x)
  }
  strands <- as.character(strand(x))
  if (all(strands == "*")) {
    file <- gsub("\\.bigWig", "", file)
    file <- paste0(file, ".bigWig")
    export.wig(x, file)
  } else {
    file <- gsub("\\.bigWig", "", file)
    forward_file <- paste0(file, "_forward.bigWig")
    reverse_file <- paste0(file, "_reverse.bigWig")
    export.bw(x[strandBool(x)], forward_file)
    export.bw(x[!strandBool(x)], reverse_file)
  }
  return(invisible(NULL))
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
#' A fast way to store, load and use bam files.
#' (we now recommend using \code{link{export.ofst}} instead!)\cr
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
#' @examples
#' ## GRanges
#' gr <- GRanges("1:1-3:-")
#' # export.ofst(gr, file = "path.ofst")
#' ## GAlignment
#' # Make input data.frame
#' df <- data.frame(seqnames = "1", cigar = "3M", start = 1L, strand = "+")
#' ga <- ORFik:::getGAlignments(df)
#' # export.ofst(ga, file = "path.ofst")
setGeneric("export.ofst", function(x, ...) standardGeneric("export.ofst"))

#' @inherit export.ofst
#' @param file a character, location on disc (full path)
setMethod("export.ofst", "GRanges",
          function(x, file, ...) {
            df <- data.frame(seqnames = x@seqnames,
                             start = x@ranges@start,
                             width = x@ranges@width,
                             strand = x@strand)
            if (!is.null(x@ranges@NAMES)) df$NAMES <- x@ranges@NAMES
            if (ncol(x@elementMetadata) > 0)
              df <- as.data.frame(cbind(df, x@elementMetadata))
            write_fst(df, file,...)
          })

#' @inherit export.ofst
#' @param file a character, location on disc (full path)
setMethod("export.ofst", "GAlignments",
          function(x, file, ...) {
            df <- data.frame(seqnames = x@seqnames,
                             start = x@start,
                             cigar = factor(x@cigar, levels = unique(x@cigar)),
                             strand = x@strand)
            if (!is.null(x@NAMES)) df$NAMES <- x@ranges@NAMES
            if (ncol(x@elementMetadata) > 0)
              df <- as.data.frame(cbind(df, x@elementMetadata))
            write_fst(df, file,...)
          })

#' @inherit export.ofst
#' @param file a character, location on disc (full path)
setMethod("export.ofst", "GAlignmentPairs",
          function(x, file, ...) {
            # There is always equal seqname in a pair,
            # strand is always reverse of the other
            df <- data.frame(seqnames = x@first@seqnames,
                             start1 = x@first@start,
                             start2 = x@last@start,
                             cigar1 = factor(x@first@cigar,
                                             levels = unique(x@first@cigar)),
                             cigar2 = factor(x@last@cigar,
                                             levels = unique(x@last@cigar)),
                             strand = x@first@strand)
            if (!is.null(x@NAMES)) df$NAMES <- x@NAMES # Check that this is correct
            if (ncol(x@elementMetadata) > 0)
              df <- as.data.frame(cbind(df, x@elementMetadata))
            write_fst(df, file,...)
          })

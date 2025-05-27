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
#' @references
#' \url{https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bed-format}
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
#'  Will be converted to 5' end position of original range. If score column
#'  does not exist, will group ranges and give replicates as score column.
#'  Since bigWig needs a score column to represent counts!
#' @param file a character path to valid output file name
#' @param is_pre_collapsed logical, default FALSE. Have you already
#'  collapsed reads with collapse.by.scores,
#'  so each positions is only in 1 GRanges object with
#'  a score column per readlength?
#'  Set to TRUE, only if you are sure, will give a speedup.
#' @param split.by.strand logical, default TRUE. Split bigWig into 2 files,
#'  one for each strand.
#' @param seq_info a Seqinfo object, default seqinfo(x).
#'  Must have non NA seqlengths defined!
#' @return invisible(NULL) (File is saved as 2 .bigWig files)
#' @importFrom rtracklayer export.bw
#' @export
#' @family utils
#' @examples
#' x <- c(GRanges("1", c(1,3,5), "-"), GRanges("1", c(1,3,5), "+"))
#' seqlengths(x) <- 10
#' file <- file.path(tempdir(), "rna.bigWig")
#' # export.bigWig(x, file)
#' # export.bigWig(covRleFromGR(x), file)
export.bigWig <- function(x, file, split.by.strand = TRUE,
                          is_pre_collapsed = FALSE, seq_info = seqinfo(x)) {
  if(anyNA(seqlengths(seq_info))) stop("seqinfo of x must be defined and have defined seqlengths!")
  if (!(is(x, "GRanges") | is(x, "GAlignmentPairs") | is(x, "GAlignments") | is(x, "covRle")))
    stop("x must be GRanges, GAlignments, GAlignmentPairs or covRLE")

  if (!is(x, "covRle")) {
    if (!is(x, "GRanges")) { x <- GRanges(x, seqinfo = seq_info)
    } else {seqlevels(x) <- seqlevels(seq_info); seqinfo(x) <- seq_info}

    if (!all(width(x) == 1)) x <- resize(x, width = 1, fix = "start")
    if (!("score" %in% colnames(mcols(x)))) {
      x <- convertToOneBasedRanges(x, method = "None",
                                   addScoreColumn = TRUE,
                                   addSizeColumn = FALSE)
    } else { # merge reads by sum of existing scores
      if (!is_pre_collapsed) x <- collapse.by.scores(x)
    }
    strands <- as.character(strand(x))
    can_split_by_strand <- !all(strands == "*")
    make_single_file <- !can_split_by_strand | !split.by.strand
    x <- covRleFromGR(x, ignore.strand = make_single_file)
  } else make_single_file <- !strandMode(x) | !split.by.strand


  if (make_single_file) {
    file <- gsub("\\.bigWig", "", file, ignore.case = TRUE)
    file <- paste0(file, ".bigWig")
    if (length(r(x)) > 0) {
      x@forward <- f(x) + r(x)
    }
    export.bw(GRanges(f(x)), file)
  } else {
    file <- gsub("\\.bigWig", "", file, ignore.case = TRUE)
    forward_file <- paste0(file, "_forward.bigWig")
    reverse_file <- paste0(file, "_reverse.bigWig")
    export.bw(GRanges(f(x)), forward_file)
    export.bw(GRanges(r(x)), reverse_file)
  }
  return(invisible(NULL))
}


#' Export as fstwig (fastwig) format
#'
#' Will create 2 files, 1 for + strand (*_forward.fstwig)
#' and 1 for - strand (*_reverse.fstwig). If all
#' ranges are * stranded, will output 1 file.
#'
#' @references "TODO"
#' @param x A GRangesList, GAlignment GAlignmentPairs with score column
#'  or coverage RLElist
#' Will be converted to 5' end position of original range. If score column
#' does not exist, will group ranges and give replicates as score column.
#' @inheritParams fst::write_fst
#' @param file a character path to valid output file name
#' @param by.readlength logical, default TRUE
#' @param by.chromosome logical, default TRUE
#' @return invisible(NULL) (File is saved as 2 .fstwig files)
#' @export
#' @family utils
#' @examples
#' x <- c(GRanges("1", c(1,3,5), "-"), GRanges("1", c(1,3,5), "+"))
#' x$size <- rep(c(28, 29), length.out = length(x))
#' x$score <- c(5,1,2,5,1,6)
#' seqlengths(x) <- 5
#' # export.fstwig(x, "~/Desktop/ribo")
export.fstwig <- function(x, file, by.readlength = TRUE,
                          by.chromosome = TRUE, compress = 50) {
  # TODO: Implement split.by.strand argument!
  if (length(x) == 0) stop("length of x is 0, no values to export!")
  if (anyNA(seqlengths(x))) stop("Define seqlength to make sure tails are correct")

  if ((is(x, "GRanges") | is(x, "GAlignmentPairs") | is(x, "GAlignments"))) {
    if (!is(x, "GRanges")) x <- GRanges(x)
    message("-- Creating covRle object from GRanges")
    strands <- as.character(strand(x))
    strandBool <- strandBool(x)
    strand_mode <- ifelse(all(strands == "*"), "single", "double")
    # Get coverage
    if (by.readlength) {
      readlengths <- readWidths(x)
    } else readlengths <- rep.int(1, length(x))
    unique_lengths <- sort(unique(readlengths))
    x_b <- x_p <- x_m <- list()
    message("- Read length: ", appendLF = FALSE)
    for (length in unique_lengths) {
      message(", ", length, appendLF = FALSE)
      if (strand_mode == "single") { # b = both strands
        x_b <- c(x_b, coverage(x[readlengths == length], weight = "score"))
      } else { # p = pluss strand, m = minus strand
        hits <- which(readlengths == length)
        x_p <- c(x_p, coverage(x[hits][strandBool[hits]], weight = "score"))
        x_m <- c(x_m, coverage(x[hits][!strandBool[hits]], weight = "score"))
      }
    }

    if (strand_mode != "single") { # Add check for readlength
      names(x_p) <- names(x_m) <- unique_lengths
    } else {
      names(x_b) <-  unique_lengths
    }

  } else if (is(x, "covRle")) {
    if (strandMode(x)) {
      x_p <- f(x)
      x_m <- r(x)
      strand_mode <- "double"
    }
  } else if (is(x, "RleList")) {
    x_b <- x
    strand_mode <- "single"
  } else stop("x must be GRanges, GAlignments, GAlignmentPairs, RleList or covRle")
  chrs <- seqlevels(x)
  if (!dir.exists(dirname(file))) dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
  if (by.readlength) {
    message("- Output fstwig for chr:")
    for (chr in chrs) {
      message(chr)
      if (strand_mode == "single") {
        stop("Not implemented")
      } else {
        if (by.readlength) {
          lengths_to_use <- as.character(sort(unique_lengths))
          dtt_p <- data.table()
          dtt_m <- data.table()
          message("- Inserting column for read length:", appendLF = FALSE)
          for(length in lengths_to_use) {
            message(", ", length, appendLF = FALSE)
            dtt_p <-cbind(dtt_p, data.table(unlist(IntegerList(x_p[[length]][chr]))))
            dtt_m <-cbind(dtt_m, data.table(unlist(IntegerList(x_m[[length]][chr]))))
          }
          colnames(dtt_p) <- colnames(dtt_m) <- lengths_to_use
          dt <- list(dtt_p, dtt_m)
        } else {
          dt <- list(data.table(unlist(IntegerList(x_p[chr]))), data.table(unlist(IntegerList(x_m[chr]))))
        }

      }
      file_chr <- paste0(file, "_", chr)
      save.fstwig(dt, file_chr, compress = compress)
      rm(dt)
    }
  } else {
    stop("Not implemented")
  }

  return(invisible(NULL))
}

save.fstwig <- function(x, file, compress = 50) {
  mode <- ifelse(is(x, "data.table"), "single", "double")
  file <- gsub("\\.fstwig", "", file, ignore.case = TRUE)
  if (mode == "single") {
    file <- paste0(file, ".fstwig")
    fst::write_fst(x, file, compress = compress)
  } else {
    file <- gsub("\\.fstwig", "", file, ignore.case = TRUE)
    forward_file <- paste0(file, "_forward.fstwig")
    reverse_file <- paste0(file, "_reverse.fstwig")

    fst::write_fst(x[[1]], path = forward_file, compress = compress)
    fst::write_fst(x[[2]], path = reverse_file, compress = compress)
  }
  return(invisible(NULL))
}

export.cov <- function(x, file, seqinfo, split.by.strand = TRUE,
                       weight = "score", format = "qs") {
  stopifnot(is(seqinfo, "Seqinfo"))
  stopifnot(format %in% c("qs", "rds"))
  if (!is(x, "covRle") & !is(x, "RleList")) {
    seqlevels(x) <- seqlevels(seqinfo)
    seqinfo(x) <- seqinfo
    x <- covRleFromGR(x, weight = weight, ignore.strand = !split.by.strand)
  }
  seqinfo(x) <- seqinfo

  format <- paste0("cov", format)
  file <- paste0(gsub(paste0("\\.", format, "$"), "", file, ignore.case = TRUE),
                 ".", format)
  save_RDSQS(x, file = file)
}

export.covlist <- function(x, file, seqinfo, split.by.strand = TRUE,
                       weight = "score", verbose = TRUE, format = "qs") {

  if (!is(x, "covRleList")) {
    stopifnot(is(seqinfo, "Seqinfo"))
    seqlevels(x) <- seqlevels(seqinfo)
    seqinfo(x) <- seqinfo
    x <- covRleListFromGR(x, weight = weight, ignore.strand = !split.by.strand,
                          verbose = verbose)
  }

  format <- paste0("cov", format)
  file <- paste0(gsub(paste0("\\.", format, "$"), "", file, ignore.case = TRUE),
                 ".", format)
  save_RDSQS(x, file = file)
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
#' @param file a character, location on disc (full path)
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
setGeneric("export.ofst", function(x, file, ...) standardGeneric("export.ofst"))

#' @inherit export.ofst
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

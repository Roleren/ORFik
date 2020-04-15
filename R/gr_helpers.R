#' Group GRanges
#'
#' It will group / split the GRanges object by the argument `other`.
#' For example if you would like to to group GRanges object by gene,
#' set other to gene names.
#'
#' If `other` is not specified function will try to use the names of the
#' GRanges object. It will then be similar to `split(gr, names(gr))`.
#'
#' It is important that all groups in `other` are unique, otherwise
#' duplicates will be grouped together.
#' @param gr a GRanges object
#' @param other a vector of unique names to group by (default: NULL)
#' @return a GRangesList named after names(Granges) if other is NULL, else
#' names are from unique(other)
#' @export
#' @importFrom S4Vectors nrun
#' @examples
#' ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
#'                      ranges = IRanges(start = c(1, 10, 20),
#'                                       end = c(5, 15, 25)),
#'                      strand = "+")
#' ORFranges2 <- GRanges("1",
#'                       ranges = IRanges(start = c(20, 30, 40),
#'                                        end = c(25, 35, 45)),
#'                       strand = "+")
#' names(ORFranges) = rep("tx1_1", 3)
#' names(ORFranges2) = rep("tx1_2", 3)
#' grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2)
#' gr <- unlist(grl, use.names = FALSE)
#' ## now recreate the grl
#' ## group by orf
#' grltest <- groupGRangesBy(gr) # using the names to group
#' identical(grl, grltest) ## they are identical
#'
#' ## group by transcript
#' names(gr) <- txNames(gr)
#' grltest <- groupGRangesBy(gr)
#' identical(grl, grltest) ## they are not identical
#'
groupGRangesBy <- function(gr, other = NULL) {
  if (!is(gr, "GRanges")) stop("gr must be GRanges Object")
  if (is.null(other)) { # if not using other
    if (is.null(names(gr))) stop("gr object have no names")
    l <- S4Vectors::Rle(names(gr))
  } else { # else use other
    if (length(gr) != length(other))
      stop(" in GroupGRangesByOther: lengths of gr and other does not match")
    l <- S4Vectors::Rle(other)
  }
  grouping <- rep.int(seq.int(nrun(l)), runLength(l))
  grl <- split(gr, grouping)
  if (is.null(other)) {
    names(grl) <- unique(names(gr))
  } else {
    names(grl) <- unique(other)
  }
  return(grl)
}


#' Get read widths
#'
#' Input any reads, e.g. ribo-seq object and get width of reads, this is to
#' avoid confusion between width, qwidth and meta column containing original
#' read width.
#'
#' If input is p-shifted and GRanges, the "$size" or "$score" colum" must
#' exist, and the column must contain the original read widths. In ORFik
#' "$size" have higher priority than "$score" for defining length.
#' ORFik P-shifting creates a $size column, other softwares like shoelaces
#' creates a score column.
#'
#' Remember to think about how you define length. Like the question:
#' is a Illumina error mismatch sufficient to reduce size of read and how
#' do you know what is biological variance and what are Illumina errors?
#' @param reads a GRanges or GAlignment object.
#' @param after.softclips logical (TRUE), include softclips in width. Does not
#' apply if along.reference is TRUE.
#' @param along.reference logical (FALSE), example: The cigar "26MI2" is
#' by default width 28, but if along.reference is TRUE, it will be 26.
#' The length of the read along the reference. Also "1D20M" will be
#' 21 if by along.reference is TRUE. Intronic regions (cigar: N) will
#' be removed. So: "1M200N19M" is 20, not 220.
#' @return an integer vector of widths
#' @export
#' @examples
#' gr <- GRanges("chr1", 1)
#' readWidths(gr)
#'
#' # GAlignment with hit (1M) and soft clipped base (1S)
#' ga <- GAlignments(seqnames = "1", pos = as.integer(1), cigar = "1M1S",
#'  strand = factor("+", levels = c("+", "-", "*")))
#' readWidths(ga) # Without soft-clip bases
#'
#' readWidths(ga, after.softclips = FALSE) # With soft-clip bases
#'
readWidths <- function(reads, after.softclips = TRUE, along.reference = FALSE) {

  if (is(reads, "GRanges")) {
    readWidth <- width(reads)
    is.one_based <- all(as.integer(readWidth) == rep(1, length(readWidth)))
    if (is.one_based ) {
      if (is.null(reads$size)) {
        if (is.null(reads$score)) {
          message("All widths are 1, If ribo-seq is p-shifted, ",
                  "score or size meta column should contain widths of read, ",
                  "will continue using 1-widths")
        } else {
          message("All widths are 1, using score column for widths, remove ",
                  "score column and run again if this is wrong.")
          readWidth <- reads$score
        }

      } else {
        readWidth <- reads$size
      }
    }
  } else {
    readWidth <- if (along.reference) {
      cigarWidthAlongReferenceSpace(cigar(reads), N.regions.removed = TRUE)
    } else {
      cigarWidthAlongQuerySpace(cigar(reads),
                                after.soft.clipping =
                                  after.softclips)
    }
  }
  return(readWidth)
}

#' Get weights from a subject GenomicRanges object
#' @param subject a GRanges, IRanges or GAlignment object
#' @param weight a vector (default: 1L, if 1L it is identical to
#' countOverlaps()),
#' if single number (!= 1), it applies for all,
#' if more than one must be equal size of 'reads'.
#' else it must be the string name of a defined meta column in subject
#' "reads", that gives number of times a read was found.
#' GRanges("chr1", 1, "+", score = 5), would mean "score" column tells
#' that this alignment region was found 5 times.
#' @return a numeric vector of weights of equal size to subject
getWeights <- function(subject, weight = 1L) {
  weight <- if (is.numeric(weight)) {
    if (length(weight) == 1) {
      rep.int(weight, length(subject))
    } else if (length(weight) != length(subject)) {
      stop("weight does not have equal length to subject")
    } else weight
  } else if (is.character(weight)) {
    error <- paste("weights:", weight, "is not a mcol in 'reads'")
    if (!(weight %in% colnames(mcols(subject)))) stop(error)
    mcols(subject)[, weight]
  } else stop("weight must be numeric or character",
              "name of valid mcol in subject")
  return(weight)
}

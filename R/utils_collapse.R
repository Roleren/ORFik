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
#' @param method character, default \code{"5prime"},
#' the method to reduce ranges, see NOTE for more info.
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
#' @keywords internal
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
  return(makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE, seqinfo = seqinfo(x)))
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
#' only collapse and not keep score column of counts for collapsed reads.
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
            return(makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE, seqinfo = seqinfo(x)))
          })

#' @inherit collapseDuplicatedReads
#' @param addScoreColumn = TRUE, if FALSE,
#' only collapse and not keep score column of counts for collapsed reads.
#' Returns directly without collapsing if reuse.score.column is FALSE and
#' score is already defined.
#' @inheritParams convertToOneBasedRanges
setMethod("collapseDuplicatedReads", "GAlignments",
          function(x, addScoreColumn = TRUE, reuse.score.column = TRUE) {
            if (("score" %in% colnames(mcols(x))) & !reuse.score.column) return(x)

            dt <- data.table(seqnames = factor(seqnames(x)),
                             start = start(ranges(x)),
                             cigar = cigar(x),
                             strand = factor(strand(x)))
            if (reuse.score.column & ("score" %in% colnames(mcols(x)))) { # reuse
              dt[, score := mcols(x)$score]
              dt <- dt[, .(score = sum(score)), .(seqnames, start, cigar, strand)]
            } else { # Do not reuse or "score" does not exist
              dt <- dt[, .(score = .N), .(seqnames, start, cigar, strand)]
            }

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

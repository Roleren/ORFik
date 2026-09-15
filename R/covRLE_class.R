#' Coverage Rle for both strands or single
#'
#' Given a run of coverage(x) where x are reads,
#' this class combines the 2 strands into 1 object
#' @importFrom methods new
#' @return a covRLE object
#' @export
#' @family covRLE
setClass("covRle",
         contains="RleList",
         representation(
           forward="RleList",          # of length N, no names
           reverse="RleList",           # of length N, no names
           strandMode="integer"
         ),
         prototype(
           forward="RleList",
           reverse="RleList",
           strandMode=integer()
         )
)

#' List of covRle
#'
#' Given a run of coverage(x) where x are reads,
#' this covRle combines the 2 strands into 1 object
#' This list can again combine these into 1 object, with accession functions
#' and generalizations.
#' @importFrom methods new
#' @export
#' @return a covRleList object
#' @family covRLE
setClass("covRleList",
         contains=c("List"),
         representation(
           list="List",
           strandMode="integer",
           fraction="character"
         ),
         prototype(
           list=List(),
           strandMode=integer(),
           fraction=character()
         )
)

#' Coverage Rlelist for both strands
#'
#' @param forward a RleList with defined seqinfo for forward strand counts
#' @param reverse a RleList with defined seqinfo for reverse strand counts
#' @return a covRle object
#' @export
#' @family covRLE
#' @examples
#' covRle()
#' covRle(RleList(), RleList())
#' chr_rle <- RleList(chr1 = Rle(c(1,2,3), c(1,2,3)))
#' covRle(chr_rle, chr_rle)
covRle <- function(forward = RleList(), reverse = RleList()) {
  strandMode <- ifelse(length(reverse) > 0, 1L, 0L)
  if (strandMode) {
    if (length(forward) != length(reverse))
      stop("Number of chromosomes must match (length)")
    if (any(lengths(forward) != lengths(reverse)))
      stop("Lengths of each chromosome must match (lengths)")
  }
  new("covRle", forward = forward, reverse = reverse, strandMode = strandMode)
}

#' Coverage Rlelist for both strands
#'
#' @param list a list or List of covRle objects of equal length and lengths
#' @param fraction character, default \code{names(list)}.
#' Names to elements of list, can be integers, as readlengths etc.
#' @return a covRleList object
#' @export
#' @family covRLE
#' @examples
#' covRleList(List(covRle()))
covRleList <- function(list, fraction = names(list)) {
  if (is.null(fraction)) {
    fraction <- rep("NA", length(list))
  }
  fraction <- as.character(fraction)

  strandMode <- as.integer(NA)
  if (length(list) > 0) {
    if (is.null(names(list))) {
      names(list) <- fraction
    }
    if (!is(list[[1]], "covRle")) stop("'list' must only contain covRle objects")
    strandMode <- ifelse(length(r(list[[1]])) > 0, 1L, 0L)
  }

  if (length(list) > 1) {
    stopifnot(all(unlist(lapply(list, function(x) is(x, "covRle")))))
  }
  new("covRleList", list = List(list), strandMode = strandMode, fraction = fraction)
}

#' covRle show definition
#'
#' Show a simplified version of the covRle
#' @param object a\code{\link{covRle}}
#' @export
#' @return print state of covRle
setMethod("show", "covRle",
          function(object) {
            if (length(object) == 0) {
              cat(paste0("CovRle object with 0 sequences"))
            } else {
              if (strandMode(object)) {
                cat(paste0("CovRle object with ", length(object), " sequences (Stranded)\n"))
                print(c(forward = object@forward, reverse = object@reverse))
              } else {
                cat(paste0("CovRle object with ", length(object@forward), " sequences (Unstranded)\n"))
                print(object@forward)
              }
            }
          }
)

#' covRleList show definition
#'
#' Show a simplified version of the covRleList.
#' @param object a\code{\link{covRleList}}
#' @export
#' @return print state of covRleList
setMethod("show", "covRleList",
          function(object) {
            if (length(object) == 0) {
              cat(paste0("CovRleList object with 0 covRle objects\n"))
            } else {
              if (is.na(strandMode(object))) {
                cat(paste0("CovRleList object with ", length(object)),
                " covRle objects\n")
              } else if (strandMode(object)) {
                  cat(paste0("CovRleList object with ", length(object),
                             " covRle objects (Stranded)\n"))
              } else {
                  cat(paste0("CovRleList object with ", length(object),
                             " covRle objects (Unstranded)\n"))
              }
              if (length(object@list) > 3) {
                print(head(as.list(object@list)))
                cat("...\n")
                cat(paste0("<", length(object@list) - 3,
                           " more covRle elements>\n"))
              } else print(as.list(object@list))
            }
          }
)

#' Seqlevels covRle
#' Extracted from forward RleList
#' @param x a covRle object
#' @return integer vector with names
#' @export
setMethod("seqlevels",
          "covRle",
          function(x) {
            seqlevels(x@forward)
          }
)

#' Seqinfo covRle
#' Extracted from forward RleList
#' @param x a covRle object
#' @return integer vector with names
#' @export
setMethod("seqinfo",
          "covRle",
          function(x) {
            seqinfo(x@forward)
          }
)

#' Seqlevels covRleList
#' Extracted from forward RleList
#' @param x a covRle object
#' @return integer vector with names
#' @export
setMethod("seqlevels",
          "covRleList",
          function(x) {
            if (length(x) != 0) {
              seqlevels(x@list[[1]])
            } else character()
          }
)

#' Seqinfo covRle
#' Extracted from forward RleList
#' @param x a covRle object
#' @return integer vector with names
#' @export
setMethod("seqinfo",
          "covRleList",
          function(x) {
            if (length(x) != 0) {
              seqinfo(x@list[[1]])
            } else Seqinfo()
          }
)

#' strandMode covRle
#' @param x a covRle object
#' @return integer vector with names
#' @export
setMethod("strandMode",
          "covRle",
          function(x) {
            x@strandMode
          }
)

#' strandMode covRle
#' @param x a covRle object
#' @return integer vector with names
#' @export
setMethod("strandMode",
          "covRleList",
          function(x) {
            x@strandMode
          }
)

#' strandMode covRle
#' @param x a covRle object
#' @return the forward RleList
#' @export
setGeneric("f", function(x) standardGeneric("f"))


#' @inherit f
setMethod("f",
          "covRle",
          function(x) {
            x@forward
          }
)

#' strandMode covRle
#' @param x a covRle object
#' @return the forward RleList
#' @export
setGeneric("r", function(x) standardGeneric("r"))

#' @inherit r
setMethod("r",
          "covRle",
          function(x) {
            x@reverse
          }
)

#' lengths covRle
#'
#' Lengths of each chromosome
#' @param x a covRle object
#' @return a named integer vector of chromosome lengths
#' @export
setMethod("lengths",
          "covRle",
          function(x) {
            lengths(x@forward)
          }
)

#' lengths covRleList
#'
#' Lengths of each chromosome
#' @param x a covRle object
#' @return a named integer vector of chromosome lengths
#' @export
setMethod("lengths",
          "covRleList",
          function(x) {
            if (length(x) > 0) {
              lengths(x@list[[1]])
            } else {
              0L
            }
          }
)

#' sum covRle
#'
#' Sum coverage per chromosome
#' @param x a covRle object
#' @return an integer, sum of coverage per chromosomes in covRle object
#' @export
setMethod("sum",
          "covRle",
          function(x) {
            f <- sum(f(x))
            r <- sum(r(x))
            return(matrix(data = c(f, r), ncol = 2,
                          dimnames = list(names(f), c("f", "r"))))
          }
)

#' length covRle
#'
#' Number of chromosomes
#' @param x a covRle object
#' @return an integer, number of chromosomes in covRle object
#' @export
setMethod("length",
          "covRle",
          function(x) {
            length(x@forward)
          }
)

#' length covRleList
#'
#' Number of covRle objects
#' @param x a covRleList object
#' @return an integer, number of covRle objects
#' @export
setMethod("length",
          "covRleList",
          function(x) {
            length(x@list)
          }
)

setMethod("countOverlaps",
          c("GRangesList", "covRle"),
          function(query, subject) {
            if (length(query) > 5000) {
              coverageByTranscriptSum(subject, query)
            } else {
              sum(coveragePerTiling(query, subject, is.sorted = TRUE))
            }
          }
)

setMethod("countOverlaps",
          c("GRanges", "covRle"),
          function(query, subject) {
            old_query_names <- names(query)
            query <- split(query, seq(length(query)))
            names(query) <- old_query_names
            countOverlaps(query, subject) # Forward as GRangesList
          }
)

setMethod("countOverlaps",
          c("GRangesList", "character"),
          function(query, subject) {
            sum(coveragePerTiling(query, subject, is.sorted = TRUE))
          }
)

setMethod("countOverlaps",
          c("GRanges", "character"),
          function(query, subject) {
            old_query_names <- names(query)
            query <- split(query, seq(length(query)))
            names(query) <- old_query_names
            countOverlaps(query, subject) # Forward as GRangesList
          }
)



#' Check for integer overflow in covRle
#' If any sum is > 2^31-1, it will give NA, convert those Rles to numeric
#' @noRd
check_for_na_covRle <- function(cov, x, weight) {

  cov@forward <- check_for_na_coverage(f(cov), x, weight)
  cov@reverse <- check_for_na_coverage(r(cov), x, weight)
  return(cov)
}

#' Check for integer overflow in RleList
#' If any sum is > 2^31-1, it will give NA, convert those Rles to numeric
#' @noRd
check_for_na_coverage <- function(cov, x, weight) {
  if (sum(lengths(cov)) == 0) return(cov)

  had_names <- ifelse(is.null(names(RleList)), TRUE, FALSE)
  if (!had_names) names(cov) <- unique(as.character(seqnames(x)))
  which_na <- which(is.na(sum(cov)))
  if (length(which_na) > 0) {
    message("ORFik internally fixes integer overflow, ignore overflow warning")
    chr_na <- names(which_na)
    x_numeric <- x[as.character(seqnames(x)) %in% chr_na]
    if (is.character(weight)) {
      weights <- as.numeric(mcols(x_numeric)[[weight]])
    }
    cov_subset_numeric <- coverage(x_numeric, weight = weights)
    cov[chr_na] <- cov_subset_numeric[chr_na]
  }
  return(cov)
}


#' Convert GRanges to covRle
#' @param x a GRanges, GAlignment or GAlignmentPairs object.
#' Note that coverage calculation for GAlignment is slower, so usually best
#' to call convertToOneBasedRanges on GAlignment object to speed it up.
#' @param weight default "AUTO", pick 'score' column if exist, else all are 1L.
#' Can also be a manually assigned meta column like 'score2' etc.
#' @param ignore.strand logical, default FALSE.
#' @param chunk.size NULL (default) tries coverage directly, then retries in
#' smaller alignment chunks on recognized expansion/allocation limit errors.
#' A positive integer forces chunks of at most this many alignments (or pairs).
#' Can also be set globally with option \code{ORFik.coverage.chunk.size}.
#' Chunked coverage uses numeric counts to avoid integer overflow when adding
#' chunks. No reads are discarded. A single alignment that cannot be expanded,
#' or a final coverage object that cannot fit, still causes an error.
#' Single-end alignments use direct CIGAR-to-IRanges coverage with numeric
#' counts; no intermediate GRangesList is constructed.
#' Chunking may slightly change rounding for fractional weights.
#' @return covRle object
#' @family covRLE
#' @export
#' @examples
#' seqlengths <- as.integer(c(200, 300))
#' names(seqlengths) <- c("chr1", "chr2")
#' gr <- GRanges(seqnames = c("chr1", "chr1", "chr2", "chr2"),
#'                ranges = IRanges(start = c(10, 50, 100, 150), end = c(40, 80, 129, 179)),
#'                strand = c("+", "+", "-", "-"), seqlengths = seqlengths)
#' cov_both_strands <- covRleFromGR(gr)
#' cov_both_strands
#' cov_ignore_strand <- covRleFromGR(gr, ignore.strand = TRUE)
#' cov_ignore_strand
#' strandMode(cov_both_strands)
#' strandMode(cov_ignore_strand)
covRleFromGR <- function(x, weight = "AUTO",
                         ignore.strand = FALSE,
                         chunk.size = getOption("ORFik.coverage.chunk.size", NULL)) {
  if (!is.null(chunk.size) && (!is.numeric(chunk.size) ||
      length(chunk.size) != 1L || is.na(chunk.size) || !is.finite(chunk.size) ||
      chunk.size < 1 || chunk.size > .Machine$integer.max ||
      chunk.size != floor(chunk.size)))
    stop("ORFik coverage: chunk.size must be NULL or a positive integer <= 2^31-1.",
         call. = FALSE)
  if (!(is(x, "GAlignments") || is(x, "GAlignmentPairs")))
    return(.covRleFromGR_once(x, weight, ignore.strand))
  # Validate before retrying: a malformed weight must not become valid merely
  # because its length happens to equal that of a smaller chunk.
  weight <- .coverage_alignment_weight(x, weight)
  if (is.null(chunk.size)) {
    result <- .coverage_try(function() .covRleFromGR_once(x, weight, ignore.strand))
    if (!inherits(result, "error")) return(result)
    message("ORFik coverage: CIGAR expansion/coverage exceeded a size limit; ",
            "retrying in chunks. Original error: ", conditionMessage(result))
    chunk.size <- min(5000000, max(1, floor(length(x) / 2)))
    rm(result)
    gc(FALSE)
  }
  .coverage_alignment_chunks(x, weight, ignore.strand, chunk.size)
}

#' @noRd
.covRleFromGR_once <- function(x, weight, ignore.strand) {
  is_GAlignment <- is(x, "GAlignments") | is(x, "GAlignmentPairs")
  stopifnot(is(x, "GRanges") | is_GAlignment)
  seq_info <- seqinfo(x)
  if (anyNA(seqlengths(seq_info))) stop("Seqlengths of x contains NA values!")
  if (is(x, "GAlignments"))
    return(.coverage_alignment_irl(x, .coverage_alignment_weight(x, weight), ignore.strand))

  # Make sure weight argument is valid for all input types
  if (is.character(weight) & length(weight) == 1) {
    if (weight == "AUTO") {
      if ("score" %in% colnames(mcols(x))) {
        weight <- "score"
      } else weight <- 1L
    }
  }
  if (is_GAlignment) {
    weight <- .coverage_alignment_weight(x, weight)
    x <- grglist(x)
    # Repeat per alignment block, not per base, and avoid an extra full-sized
    # integer grouping-index vector. Also handles manually supplied weights.
    if (length(weight) != 1L) weight <- rep(weight, lengths(x, use.names = FALSE))
  }

  if (ignore.strand) {
    both <- coverage(x, weight = weight)
    seqinfo(both) <- seq_info
    return(check_for_na_covRle(covRle(both), x, weight))
  } else {
    pluss <- BiocGenerics::`%in%`(strand(x), c("+", "*"))
    minus <- BiocGenerics::`%in%`(strand(x), c("-", "*"))
    if (length(weight) > 1) {
      f <- coverage(x[pluss], weight = weight[as.logical(unlist(pluss))])
      r <- coverage(x[minus], weight = weight[as.logical(unlist(minus))])
    } else {
      strands <- strand(x)
      f <- coverage(x[pluss], weight = weight)
      r <- coverage(x[minus], weight = weight)
    }
    seqinfo(f) <- seqinfo(r) <- seq_info
    return(check_for_na_covRle(covRle(f, r), x, weight))
  }
}

covRleListFromGR <- function(x, weight = "AUTO",
                         ignore.strand = FALSE, verbose = TRUE) {
  all_readl_lengths <- readWidths(x)
  read_lengths <- sort(unique(all_readl_lengths))

  if (verbose) message("Readlength:", appendLF = FALSE)
  list <- list()
  for (i in read_lengths) {
    if (verbose) message(", ", i, appendLF = FALSE)
    list <- c(list, covRleFromGR(x[all_readl_lengths == i],
                                 weight = weight,
                                 ignore.strand = ignore.strand))
  }
  return(covRleList(list, fraction = read_lengths))
}

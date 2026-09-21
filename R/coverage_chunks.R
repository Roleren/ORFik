# Internal helpers for coverage of CIGAR expansions beyond short-vector limits.

.coverage_cigar_ranges <- function(x) {
  # Same parser/default deletion handling as GenomicAlignments::rglist().
  GenomicAlignments::extractAlignmentRangesOnReference(cigar(x), start(x))
}

.coverage_alignment_irl <- function(x, weight, ignore.strand) {
  irl <- .coverage_cigar_ranges(x)
  counts <- lengths(irl, use.names = FALSE)
  # Group row indices once instead of rescanning every read for every contig
  # and strand (important for references with many alternate scaffolds).
  rows <- split(seq_along(x), as.character(seqnames(x)))
  s <- as.character(strand(x))
  si <- seqinfo(x)
  # Coverage below has explicit chromosome widths. Summing Seqinfo lengths
  # is O(number of chromosomes), independent of reads and coverage runs.
  compress_list <- sum(as.double(seqlengths(si))) <= .Machine$integer.max
  one_strand <- function(keep) {
    ans <- lapply(seqlevels(si), function(ch) {
      i <- rows[[ch]]
      if (is.null(i)) i <- integer()
      i <- i[keep[i]]
      w <- if (length(weight) == 1L) as.numeric(weight) else
        rep(as.numeric(weight[i]), counts[i])
      coverage(unlist(irl[i], use.names = FALSE),
               width = seqlengths(si)[[ch]], weight = w)
    })
    names(ans) <- seqlevels(si)
    # Above INT_MAX, each chromosome still remains run-length encoded.
    ans <- RleList(ans, compress = compress_list)
    seqinfo(ans) <- si
    ans
  }
  if (ignore.strand) covRle(one_strand(rep(TRUE, length(x)))) else
    covRle(one_strand(s %in% c("+", "*")), one_strand(s %in% c("-", "*")))
}

.coverage_alignment_weight <- function(x, weight) {
  if (is.character(weight) && length(weight) == 1L && !is.na(weight)) {
    if (weight == "AUTO") weight <- if ("score" %in% colnames(mcols(x))) "score" else 1L
    if (is.character(weight)) {
      if (!(weight %in% colnames(mcols(x))))
        stop("ORFik coverage: weight column not found: ", weight, call. = FALSE)
      weight <- mcols(x)[[weight]]
    }
  }
  if (!(is.integer(weight) || is.double(weight)) || is.object(weight) ||
      !is.null(dim(weight)) || !(length(weight) %in% c(1L, length(x))))
    stop("ORFik coverage: alignment weights must be numeric, with length 1 or one per alignment/pair.",
         call. = FALSE)
  weight
}

.coverage_size_error <- function(e) {
  grepl(paste(c("long vectors? (are )?not supported", "negative length vectors",
                "cannot allocate (vector|memory)",
                "subsetting a Vector derivative of length 2\\^31",
                "cumulated length of its list elements",
                "too many (ranges|elements)",
                "maximum.*buffer.*size", "buffer.*size.*maximum"), collapse = "|"),
        gsub("[[:space:]]+", " ", conditionMessage(e)), ignore.case = TRUE)
}

.coverage_try <- function(fun) {
  tryCatch(fun(), error = function(e) {
    if (.coverage_size_error(e)) return(e)
    stop(e)
  })
}

.coverage_alignment_chunks <- function(x, weight, ignore.strand, chunk.size) {
  n <- length(x)
  if (!n) return(.covRleFromGR_once(x, weight, ignore.strand))
  from <- 1
  answer <- NULL
  while (from <= n) {
    to <- min(n, from + chunk.size - 1)
    # Keep failed chunk allocations out of this frame: the error handler runs
    # after unwinding .covRleFromGR_once, so a retry need not retain its ranges.
    part <- .coverage_try(function() {
      i <- seq.int(from, to)
      w <- if (length(weight) == 1L) weight else weight[i]
      .covRleFromGR_once(x[i], as.numeric(w), ignore.strand)
    })
    if (inherits(part, "error")) {
      if (to == from)
        stop("ORFik coverage: cannot expand even one alignment/pair at index ",
             from, "; chunking cannot resolve this limit. ",
             conditionMessage(part), call. = FALSE)
      chunk.size <- max(1, floor((to - from + 1) / 2))
      message("ORFik coverage: reducing chunk size to ", chunk.size,
              " alignments/pairs after size-limit error: ", conditionMessage(part))
      rm(part)
      gc(FALSE)
      next
    }
    if (is.null(answer)) answer <- part else {
      # Numeric weights make both operands numeric Rles: adding large counts
      # must not silently overflow at 2^31-1. Retain the original Seqinfo.
      answer@forward <- answer@forward + part@forward
      if (!ignore.strand) answer@reverse <- answer@reverse + part@reverse
    }
    message("ORFik coverage: completed alignments/pairs ", from, "-", to, " of ", n, ".")
    rm(part)
    from <- to + 1
  }
  seqinfo(answer@forward) <- seqinfo(x)
  if (!ignore.strand) seqinfo(answer@reverse) <- seqinfo(x)
  answer
}

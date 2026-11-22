# bigWig export of coverage (RleList) with ORFik's C++ writer (src/bigwig_writer.cpp).
#
# Equivalent to rtracklayer::export.bw(GRanges(<RleList>)): every run, zero runs included, is one
# bedGraph item, so rtracklayer::import.bw() returns the same GRanges, RleList and NumericList.
# Data blocks are deflated in parallel (OpenMP); the result is about 7x faster at 8 threads.
# Sequences keep the seqlevels order (as rtracklayer); there is no limit on their number.
# Zoom levels are fixed (256 bp * 4^k, k = 0..6), so zoomed-out summaries can differ slightly
# from rtracklayer's files; full-resolution values are identical.

# Default number of threads for the bigWig writer: option ORFik.bigwig_threads, else the
# data.table thread count (which honours OMP limits and R CMD check), at most 8.
.bigwig_threads <- function() {
  th <- getOption("ORFik.bigwig_threads", NULL)
  if (!is.null(th)) return(max(1L, as.integer(th)))
  max(1L, min(8L, data.table::getDTthreads()))
}

#' Write one RleList as a bigWig file with ORFik's writer
#' @param x RleList (named, one Rle per sequence)
#' @param file output path
#' @param threads integer, threads for block compression
#' @param compression_level integer 0-9, zlib level (6 as rtracklayer)
#' @param seq_info Seqinfo or NULL; sequences of the file (default: names and lengths of x).
#' Sequences of seq_info that are missing in x are written without data.
#' @return invisible numeric c(size, items, blocks)
#' @noRd
.export_bw_rlelist <- function(x, file, threads = .bigwig_threads(), compression_level = 6L,
                               seq_info = NULL) {
  stopifnot(is(x, "RleList"))
  rles <- as.list(x)
  nm <- names(rles)
  if (is.null(nm) || anyNA(nm) || any(!nzchar(nm)) || anyDuplicated(nm))
    stop("bigWig export: coverage must have unique sequence names")
  sizes <- lengths(rles)
  if (!is.null(seq_info)) {
    extra <- setdiff(seqlevels(seq_info), nm)
    if (length(extra)) {
      nm <- c(nm, extra); sizes <- c(sizes, seqlengths(seq_info)[extra])
      rles <- c(rles, vector("list", length(extra)))
    }
  }
  if (anyNA(sizes) || any(sizes > .Machine$integer.max)) stop("bigWig export: invalid sequence lengths")
  # chromosome ids in seqlevels order, as rtracklayer (one B+ tree leaf); above 65,535 sequences
  # the tree has several levels and needs bytewise name order
  ord <- if (length(nm) > 65535L) order(nm, method = "radix") else seq_along(nm)
  values <- lapply(rles[ord], function(r) {
    if (is.null(r)) return(NULL)
    v <- S4Vectors::runValue(r)
    if (!is.integer(v) && !is.double(v)) v <- as.numeric(v)
    if (anyNA(v)) stop("bigWig export: coverage has NA values")
    v
  })
  run_lengths <- lapply(rles[ord], function(r) if (is.null(r)) NULL else S4Vectors::runLength(r))
  dir <- dirname(file)
  if (!dir.exists(dir)) stop("bigWig export: directory does not exist: ", dir)
  tmp <- tempfile(paste0(".", basename(file), "_"), tmpdir = dir)
  on.exit(unlink(tmp), add = TRUE)
  res <- bw_write_runs(path.expand(tmp), nm[ord], as.integer(sizes[ord]), values, run_lengths,
                       threads = as.integer(threads), level = as.integer(compression_level),
                       zero_runs = TRUE)
  if (!file.rename(tmp, file)) stop("bigWig export: cannot write ", file)
  invisible(res)
}

# One coverage track to bigWig with the chosen writer.
.export_bw_track <- function(x, file, writer, threads, compression_level) {
  if (writer == "rtracklayer") {
    rtracklayer::export.bw(GRanges(x), file)
  } else .export_bw_rlelist(x, file, threads = threads, compression_level = compression_level)
  invisible(NULL)
}

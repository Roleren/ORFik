# Single-pass derived formats of a sorted point OFST file (merged master tracks, P-site
# libraries): covRleList, covRle, bigWig and count tables from one read of each chromosome block.
#
# The file is read in waves of whole chromosome blocks (sidecar index <file>.idx, else one scan
# of seqnames). ofst_rle_blocks() (src/ofst_derive.cpp) builds the Rle runs of every
# (read length, strand) group and of all lengths merged, in parallel over position ranges of
# the blocks, and the table chunk is dropped before the next wave. The merged runs go straight
# to the bigWig writer, which deflates them on a background thread while the next waves are
# read. Count tables are range sums over the merged runs. Outputs are identical to
# .convert_point_ofst_cov() (old path), export.bigWig(<covRle>) and countTable_regions(lib.type
# = "cov"). The kernel validates every row (sorted, positions inside the chromosome, known
# strand, read lengths, no NA score), so the R validation passes of .ofst_point_table() are not
# needed; any violation makes ofst_derive() return FALSE before anything is written, and
# callers use the old path.

# Default threads: option ORFik.derive_threads, else the data.table thread count, at most 16.
.derive_threads <- function() {
  th <- getOption("ORFik.derive_threads", NULL)
  if (!is.null(th)) return(max(1L, as.integer(th)))
  max(1L, min(16L, data.table::getDTthreads()))
}

#' Derived formats of a sorted point OFST file in one pass
#'
#' Builds the coverage formats of a point OFST file (one position per row, e.g. P-sites of a
#' merged master track, with columns seqnames, start, strand, score and optionally size)
#' from one read of each chromosome block: the per read length coverage (covRleList), the
#' merged coverage (covRle), the bigWig files of both strands and count tables over regions.
#' Chromosome blocks come from the sidecar index (\code{<in_file>.idx}) when present, so the
#' file is never scanned. Work is parallel over position ranges of the chromosomes, and the
#' bigWig files are compressed on a background thread while the next chromosomes are read.
#'
#' The outputs are identical to the separate steps: covRleList and covRle files as written by
#' \code{convertLibs(type = "covl")} for point data, the bigWig files of
#' \code{\link{export.bigWig}} on the covRle, and the counts of
#' \code{\link{countTable_regions}} with \code{lib.type = "cov"}.
#' @param in_file path to a point OFST file, rows sorted by chromosome blocks and start
#'  (merge and collapse output order).
#' @param seq_info Seqinfo with the lengths of all chromosomes; the coverage covers all of them.
#' @param covlist_file path of the covRleList file (.covqs / .covrds) or NULL to skip it.
#' @param cov_file path of the covRle file or NULL to skip it.
#' @param bigwig_file NULL to skip, else the bigWig prefix: files
#'  \code{<prefix>_forward.bigWig} and \code{<prefix>_reverse.bigWig} are written.
#' @param regions NULL or a named list of GRangesList (e.g. mrna, cds): their count vectors
#'  are returned.
#' @param threads integer, threads for the coverage kernel and the fst reads.
#' @param bigwig_threads integer, threads of each bigWig writer (two run at once, mostly while
#'  the coverage kernel waits for reads or after it finished).
#' @param compression_level integer 0-9, zlib level of the bigWig files (6 as rtracklayer).
#' @param qs_threads integer, threads for the qs2 saves.
#' @param wave_rows numeric, about how many rows are read at once (whole chromosome blocks;
#'  a larger block is read alone). Lower it to save memory.
#' @param task_rows numeric, rows per parallel task inside a chromosome block.
#' @return FALSE (nothing written) if the file is not sorted point data or a chromosome is
#'  missing from seq_info, so callers can use the general converters; else invisibly a list
#'  with \code{sizes} (read lengths), \code{counts} (named list of count vectors, one per
#'  region) and \code{seconds} (time per step).
#' @export
#' @examples
#' # ofst_derive("all.ofst", seqinfo(df), "cov_RLE_List/all.covqs", "cov_RLE/all.covqs",
#' #             "bigwig/all", regions = list(cds = loadRegion(df, "cds")))
ofst_derive <- function(in_file, seq_info, covlist_file = NULL, cov_file = NULL, bigwig_file = NULL,
                        regions = NULL, threads = .derive_threads(), bigwig_threads = threads,
                        compression_level = 6L, qs_threads = threads, wave_rows = 3e8, task_rows = 2e6) {
  secs <- list(); t0 <- proc.time()[[3]]
  tick <- function(name) {   # seconds per step, summed over the waves
    t1 <- proc.time()[[3]]; secs[[name]] <<- (if (is.null(secs[[name]])) 0 else secs[[name]]) + t1 - t0; t0 <<- t1 }
  if (!is.null(regions) && (!is.list(regions) || is.null(names(regions)))) stop("regions must be a named list")
  if (length(in_file) != 1L || !grepl("\\.ofst$", in_file, ignore.case = TRUE) || !file.exists(in_file)) return(FALSE)
  meta <- fst::metadata_fst(in_file)
  cols <- meta$columnNames
  if (!all(c("seqnames", "start", "strand", "score") %in% cols) ||
      any(c("cigar", "cigar1", "cigar2", "start1", "start2", "end") %in% cols)) return(FALSE)
  if (meta$nrOfRows == 0 || anyNA(seqlengths(seq_info))) return(FALSE)
  chrs <- seqlevels(seq_info)
  lens <- as.integer(seqlengths(seq_info))
  blocks <- .ofst_chr_blocks(in_file)
  if (is.null(blocks) || !nrow(blocks) || anyNA(blocks$chr) || !all(blocks$chr %in% chrs)) return(FALSE)
  first <- fst::read_fst(in_file, intersect(c("strand", "size", "score"), cols), from = 1, to = 1)
  if (!is.factor(first$strand) || !(is.integer(first$score) || is.double(first$score)) ||
      ("size" %in% cols && !is.integer(first$size))) return(FALSE)
  strand_map <- match(levels(first$strand), c("+", "-", "*"))
  integer_weights <- is.integer(first$score)
  zero <- if (integer_weights) 0L else 0
  read_cols <- intersect(c("start", "strand", "size", "score", "width"), cols)
  blocks[, ci := match(chr, chrs)]
  tick("blocks")

  old_fst <- fst::threads_fst(max(1L, as.integer(threads)))
  on.exit(fst::threads_fst(old_fst), add = TRUE)

  # bigWig writers (forward, reverse): sequence ids in seqlevels order, as export.bigWig()
  bw <- NULL
  if (!is.null(bigwig_file)) {
    if (any(lens > .Machine$integer.max)) return(FALSE)
    ord <- if (length(chrs) > 65535L) order(chrs, method = "radix") else seq_along(chrs)
    bw_id <- integer(length(chrs)); bw_id[ord] <- seq_along(ord)
    prefix <- gsub("\\.bigWig$", "", bigwig_file, ignore.case = TRUE)
    bw <- lapply(c(forward = "_forward.bigWig", reverse = "_reverse.bigWig"), function(sfx) {
      file <- paste0(prefix, sfx)
      if (!dir.exists(dirname(file))) stop("bigWig export: directory does not exist: ", dirname(file))
      tmp <- tempfile(paste0(".", basename(file), "_"), tmpdir = dirname(file))
      list(file = file, tmp = tmp,
           ptr = bw_stream_open(path.expand(tmp), chrs[ord], lens[ord], threads = as.integer(bigwig_threads),
                                level = as.integer(compression_level), zero_runs = TRUE, async = TRUE))
    })
    # on any exit before both files are closed: stop the writers and remove the partial files
    on.exit(for (b in bw) { bw_stream_abort(b$ptr); unlink(b$tmp) }, add = TRUE)
  }
  bw_add <- function(i, fwd, rev) {
    if (is.null(bw)) return(invisible())
    bw_stream_add(bw$forward$ptr, bw_id[i], fwd$values, fwd$lengths)
    bw_stream_add(bw$reverse$ptr, bw_id[i], rev$values, rev$lengths)
  }
  zero_runs <- function(i) list(values = zero, lengths = lens[i])

  # waves: consecutive blocks of the file, about wave_rows rows each
  n_rows <- blocks$to - blocks$from + 1
  wave <- integer(nrow(blocks)); w <- 1L; acc <- 0
  for (b in seq_len(nrow(blocks))) {
    if (acc > 0 && acc + n_rows[b] > wave_rows) { w <- w + 1L; acc <- 0 }
    wave[b] <- w; acc <- acc + n_rows[b]
  }
  per_chr <- vector("list", length(chrs))   # list(sizes, runs = list of group runs or NULL)
  for (w in unique(wave)) {
    bl <- blocks[wave == w]
    from <- bl$from[1]; to <- bl$to[nrow(bl)]
    dt <- fst::read_fst(in_file, read_cols, from = from, to = to)
    tick("read")
    if (!is.null(dt$width) && (anyNA(dt$width) || any(dt$width != 1L))) return(FALSE)
    res <- ofst_rle_blocks(dt$start, dt$strand, strand_map, dt$size, dt$score,
                           bl$from - from, bl$to - from + 1, lens[bl$ci],
                           threads = as.integer(threads), task_rows = task_rows)
    dt <- NULL
    if (is.null(res)) return(FALSE)
    tick("kernel")
    S <- length(res$sizes)
    for (b in seq_len(nrow(bl))) {
      i <- bl$ci[b]
      runs <- res$blocks[[b]]
      per_chr[[i]] <- list(sizes = res$sizes, runs = runs)
      bw_add(i, if (is.null(runs[[2L * S + 1L]])) zero_runs(i) else runs[[2L * S + 1L]],
                if (is.null(runs[[2L * S + 2L]])) zero_runs(i) else runs[[2L * S + 2L]])
    }
    res <- NULL
    tick("bigwig_queue")
  }
  for (i in which(vapply(per_chr, is.null, logical(1)))) bw_add(i, zero_runs(i), zero_runs(i))

  # Rle objects of every group, in the group order of ofst_point_sums()
  sizes <- sort(unique(unlist(lapply(per_chr, `[[`, "sizes"), use.names = FALSE)))
  n_size <- length(sizes)
  local_index <- lapply(per_chr, function(p) if (is.null(p)) NULL else match(sizes, p$sizes))
  rle_of <- function(i, g) {   # g: global group (2 * n_size + 1/2 = merged)
    p <- per_chr[[i]]
    if (is.null(p)) return(S4Vectors::Rle(zero, lens[i]))
    S <- length(p$sizes)
    k <- (g + 1L) %/% 2L; s <- 2L - g %% 2L      # size group, strand (1 forward, 2 reverse)
    lg <- if (k > n_size) 2L * S + s else { kl <- local_index[[i]][k]; if (is.na(kl)) NA else 2L * kl - 2L + s }
    r <- if (is.na(lg)) NULL else p$runs[[lg]]
    # kernel runs are normalised (no empty runs, no equal neighbours): no Rle() re-check, same object
    if (is.null(r)) S4Vectors::Rle(zero, lens[i]) else S4Vectors::new2("Rle", values = r$values, lengths = r$lengths, check = FALSE)
  }
  groups <- lapply(seq_len(2L * (n_size + 1L)), function(g) {
    rles <- lapply(seq_along(chrs), rle_of, g = g)
    names(rles) <- chrs
    IRanges::RleList(rles, compress = FALSE)
  })
  per_chr <- NULL
  build_g <- function(k) {
    f <- groups[[2L * k - 1L]]; r <- groups[[2L * k]]
    seqinfo(f) <- seqinfo(r) <- seq_info
    covRle(f, r)
  }
  merged <- build_g(n_size + 1L)
  tick("rle_objects")

  counts <- NULL
  if (!is.null(regions)) {
    counts <- .cov_regions_counts(regions, merged)
    tick("counts")
  }
  if (!is.null(covlist_file)) {
    covs <- lapply(seq_len(n_size), build_g)
    save_RDSQS(covRleList(covs, fraction = sizes), file = covlist_file, nthread = qs_threads)
    covs <- NULL
  }
  groups <- NULL
  if (!is.null(cov_file)) {
    seqinfo(merged) <- seq_info   # as export.cov
    save_RDSQS(merged, file = cov_file, nthread = qs_threads)
  }
  tick("qs_save")
  if (!is.null(bw)) {
    for (b in bw) {
      bw_stream_close(b$ptr)
      if (!file.rename(b$tmp, b$file)) stop("bigWig export: cannot write ", b$file)
    }
    tick("bigwig_close")
  }
  invisible(list(sizes = sizes, counts = counts, seconds = unlist(secs)))
}

# Counts of regions (named list of GRangesList) on a covRle, each identical to
# countOverlapsW(region, cov): per exon sums of the strand coverage ("*" exons count 0), summed
# per transcript. Range sums over the runs (rle_range_sums, one sweep per chromosome strand for
# the exons of all regions); exact for whole-number coverage. The result type follows
# coverageByTranscriptSum(): double if any exon of the region lies on a double coverage
# strand or its sum leaves the integer range (viewSums then returns double), else integer.
# Anything else (regions of <= 5000 transcripts, fractional coverage, integer totals above the
# integer range, exons outside the coverage) uses countOverlapsW.
.cov_regions_counts <- function(regions, cov) {
  slow <- function(tx) countOverlapsW(tx, cov, weight = NULL)
  fast <- vapply(regions, function(tx) length(tx) > 5000 && length(unlist(tx, use.names = FALSE)) > 0, logical(1))
  if (!strandMode(cov) || identical(attr(cov, "rna_strand"), "unstranded")) fast[] <- FALSE
  out <- vector("list", length(regions)); names(out) <- names(regions)
  if (any(fast)) out[fast] <- .cov_regions_counts_fast(regions[fast], cov)
  for (i in which(!fast | vapply(out, is.null, logical(1)))) out[[i]] <- slow(regions[[i]])
  out
}

# list of count vectors, NULL elements where the fast path does not apply
.cov_regions_counts_fast <- function(regions, cov) {
  fw <- f(cov); rv <- r(cov)
  ex <- lapply(regions, unlist, use.names = FALSE)
  cm <- unlist(lapply(ex, function(x) match(seqlevels(x), names(fw))[as.integer(seqnames(x))]), use.names = FALSE)
  st <- unlist(lapply(ex, function(x) as.integer(strand(x))), use.names = FALSE)   # levels "+", "-", "*"
  s <- unlist(lapply(ex, start), use.names = FALSE); e <- unlist(lapply(ex, end), use.names = FALSE)
  n_ex <- lengths(ex); ex <- NULL
  none <- vector("list", length(regions))
  keep <- st <= 2L
  if (anyNA(cm[keep])) return(none)
  sums <- numeric(length(st))
  dbl <- logical(length(st))   # exon on a double coverage strand
  idx <- split(which(keep), cm[keep] * 2L + st[keep] - 1L)   # key = 2 * chromosome + strand - 1
  for (key in names(idx)) {
    rows <- idx[[key]]
    k <- as.integer(key)
    x <- if (k %% 2L == 0L) fw[[k %/% 2L]] else rv[[k %/% 2L]]
    v <- S4Vectors::runValue(x)
    if (is.double(v)) {
      if (anyNA(v) || any(v != trunc(v)) || any(abs(v) >= 2^53)) return(none)
      dbl[rows] <- TRUE
    } else if (!is.integer(v)) return(none)
    res <- rle_range_sums(v, S4Vectors::runLength(x), s[rows], e[rows])
    if (is.null(res)) return(none)
    sums[rows] <- res
  }
  region_of <- rep.int(seq_along(regions), n_ex)
  lapply(seq_along(regions), function(i) {
    tx <- regions[[i]]
    in_region <- region_of == i
    x <- sums[in_region]
    as_double <- any(dbl[in_region]) || any(x > .Machine$integer.max)
    ends <- cumsum(as.numeric(lengths(tx)))
    cs <- c(0, cumsum(x))
    tot <- cs[ends + 1] - cs[c(0, ends[-length(ends)]) + 1]
    if (any(abs(cs) >= 2^53)) return(NULL)
    if (!as_double) {
      if (any(tot > .Machine$integer.max)) return(NULL)   # integer sums of today's path wrap around
      tot <- as.integer(tot)
    }
    names(tot) <- names(tx)
    tot
  })
}

#' Derived formats of a one-library experiment (merged master) in one pass
#'
#' Writes the files of convertLibs(type = "covl") / make_point_covrle (cov_RLE_List, cov_RLE),
#' export.bigWig(<covRle>) (bigwig/<lib>_forward/_reverse.bigWig) and
#' countTable_regions(df, lib.type = "cov", forceRemake = TRUE) (QC_STATS/countTable_<region>)
#' with ofst_derive(). Returns FALSE (nothing written) when the library is not sorted point
#' data or is unstranded RNA; callers then run the separate steps.
#' @param df an ORFik experiment with one library (pshifted file = the point OFST)
#' @param regions character, count table regions
#' @param format "qs" or "rds"
#' @param in_file the point OFST file, default the pshifted file of df
#' @param out_dir library folder of the outputs, default libFolder(df)
#' @param ... further arguments to ofst_derive (threads, bigwig_threads, ...)
#' @return FALSE if nothing was written, else invisibly the ofst_derive() result
#' @noRd
ofst_derive_experiment <- function(df, regions = c("mrna", "leaders", "cds", "trailers"), format = "qs",
                                   in_file = filepath(df, "pshifted"), out_dir = libFolder(df), ...) {
  if (nrow(df) != 1) return(FALSE)
  lib_name <- remove.file_ext(basename(in_file))
  folder <- out_dir
  cov_dir <- file.path(folder, "cov_RLE"); covl_dir <- file.path(folder, "cov_RLE_List")
  bw_dir <- file.path(folder, "bigwig")
  if (rna_strand_resolve(in_file) == "unstranded") return(FALSE)
  library.names <- bamVarName(df)
  txdb <- loadTxdb(df)
  tx <- lapply(regions, function(region) loadRegion(txdb, region))
  names(tx) <- regions
  for (d in c(cov_dir, covl_dir, bw_dir)) dir.create(d, showWarnings = FALSE)
  res <- ofst_derive(in_file, seqinfo(df),
                     covlist_file = file.path(covl_dir, paste0(lib_name, ".cov", format)),
                     cov_file = file.path(cov_dir, paste0(lib_name, ".cov", format)),
                     bigwig_file = file.path(bw_dir, lib_name), regions = tx, ...)
  if (isFALSE(res)) return(FALSE)
  t0 <- proc.time()[[3]]
  prefix <- pasteDir(file.path(folder, "QC_STATS", "countTable_"))
  for (region in regions) {   # as makeSummarizedExperimentFromBam()
    rawCounts <- data.table(matrix(0, ncol = length(library.names), nrow = length(tx[[region]])))
    colnames(rawCounts) <- library.names
    rawCounts[, (paste0(library.names)) := res$counts[[region]]]
    se <- SummarizedExperimentExp(df, rawCounts, tx[[region]], library.names)
    saveName <- paste0(prefix, region, ".", format)
    if (!dir.exists(dirname(saveName))) dir.create(dirname(saveName), showWarnings = FALSE, recursive = TRUE)
    save_RDSQS(se, file = saveName)
  }
  res$seconds <- c(res$seconds, count_tables_save = proc.time()[[3]] - t0)
  invisible(res)
}

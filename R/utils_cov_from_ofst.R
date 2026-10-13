# Fast covRle / covRleList for point OFST files (one position per row: P-sites/5' ends, no cigar,
# width absent or 1), e.g. merged master tracks. Coverage is built from positions summed per
# (read length, strand, chromosome) with keyed data.table group-bys and one vectorised Rle per
# chromosome, instead of a GRanges import and one coverage() pass per read length.
# Results are identical to covRleFromGR() / covRleListFromGR() on the imported GRanges.

# Point OFST table with a score column, or NULL if the file is not point data.
.ofst_point_table <- function(path, seq_info) {
  if (length(path) != 1L || !grepl("\\.ofst$", path, ignore.case = TRUE)) return(NULL)
  cols <- fst::metadata_fst(path)$columnNames
  if (!all(c("seqnames", "start", "strand", "score") %in% cols) ||
      any(c("cigar", "cigar1", "cigar2", "start1", "start2", "end") %in% cols)) return(NULL)
  dt <- fst::read_fst(path, intersect(c("seqnames", "start", "strand", "size", "width", "score"), cols),
                      as.data.table = TRUE)
  if ("width" %in% cols) {
    if (anyNA(dt$width) || any(dt$width != 1L)) return(NULL)
    dt[, width := NULL]
  }
  if (!is.factor(dt$seqnames)) data.table::set(dt, j = "seqnames", value = factor(dt$seqnames))
  if (!is.factor(dt$strand)) data.table::set(dt, j = "strand", value = factor(dt$strand))
  if (!is.numeric(dt$score) || anyNA(dt$score) || anyNA(dt$start) || anyNA(dt$seqnames) || anyNA(dt$strand))
    return(NULL)
  if (any(!levels(dt$strand)[tabulate(dt$strand, nlevels(dt$strand)) > 0] %in% c("+", "-", "*"))) return(NULL)
  used <- levels(dt$seqnames)[tabulate(dt$seqnames, nlevels(dt$seqnames)) > 0]
  if (any(!used %in% seqlevels(seq_info)) || anyNA(seqlengths(seq_info))) return(NULL)
  lens <- seqlengths(seq_info)[levels(dt$seqnames)]
  if (any(dt$start < 1L) || any(dt$start > lens[as.integer(dt$seqnames)])) return(NULL)
  dt
}

# RleList over all chromosomes of seq_info from positions summed per chromosome.
# agg: data.table(seqnames = chromosome index into seqlevels(seq_info), start, score), sorted.
.rlelist_from_points <- function(agg, seq_info, integer_weights) {
  chrs <- seqlevels(seq_info)
  lens <- as.integer(seqlengths(seq_info))
  first <- last <- integer(length(chrs))
  if (nrow(agg)) {   # chromosomes are contiguous blocks of agg
    r <- ofst_code_runs(agg$seqnames)
    first[r$value] <- as.integer(c(1, head(r$end, -1) + 1))
    last[r$value] <- as.integer(r$end)
  }
  zero <- if (integer_weights) 0L else 0
  rles <- lapply(seq_along(chrs), function(i) {
    if (!first[i]) return(S4Vectors::Rle(zero, lens[i]))
    rows <- first[i]:last[i]
    p <- agg$start[rows]
    v <- agg$score[rows]
    # integer coverage overflows to NA in coverage(); ORFik then recomputes that chromosome as numeric
    if (integer_weights && all(v <= .Machine$integer.max)) v <- as.integer(v) else v <- as.numeric(v)
    z <- if (is.integer(v)) 0L else 0
    gaps <- diff(c(0L, p)) - 1L
    S4Vectors::Rle(c(as.vector(rbind(rep(z, length(v)), v)), z),
                   c(as.vector(rbind(gaps, 1L)), lens[i] - p[length(p)]))
  })
  names(rles) <- chrs
  rl <- IRanges::RleList(rles, compress = FALSE)
  rl
}

# Sum score per (chromosome, position). Input rows of each chromosome are one contiguous block
# with nondecreasing start (OFST merge/collapse output order): one O(n) pass with cumulative
# sums, no sort. Otherwise a keyed group-by.
.point_sums <- function(chr, start, score) {
  n <- length(start)
  if (!n) return(data.table::data.table(seqnames = integer(), start = integer(), score = numeric()))
  same_chr <- chr[-1L] == chr[-n]
  chr_first <- c(1L, which(!same_chr) + 1L)
  if (!anyDuplicated(chr[chr_first]) && !any(same_chr & start[-1L] < start[-n])) {
    ends <- c(which(!same_chr | start[-1L] != start[-n]), n)
    cs <- cumsum(as.numeric(score))
    return(data.table::data.table(seqnames = chr[ends], start = start[ends], score = diff(c(0, cs[ends]))))
  }
  data.table::data.table(seqnames = chr, start = start, score = as.numeric(score))[
    , .(score = sum(score)), keyby = .(seqnames, start)]
}

# covRle (both strands) of point rows; "*" rows count on both strands, as covRleFromGR().
.covRle_from_points <- function(chr, start, strand, score, strand_lv, seq_info, integer_weights) {
  one <- function(keep) {
    rows <- which(strand %in% which(strand_lv %in% keep))
    .rlelist_from_points(.point_sums(chr[rows], start[rows], score[rows]), seq_info, integer_weights)
  }
  f <- one(c("+", "*"))
  r <- one(c("-", "*"))
  seqinfo(f) <- seqinfo(r) <- seq_info
  covRle(f, r)
}

# Build and save covRleList (per read length) and merged covRle for one point OFST file.
# Returns FALSE (nothing written) when the file is not point data, so callers fall back.
.convert_point_ofst_cov <- function(in_file, out_file, out_file_merged, seq_info, format) {
  dt <- .ofst_point_table(in_file, seq_info)
  if (is.null(dt)) return(FALSE)
  integer_weights <- is.integer(dt$score)
  chr <- match(levels(dt$seqnames), seqlevels(seq_info))[as.integer(dt$seqnames)]
  strand <- as.integer(dt$strand)
  strand_lv <- levels(dt$strand)
  size <- if ("size" %in% names(dt)) dt$size else rep.int(1L, nrow(dt))
  start <- dt$start; score <- dt$score
  dt <- NULL
  if (is.integer(size) && !anyNA(size) && min(size) >= 1L && max(size) <= 1e6) {
    present <- tabulate(size, max(size)) > 0L
    sizes <- which(present)
    size_index <- cumsum(present)[size]
  } else {
    sizes <- sort(unique(size))
    size_index <- match(size, sizes)
  }
  # one C++ pass: summed positions for every (read length, strand) and all lengths merged
  sums <- ofst_point_sums(chr, start, match(strand_lv, c("+", "-", "*"))[strand], size_index, length(sizes), score)
  build <- function(fwd, rev) {
    f <- .rlelist_from_points(data.table::setDT(fwd), seq_info, integer_weights)
    r <- .rlelist_from_points(data.table::setDT(rev), seq_info, integer_weights)
    seqinfo(f) <- seqinfo(r) <- seq_info
    covRle(f, r)
  }
  if (!is.null(sums)) {
    covs <- lapply(seq_along(sizes), function(k) build(sums[[2 * k - 1]], sums[[2 * k]]))
    merged <- build(sums[[2 * length(sizes) + 1]], sums[[2 * length(sizes) + 2]])
  } else {   # rows not in chromosome/position order
    covs <- lapply(seq_along(sizes), function(k) {
      rows <- which(size_index == k)
      .covRle_from_points(chr[rows], start[rows], strand[rows], score[rows], strand_lv, seq_info, integer_weights)
    })
    merged <- .covRle_from_points(chr, start, strand, score, strand_lv, seq_info, integer_weights)
  }
  sums <- NULL
  save_RDSQS(covRleList(covs, fraction = sizes),
             file = paste0(gsub(paste0("\\.cov", format, "$"), "", out_file, ignore.case = TRUE), ".cov", format))
  if (!is.null(out_file_merged)) {
    seqinfo(merged) <- seq_info   # as export.cov
    save_RDSQS(merged, file = paste0(gsub(paste0("\\.cov", format, "$"), "", out_file_merged, ignore.case = TRUE),
                                     ".cov", format))
  }
  TRUE
}

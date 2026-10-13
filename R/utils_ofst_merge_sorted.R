# Sorted, chromosome-blocked pooled merge (keep_all_scores = FALSE).
#
# Merged/collapsed OFST files keep each chromosome as one contiguous row block (in label or
# factor-code order). The pooled merge is then done per unit of consecutive
# chromosomes: read that row block from every input, collapse with a sorted (keyed) group-by,
# append. This replaces "bind everything, hash group-by, then sort" with bounded,
# cache-friendly units, needs no scratch spill, and its output is identical to the lossless
# ofst_merge() result. Inputs above 2^31 rows in total are fine as long as every unit and the
# output fit (the output size is checked against filter_target_rows; above it the caller falls
# back to the filtering path).

# Per-file chromosome blocks, or NULL if a chromosome is split over several row blocks.
.ofst_chr_blocks <- function(path) {
  s <- fst::read_fst(path, "seqnames")[[1L]]
  if (is.character(s)) s <- factor(s)
  if (!is.factor(s)) return(NULL)
  lv <- levels(s)
  if (!length(s)) return(data.table::data.table(chr = character(), from = numeric(), to = numeric()))
  runs <- ofst_code_runs(s)                       # one C++ pass over the codes, no copy
  v <- runs$value
  if (anyNA(v) || anyDuplicated(v) || any(v < 1L | v > length(lv))) return(NULL)
  data.table::data.table(chr = lv[v], from = c(1, head(runs$end, -1) + 1), to = runs$end)
}

# Rows of one or more [from, to] ranges, concatenated in range order (one multi-range read with
# the fst fork, else one read per range).
.ofst_read_ranges <- function(path, from, to, columns) {
  if (length(from) == 1L) return(.ofst_read_part(path, from = from, to = to, columns = columns))
  if (.fst_multi_range_available() && max(to) < .int32_max)
    return(tryCatch(fst::read_fst(path, columns, from = from, to = to, as.data.table = TRUE),
                    error = function(e) .ofst_abort("cannot read '", path, "': ", conditionMessage(e))))
  data.table::rbindlist(lapply(seq_along(from), function(i)
    .ofst_read_part(path, from = from[i], to = to[i], columns = columns)))
}

# Group consecutive chromosomes (global label order) into units of about unit_rows input rows.
.ofst_chr_units <- function(blocks, unit_rows) {
  totals <- blocks[, .(rows = sum(to - from + 1)), by = chr]
  totals <- totals[order(chr, method = "radix")]
  unit <- integer(nrow(totals)); u <- 1L; acc <- 0
  for (i in seq_len(nrow(totals))) {
    if (acc > 0 && acc + totals$rows[i] > unit_rows) { u <- u + 1L; acc <- 0 }
    unit[i] <- u; acc <- acc + totals$rows[i]
  }
  totals[, unit := unit][]
}

# Strand codes in label order ("*" < "+" < "-"), the order .ofst_sort produces.
.ofst_strand_levels <- c("+", "-", "*")
.ofst_strand_to_sortkey <- c(2L, 3L, 1L)
.ofst_sortkey_to_strand <- c(3L, 1L, 2L)

.ofst_merge_sorted <- function(file_paths, keys, filter_target_rows = 2^31 - 2,
                               unit_rows = 3e7, workers = .ofst_merge_workers(), verbose = TRUE) {
  stopifnot(length(file_paths) > 0L, keys[1L] == "seqnames")
  t0 <- Sys.time()
  blocks <- lapply(file_paths, .ofst_chr_blocks)
  if (any(vapply(blocks, is.null, logical(1)))) return(NULL)   # a chromosome split in a file: fall back
  for (i in seq_along(blocks)) data.table::set(blocks[[i]], j = "file", value = rep.int(i, nrow(blocks[[i]])))
  blocks <- data.table::rbindlist(blocks)
  if (!nrow(blocks)) return(NULL)   # only empty inputs: the general path handles them
  units <- .ofst_chr_units(blocks, unit_rows)
  unit_rows_total <- units[, .(rows = sum(rows)), by = unit]
  if (max(unit_rows_total$rows) >= .int32_max) return(NULL)
  chrs <- units$chr
  schemas <- lapply(file_paths, function(p) fst::metadata_fst(p)$columnNames)
  value_keys <- setdiff(keys, "seqnames")
  n_units <- nrow(unit_rows_total)
  workers <- max(1L, min(as.integer(workers), n_units))
  if (verbose) message("OFST sorted merge: ", length(file_paths), " inputs, ", .ofst_number(sum(units$rows)),
                       " rows, ", length(chrs), " chromosomes in ", n_units, " units, ", workers, " worker(s)")
  merge_unit <- function(u) {
    unit_chrs <- units[unit == u]$chr
    runs <- blocks[chr %in% unit_chrs][order(file, from)]
    sub <- data.table::copy(runs)
    # adjacent chromosome blocks of a file become one range
    sub[, grp := cumsum(c(TRUE, from[-1L] != to[-.N] + 1)), by = file]
    sub <- sub[, .(from = min(from), to = max(to)), by = .(file, grp)]
    sub <- sub[, .(from = list(from), to = list(to)), by = file]
    d <- data.table::rbindlist(lapply(seq_len(nrow(sub)), function(j) {
      f <- sub$file[j]
      has_width <- "width" %in% schemas[[f]] && !"width" %in% keys
      x <- .ofst_read_ranges(file_paths[f], sub$from[[j]], sub$to[[j]],
                             c(value_keys, "score", if (has_width) "width"))
      if (has_width) {
        if (anyNA(x$width) || any(x$width != 1L))
          .ofst_abort("input '", file_paths[f], "' has a width column with values other than 1; ",
                      "cannot merge it with inputs without width.")
        x[, width := NULL]
      }
      if ("strand" %in% value_keys) {
        s <- x$strand
        if (is.character(s)) s <- factor(s, levels = .ofst_strand_levels)
        lv <- levels(s)
        if (any(!lv %in% .ofst_strand_levels))
          .ofst_abort("OFST strand values must be '+', '-', '*' or missing.")
        data.table::set(x, j = "strand", value = .ofst_strand_to_sortkey[match(lv, .ofst_strand_levels)][as.integer(s)])
      }
      for (col in intersect(c("cigar", "cigar1", "cigar2"), names(x)))
        if (is.factor(x[[col]])) data.table::set(x, j = col, value = as.character(x[[col]]))
      x
    }), use.names = TRUE)
    # chromosome from the (file, chromosome) block layout; no seqnames column read
    data.table::set(d, j = "seqnames", value = rep.int(match(runs$chr, chrs), runs$to - runs$from + 1))
    g <- d[, .(score = sum(score, na.rm = TRUE)), keyby = c("seqnames", value_keys)]
    data.table::setattr(g, "sorted", NULL)
    g
  }
  # biggest units first, so parallel workers finish together
  order_units <- unit_rows_total[order(-rows)]$unit
  parts <- if (workers > 1L) {
    parallel::mclapply(order_units, function(u) tryCatch(merge_unit(u), error = function(e) e),
                       mc.cores = workers, mc.preschedule = FALSE)
  } else lapply(order_units, merge_unit)
  failed <- vapply(parts, function(x) inherits(x, "error") || is.null(x), logical(1))
  if (any(failed)) {
    err <- parts[[which(failed)[1L]]]
    .ofst_abort("sorted merge unit failed: ", if (inherits(err, "error")) conditionMessage(err) else
                "worker returned no result (out of memory?)")
  }
  parts <- parts[order(order_units)]
  out_rows <- sum(vapply(parts, nrow, numeric(1)))
  if (out_rows > filter_target_rows) {
    if (verbose) message("OFST sorted merge: output (", .ofst_number(out_rows), " rows) exceeds filter_target_rows; using the filtering path.")
    return(NULL)
  }
  dt <- data.table::rbindlist(parts, use.names = TRUE)
  parts <- NULL
  data.table::set(dt, j = "seqnames", value = structure(dt$seqnames, levels = chrs, class = "factor"))
  if ("strand" %in% value_keys)
    data.table::set(dt, j = "strand", value = structure(.ofst_sortkey_to_strand[dt$strand],
                                                       levels = .ofst_strand_levels, class = "factor"))
  data.table::setcolorder(dt, c(keys, "score"))
  .downcast_score_if_safe(dt)
  if (verbose) message("OFST sorted merge: ", .ofst_number(nrow(dt)), " rows in ",
                       round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), " min")
  dt[]
}

# Parallel units for the sorted merge: forked workers run single-threaded data.table/fst, so use
# about half the data.table threads (option ORFik.ofst_merge_workers overrides). Each worker holds
# one unit (~3e7 input rows, ~1-2 GB).
.ofst_merge_workers <- function() {
  w <- getOption("ORFik.ofst_merge_workers", NULL)
  if (!is.null(w)) return(as.integer(w))
  if (.Platform$OS.type != "unix") return(1L)
  max(1L, min(16L, data.table::getDTthreads() %/% 2L))
}

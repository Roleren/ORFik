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

# Incremental use: `base` is an existing merged (sorted, unique-key) OFST file, e.g. the
# current master. Its rows are not re-grouped: per unit the inputs are grouped as usual and
# the result is merged into the base block with a linear C++ pass (ofst_sorted_sum_merge).
# `subtract` (logical per input) subtracts that input's scores instead (removal of an input
# that is contained in the base); rows whose total becomes 0 are dropped, and a negative
# total is an error.
.ofst_merge_sorted <- function(file_paths, keys, filter_target_rows = 2^31 - 2,
                               unit_rows = 3e7, workers = .ofst_merge_workers(), verbose = TRUE,
                               base = NULL, subtract = rep(FALSE, length(file_paths))) {
  stopifnot(length(file_paths) > 0L, keys[1L] == "seqnames",
            is.logical(subtract), length(subtract) == length(file_paths), !anyNA(subtract))
  t0 <- Sys.time()
  incremental <- !is.null(base) || any(subtract)
  blocks <- lapply(file_paths, .ofst_chr_blocks)
  if (any(vapply(blocks, is.null, logical(1)))) return(NULL)   # a chromosome split in a file: fall back
  for (i in seq_along(blocks)) data.table::set(blocks[[i]], j = "file", value = rep.int(i, nrow(blocks[[i]])))
  blocks <- data.table::rbindlist(blocks)
  base_blocks <- NULL
  if (!is.null(base)) {
    base_blocks <- .ofst_chr_blocks(base)
    if (is.null(base_blocks)) return(NULL)
    if (anyDuplicated(base_blocks$chr)) return(NULL)
    data.table::set(base_blocks, j = "file", value = rep.int(0L, nrow(base_blocks)))
  }
  all_blocks <- data.table::rbindlist(list(blocks, base_blocks))
  if (!nrow(all_blocks)) return(NULL)   # only empty inputs: the general path handles them
  units <- .ofst_chr_units(all_blocks, unit_rows)
  unit_rows_total <- units[, .(rows = sum(rows)), by = unit]
  if (max(unit_rows_total$rows) >= .int32_max) return(NULL)
  chrs <- units$chr
  schemas <- lapply(file_paths, function(p) fst::metadata_fst(p)$columnNames)
  base_schema <- if (!is.null(base)) fst::metadata_fst(base)$columnNames
  value_keys <- setdiff(keys, "seqnames")
  n_units <- nrow(unit_rows_total)
  workers <- max(1L, min(as.integer(workers), n_units))
  if (verbose) message("OFST sorted merge: ", length(file_paths), " inputs",
                       if (!is.null(base)) " + base" else "", if (any(subtract)) paste0(" (", sum(subtract), " subtracted)") else "",
                       ", ", .ofst_number(sum(units$rows)),
                       " rows, ", length(chrs), " chromosomes in ", n_units, " units, ", workers, " worker(s)")
  read_input <- function(path, schema, from, to) {
    has_width <- "width" %in% schema && !"width" %in% keys
    x <- .ofst_read_ranges(path, from, to, c(value_keys, "score", if (has_width) "width"))
    if (has_width) {
      if (anyNA(x$width) || any(x$width != 1L))
        .ofst_abort("input '", path, "' has a width column with values other than 1; ",
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
  }
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
      x <- read_input(file_paths[f], schemas[[f]], sub$from[[j]], sub$to[[j]])
      if (subtract[f]) data.table::set(x, j = "score", value = -as.numeric(x$score))
      x
    }), use.names = TRUE)
    # chromosome from the (file, chromosome) block layout; no seqnames column read
    if (nrow(d)) data.table::set(d, j = "seqnames", value = rep.int(match(runs$chr, chrs), runs$to - runs$from + 1))
    if (!nrow(d)) {
      g <- data.table::data.table(seqnames = integer())
      for (k in value_keys) data.table::set(g, j = k, value = integer())
      data.table::set(g, j = "score", value = numeric())
    } else g <- d[, .(score = sum(score, na.rm = TRUE)), keyby = c("seqnames", value_keys)]
    d <- NULL
    data.table::setattr(g, "sorted", NULL)
    if (!is.null(base)) {
      # base block of this unit, in chromosome code order (the base is sorted within each block)
      bb <- base_blocks[chr %in% unit_chrs]
      bb <- bb[order(match(chr, chrs))]
      b <- data.table::rbindlist(lapply(seq_len(nrow(bb)), function(j) {
        x <- read_input(base, base_schema, bb$from[j], bb$to[j])
        data.table::set(x, j = "seqnames", value = rep.int(match(bb$chr[j], chrs), nrow(x)))
        x
      }), use.names = TRUE)
      if (nrow(b)) {
        kc <- c("seqnames", value_keys)
        int_keys <- all(vapply(kc, function(k) is.integer(b[[k]]) && is.integer(g[[k]]), logical(1)))
        # linear merge; a base that is not in merge order (not written by ofst_merge) is grouped instead
        r <- if (int_keys) tryCatch(
          ofst_sorted_sum_merge(.ofst_cols(b, kc), as.numeric(b$score), .ofst_cols(g, kc), as.numeric(g$score),
                                drop_zero = any(subtract), forbid_negative = any(subtract)),
          error = function(e) if (grepl("first input is not sorted", conditionMessage(e))) NULL else stop(e))
        if (!is.null(r)) {
          b <- NULL
          g <- data.table::setDT(r$keys)
          data.table::set(g, j = "score", value = r$score)
        } else {   # non-integer keys (e.g. cigar strings) or unsorted base: group base and delta together
          g <- data.table::rbindlist(list(b, g), use.names = TRUE)[
            , .(score = sum(score, na.rm = TRUE)), keyby = kc]
          data.table::setattr(g, "sorted", NULL)
        }
      }
    }
    if (any(subtract)) {
      if (any(g$score < 0)) .ofst_abort("a merged score is negative: a subtracted input is not contained in the base.")
      g <- g[score != 0]
    }
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
    if (incremental) .ofst_abort("incremental merge output (", .ofst_number(out_rows),
                                 " rows) exceeds filter_target_rows; the incremental path never filters.")
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

# Global OFST rescue. Scratch partitions contain complete chromosome/start
# groups; no filtering decision is made on a subset of the input libraries.
.ofst_abort <- function(...) stop("OFST merge: ", ..., call. = FALSE)
.ofst_number <- function(x) format(x, scientific = FALSE, trim = TRUE, big.mark = ",")

.ofst_rng_restore <- function() {
  # FST's Rcpp RNGScope can create .Random.seed even though our hash-based
  # selection never samples from R's RNG. Preserve an uninitialized session too.
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  saved <- if (existed) get(".Random.seed", envir = .GlobalEnv) + 0L else NULL
  function() {
    if (existed) assign(".Random.seed", saved, envir = .GlobalEnv) else
      if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
        rm(".Random.seed", envir = .GlobalEnv)
    invisible(NULL)
  }
}

.ofst_save_merge <- function(dt, path) {
  # Prepare both artifacts before replacing either. A NULL sidecar explicitly
  # records an unfiltered merge, including when overwriting a filtered output.
  output_tmp <- tempfile(".ofst-output-", tmpdir = dirname(path))
  summary_tmp <- tempfile(".ofst-summary-", tmpdir = dirname(path))
  preserve_summary <- FALSE
  on.exit({
    unlink(output_tmp)
    if (!preserve_summary) unlink(summary_tmp)
  }, add = TRUE)
  tryCatch({
    fst::write_fst(dt, output_tmp)
    saveRDS(attr(dt, "removal_summary", exact = TRUE), summary_tmp)
  }, error = function(e) .ofst_abort("cannot prepare merged output and summary for '", path,
                                    "'. Check disk space/permissions; existing output was not replaced. ", conditionMessage(e)))
  if (!file.rename(output_tmp, path)) .ofst_abort("cannot replace merged output '", path, "'.")
  sidecar <- paste0(path, ".removal_summary.rds")
  if (!file.rename(summary_tmp, sidecar)) {
    preserve_summary <- TRUE
    .ofst_abort("merged output was written to '", path, "', but the summary could not replace '", sidecar,
                "'. The correct summary is preserved at '", summary_tmp, "'; do not use any old sidecar.")
  }
  message("OFST merge: saved ", path, "; removal summary: ", sidecar,
          if (is.null(attr(dt, "removal_summary", exact = TRUE))) " (NULL: no filtering)." else ".")
  invisible(NULL)
}

.ofst_filter_controls <- function(allow_filtering, filter_target_rows, filter_seed,
                                  max_filter_score, filter_chunk_rows, filter_tmpdir) {
  if (!is.logical(allow_filtering) || !is.null(dim(allow_filtering)) || length(allow_filtering) != 1L || is.na(allow_filtering))
    .ofst_abort("allow_filtering must be TRUE or FALSE.")
  whole <- function(x) is.numeric(x) && !is.object(x) && is.null(dim(x)) && length(x) == 1L &&
    !is.na(x) && is.finite(x) && x == trunc(x)
  if (!whole(filter_target_rows) || filter_target_rows < 0 || filter_target_rows >= .int32_max)
    .ofst_abort("filter_target_rows must be a whole number from 0 to ", .ofst_number(.int32_max - 1),
                "; it controls final output rows, not chunk size.")
  if (!whole(filter_seed) || filter_seed < 0 || filter_seed > .Machine$integer.max)
    .ofst_abort("filter_seed must be a whole number from 0 to 2147483647.")
  if (!is.numeric(max_filter_score) || is.object(max_filter_score) || !is.null(dim(max_filter_score)) || length(max_filter_score) != 1L ||
      is.na(max_filter_score) || max_filter_score < 0)
    .ofst_abort("max_filter_score must be non-negative (Inf permits any pooled score).")
  if (!whole(filter_chunk_rows) || filter_chunk_rows < 1 || filter_chunk_rows > floor(.int32_max / 4))
    .ofst_abort("filter_chunk_rows must be a positive whole number no larger than ",
                .ofst_number(floor(.int32_max / 4)), ".")
  if (!is.character(filter_tmpdir) || !is.null(dim(filter_tmpdir)) || length(filter_tmpdir) != 1L || is.na(filter_tmpdir) || !nzchar(filter_tmpdir))
    .ofst_abort("filter_tmpdir must name an existing writable scratch directory.")
  invisible(NULL)
}

.ofst_scratch_dir <- function(parent) {
  if (!dir.exists(parent) || file.access(parent, 2L) != 0L)
    .ofst_abort("scratch directory '", parent, "' does not exist or is not writable. Set filter_tmpdir to a disk with free space.")
  path <- tempfile("orfik-ofst-filter-", tmpdir = normalizePath(parent, mustWork = TRUE))
  if (!dir.create(path)) .ofst_abort("could not create scratch directory '", path, "'.")
  path
}

.ofst_spill <- function(dt, scratch, prefix = "part-") {
  path <- tempfile(prefix, tmpdir = scratch, fileext = ".fst")
  tryCatch(fst::write_fst(dt, path), error = function(e)
    .ofst_abort("cannot write scratch file '", path, "'. Check disk space and permissions. ", conditionMessage(e)))
  path
}

.ofst_read_large_offset <- function(path, from, to, columns = NULL) {
  # fst 0.9.8's public read_fst coerces offsets with as.integer(), although
  # the native reader supports 64-bit offsets. Keep this compatibility bridge
  # isolated; only the requested short block is converted into a data.table.
  if (is.null(to) || !is.finite(from) || !is.finite(to) || from < 1 ||
      to < from || to - from + 1 >= .int32_max)
    .ofst_abort("large-offset reads require an explicit bounded from/to range below the data.table row ceiling.")
  reader <- utils::getFromNamespace("fstretrieve", "fst")
  expected <- c("fileName", "columnSelection", "startRow", "endRow")
  if (!identical(names(formals(reader)), expected))
    .ofst_abort("this fst version has an unsupported native reader interface for offsets above 2^31-1. Use input files split below that size.")
  raw <- reader(normalizePath(path, mustWork = TRUE), columns, as.double(from), as.double(to))
  if (inherits(raw, "fst_error")) stop(raw)
  if (!is.list(raw$resTable)) .ofst_abort("fst returned an unexpected large-offset block format.")
  result <- data.table::as.data.table(raw$resTable)
  if (length(raw$keyNames)) data.table::setattr(result, "sorted", raw$keyNames)
  result
}

.ofst_read_part <- function(path, from = 1, to = NULL, columns = NULL) {
  tryCatch({
    if (from > .Machine$integer.max || (!is.null(to) && to > .Machine$integer.max))
      .ofst_read_large_offset(path, from, to, columns) else
        fst::read_fst(path, as.data.table = TRUE, from = from, to = to, columns = columns)
  }, error = function(e)
    .ofst_abort("cannot read '", path, "': ", conditionMessage(e)))
}

.ofst_position_id <- function(seqnames, start) {
  chromosome <- enc2utf8(as.character(seqnames))
  paste0(nchar(chromosome, type = "bytes"), ":", chromosome, ":", start)
}

.ofst_hash <- function(id, seed, salt) {
  x <- as.double(digest::digest2int(paste0(salt, id), seed = as.integer(seed)))
  # INT_MIN has R's NA bit pattern; no missing strings are passed here.
  x[is.na(x)] <- -2147483648
  x + 2147483648
}

.ofst_priority <- function(seqnames, start, seed) {
  id <- .ofst_position_id(seqnames, start)
  # An exactly representable 53-bit integer priority. Salts separate the two
  # hashes and the partition hash. No calls to R's global RNG are made.
  .ofst_hash(id, seed, "rank-a:") * 2^21 + .ofst_hash(id, seed, "rank-b:") %% 2^21
}

.ofst_partition_id <- function(dt, depth, bins = 2L) {
  .ofst_hash(.ofst_position_id(dt$seqnames, dt$start), 0L, "partition:") %/% 2^depth %% bins
}

.ofst_reduce_scores <- function(dt, keys) {
  dt[, .(score = sum(score, na.rm = TRUE)), by = c(keys)]
}

.ofst_histogram <- function(positions, column) {
  positions[, .(weight = sum(as.double(alignment_rows))), by = .(value = get(column))]
}

.ofst_finish_partition <- function(dt, keys, source_col, scratch) {
  pooled <- .ofst_reduce_scores(dt, keys)
  positions <- pooled[, .(position_score = sum(score), alignment_rows = .N), by = .(seqnames, start)]
  if (any(!is.finite(positions$position_score)))
    .ofst_abort("pooled position scores overflowed; finite non-negative scores are required for filtering.")
  list(data = .ofst_spill(dt, scratch, "data-"),
       positions = .ofst_spill(positions, scratch, "positions-"),
       histogram = .ofst_spill(.ofst_histogram(positions, "position_score"), scratch, "scores-"),
       rows = as.double(nrow(pooled)), positions_n = as.double(nrow(positions)),
       score = sum(positions$position_score))
}

.ofst_compact_partition <- function(paths, keys, source_col, chunk_rows, scratch, depth) {
  accumulator <- NULL
  for (i in seq_along(paths)) {
    incoming <- .ofst_read_part(paths[i])
    combined <- if (is.null(accumulator)) incoming else
      data.table::rbindlist(list(accumulator, incoming), use.names = TRUE)
    accumulator <- incoming <- NULL
    accumulator <- .ofst_reduce_scores(combined, c(keys, source_col))
    combined <- NULL
    unlink(paths[i]) # Only this helper's owned scratch fragments.
    if (nrow(accumulator) > chunk_rows) {
      if (depth >= 32L)
        .ofst_abort("a chromosome/start group or partition needs more than filter_chunk_rows=",
                    .ofst_number(chunk_rows), " rows even after duplicate collapse. Increase filter_chunk_rows if RAM permits.")
      remaining <- c(.ofst_spill(accumulator, scratch), paths[seq_along(paths) > i])
      accumulator <- NULL
      children <- list(character(), character())
      for (path in remaining) {
        dt <- .ofst_read_part(path)
        bucket <- .ofst_partition_id(dt, depth)
        for (j in sort(unique(bucket))) {
          index <- which(bucket == j)
          children[[j + 1L]] <- c(children[[j + 1L]], .ofst_spill(dt[index], scratch))
        }
        dt <- NULL
        unlink(path)
      }
      if (chunk_rows >= 1e5) gc()
      answer <- list()
      for (child in children) if (length(child))
        answer <- c(answer, .ofst_compact_partition(child, keys, source_col, chunk_rows, scratch, depth + 1L))
      return(answer)
    }
  }
  if (is.null(accumulator) || !nrow(accumulator)) return(list())
  list(.ofst_finish_partition(accumulator, keys, source_col, scratch))
}

.ofst_stage_inputs <- function(paths, row_counts, keys, source_col, chunk_rows, scratch) {
  bins <- 2^ceiling(log2(min(64, max(1, sum(row_counts) / chunk_rows))))
  buckets <- rep(list(character()), bins)
  template <- NULL
  for (i in seq_along(paths)) {
    message("OFST rescue: reading input ", i, "/", length(paths), " (",
            .ofst_number(row_counts[i]), " rows): ", paths[i])
    if (!row_counts[i]) next
    first <- 1
    while (first <= row_counts[i]) {
      last <- min(row_counts[i], first + chunk_rows - 1)
      dt <- tryCatch(.normalize_dt(.ofst_read_part(paths[i], from = first, to = last)),
                     error = function(e) .ofst_abort("invalid input '", paths[i], "', rows ", first, "-", last, ": ", conditionMessage(e)))
      if (!"start" %in% names(dt))
        .ofst_abort("chromosome/start filtering requires a 'start' column. Paired start1/start2 inputs need an explicit positional convention first.")
      if (anyNA(dt$seqnames) || anyNA(dt$start))
        .ofst_abort("input '", paths[i], "' contains missing chromosome/start keys. Correct these before positional filtering.")
      if (any(!is.finite(dt$score) & !is.na(dt$score)) || any(dt$score < 0, na.rm = TRUE))
        .ofst_abort("input '", paths[i], "' has negative or infinite scores. Filtering requires non-negative finite scores; missing scores contribute zero.")
      data.table::set(dt, j = "score", value = as.double(dt$score))
      dt[is.na(score), score := 0]
      # Drop only keys explicitly excluded by keepCigar, retaining all others.
      dt <- dt[, c(keys, "score"), with = FALSE]
      if (is.null(template)) template <- dt[0]
      data.table::set(dt, j = source_col, value = rep.int(i, nrow(dt)))
      dt <- .ofst_reduce_scores(dt, c(keys, source_col))
      bucket <- .ofst_partition_id(dt, 0L, bins)
      for (j in sort(unique(bucket))) {
        index <- which(bucket == j)
        buckets[[j + 1L]] <- c(buckets[[j + 1L]], .ofst_spill(dt[index], scratch))
      }
      dt <- NULL
      first <- last + 1
    }
    if (chunk_rows >= 1e5) gc()
  }
  leaves <- list()
  for (i in seq_along(buckets)) if (length(buckets[[i]])) {
    message("OFST rescue: compacting positional partition ", i, "/", bins, ".")
    leaves <- c(leaves, .ofst_compact_partition(buckets[[i]], keys, source_col,
                                              chunk_rows, scratch, as.integer(log2(bins))))
    buckets[[i]] <- character()
    if (chunk_rows >= 1e5) gc()
  }
  list(leaves = leaves, template = template)
}

.ofst_scan_histograms <- function(paths, bound, ceiling = Inf) {
  answer <- list(eligible = 0, less = 0, le = 0, min = Inf, max = -Inf,
                 max_le = -Inf, min_above = Inf)
  for (path in paths) {
    h <- .ofst_read_part(path)
    h <- h[value <= ceiling]
    if (!nrow(h)) next
    below <- h$value < bound
    le <- h$value <= bound
    answer$eligible <- answer$eligible + sum(h$weight)
    answer$less <- answer$less + sum(h$weight[below])
    answer$le <- answer$le + sum(h$weight[le])
    answer$min <- min(answer$min, min(h$value))
    answer$max <- max(answer$max, max(h$value))
    if (any(le)) answer$max_le <- max(answer$max_le, max(h$value[le]))
    if (any(!le)) answer$min_above <- min(answer$min_above, min(h$value[!le]))
  }
  answer
}

.ofst_weighted_cutoff <- function(paths, needed, ceiling = Inf) {
  initial <- .ofst_scan_histograms(paths, Inf, ceiling)
  if (initial$eligible < needed)
    .ofst_abort("must remove at least ", .ofst_number(needed), " alignment rows, but max_filter_score=",
                ceiling, " permits only ", .ofst_number(initial$eligible),
                ". Raise max_filter_score or choose a different output strategy; no filtered result was returned.")
  low <- initial$min
  high <- initial$max
  iterations <- 0L
  while (low < high) {
    iterations <- iterations + 1L
    if (iterations > 128L) .ofst_abort("could not resolve a filtering cutoff after 128 bounded passes; inspect score range/precision.")
    middle <- low / 2 + high / 2
    if (middle >= high) middle <- low
    scan <- .ofst_scan_histograms(paths, middle, ceiling)
    if (scan$le >= needed) high <- scan$max_le else low <- scan$min_above
  }
  list(value = low, below = .ofst_scan_histograms(paths, low, ceiling)$less)
}

.ofst_priority_cutoff <- function(paths, needed, scratch, chunk_rows) {
  # Radix selection on exact 53-bit priorities. Each pass keeps only one of
  # 1024 buckets on disk, rather than rescanning every tied position 53 times.
  below <- lower <- 0
  for (shift in c(43, 33, 23, 13, 3, 0)) {
    weights <- numeric(1024L)
    rows <- 0
    for (path in paths) {
      h <- .ofst_read_part(path)
      rows <- rows + nrow(h)
      h[, bucket := as.integer((value - lower) %/% 2^shift) + 1L]
      counts <- h[, .(weight = sum(weight)), by = bucket]
      weights[counts$bucket] <- weights[counts$bucket] + counts$weight
    }
    if (rows <= chunk_rows) {
      h <- data.table::rbindlist(lapply(paths, .ofst_read_part))
      h <- h[, .(weight = sum(weight)), by = value]
      data.table::setorderv(h, "value")
      index <- which(cumsum(h$weight) >= needed)[1L]
      if (is.na(index)) .ofst_abort("internal priority selection has insufficient eligible rows.")
      return(list(value = h$value[index], below = below + sum(h$weight[seq_len(index - 1L)])))
    }
    selected <- which(cumsum(weights) >= needed)[1L]
    if (is.na(selected)) .ofst_abort("internal priority buckets do not cover the requested removal.")
    prior <- sum(weights[seq_len(selected - 1L)])
    below <- below + prior
    needed <- needed - prior
    new_lower <- lower + (selected - 1L) * 2^shift
    if (shift == 0) return(list(value = new_lower, below = below))
    narrowed <- character()
    for (path in paths) {
      h <- .ofst_read_part(path)
      h <- h[value >= new_lower & value < new_lower + 2^shift]
      if (nrow(h)) narrowed <- c(narrowed, .ofst_spill(h, scratch, "priority-window-"))
      unlink(path) # Only priority histograms created by this rescue.
    }
    paths <- narrowed
    lower <- new_lower
  }
  .ofst_abort("internal priority selection exhausted its 53-bit range.")
}

.ofst_choose_removals <- function(leaves, needed, seed, max_score, scratch, chunk_rows) {
  histograms <- vapply(leaves, `[[`, "", "histogram")
  score_cut <- .ofst_weighted_cutoff(histograms, needed, max_score)
  message("OFST filtering: lowest pooled scores first; cutoff score ", score_cut$value,
          ". Randomizing tied positions globally with seed ", seed, ".")
  ties <- histograms <- character()
  for (leaf in leaves) {
    p <- .ofst_read_part(leaf$positions)
    p <- p[position_score == score_cut$value]
    if (!nrow(p)) next
    p[, priority := .ofst_priority(seqnames, start, seed)]
    ties <- c(ties, .ofst_spill(p, scratch, "ties-"))
    histograms <- c(histograms, .ofst_spill(.ofst_histogram(p, "priority"), scratch, "priorities-"))
  }
  priority_cut <- .ofst_priority_cutoff(histograms, needed - score_cut$below, scratch, chunk_rows)
  boundary <- list()
  n_boundary <- 0
  for (path in ties) {
    p <- .ofst_read_part(path)
    p <- p[priority == priority_cut$value]
    n_boundary <- n_boundary + nrow(p)
    if (n_boundary > chunk_rows)
      .ofst_abort("random-priority collision group exceeds filter_chunk_rows. Try a different filter_seed or increase the chunk budget.")
    if (nrow(p)) boundary[[length(boundary) + 1L]] <- p
  }
  boundary <- data.table::rbindlist(boundary)
  boundary[, seqnames := as.character(seqnames)]
  data.table::setorderv(boundary, c("seqnames", "start"))
  missing_rows <- needed - score_cut$below - priority_cut$below
  last <- which(cumsum(as.double(boundary$alignment_rows)) >= missing_rows)[1L]
  if (is.na(last)) .ofst_abort("internal filtering selection did not cover the required row reduction.")
  selected <- boundary[seq_len(last)]
  list(score = score_cut$value, priority = priority_cut$value,
       boundary_ids = .ofst_position_id(selected$seqnames, selected$start), seed = seed)
}

.ofst_removal_mask <- function(p, selection) {
  removed <- p$position_score < selection$score
  tied <- which(p$position_score == selection$score)
  if (length(tied)) {
    priority <- .ofst_priority(p$seqnames[tied], p$start[tied], selection$seed)
    boundary <- priority == selection$priority
    selected <- priority < selection$priority
    if (any(boundary)) selected[boundary] <-
      .ofst_position_id(p$seqnames[tied[boundary]], p$start[tied[boundary]]) %in% selection$boundary_ids
    removed[tied] <- selected
  }
  removed
}

.ofst_chromosome_removal_stats <- function(p, removed) {
  p[, .(positions_before = as.double(.N), positions_removed = sum(removed[.I]),
         rows_before = sum(as.double(alignment_rows)), rows_removed = sum(as.double(alignment_rows[removed[.I]])),
         score_before = sum(position_score), score_removed = sum(position_score[removed[.I]])), by = seqnames]
}

.ofst_leaf_output <- function(dt, source_col, keys, lib_names, keep_all_scores) {
  if (!keep_all_scores) return(.ofst_reduce_scores(dt, keys))
  indices <- sort(unique(dt[[source_col]]))
  if (!length(indices)) return(dt[0, c(keys, "score"), with = FALSE])
  labels <- unique(lib_names[indices])
  tables <- lapply(labels, function(label) {
    sources <- which(lib_names == label)
    index <- which(dt[[source_col]] %in% sources)
    d <- dt[index, c(keys, "score"), with = FALSE]
    d <- .ofst_reduce_scores(d, keys)
    data.table::setnames(d, "score", label)
    d
  })
  result <- .merge_by_keys_reduce(tables, keys)
  .recompute_total_score(result, labels)
}

.ofst_merge_filtered <- function(file_paths, lib_names, row_counts, keys,
                                  keep_all_scores, controls) {
  if (!all(c("seqnames", "start") %in% keys))
    .ofst_abort("positional rescue requires chromosome and 'start' columns; paired start1/start2 coordinates are not implicitly substituted.")
  scratch <- .ofst_scratch_dir(controls$filter_tmpdir)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  message("OFST rescue: ", .ofst_number(sum(row_counts)), " input rows may exceed the final target of ",
          .ofst_number(controls$filter_target_rows), ". Checking the distinct merged rows before removing anything.")
  message("OFST rescue: scratch directory ", scratch, "; block budget ", .ofst_number(controls$filter_chunk_rows), " rows.")
  source_col <- ".ofst_source"
  while (source_col %in% c(keys, lib_names, "score")) source_col <- paste0("_", source_col)
  staged <- .ofst_stage_inputs(file_paths, row_counts, keys, source_col, controls$filter_chunk_rows, scratch)
  leaves <- staged$leaves
  before <- sum(vapply(leaves, `[[`, 0, "rows"))
  positions_before <- sum(vapply(leaves, `[[`, 0, "positions_n"))
  score_before <- sum(vapply(leaves, `[[`, 0, "score"))
  if (!is.finite(score_before)) .ofst_abort("the aggregate input score overflowed numeric precision; finite score totals are required for a meaningful removal summary.")
  message("OFST rescue: distinct merged rows = ", .ofst_number(before), "; chromosome/start positions = ",
          .ofst_number(positions_before), ".")
  selection <- NULL
  if (before > controls$filter_target_rows) {
    if (!controls$allow_filtering)
      .ofst_abort(.ofst_number(before), " distinct merged rows exceed filter_target_rows=",
                  .ofst_number(controls$filter_target_rows), ", but allow_filtering=FALSE. Enable filtering or retain partitioned output.")
    selection <- .ofst_choose_removals(leaves, before - controls$filter_target_rows,
                                       controls$filter_seed, controls$max_filter_score,
                                       scratch, controls$filter_chunk_rows)
  } else message("OFST rescue: duplicates collapsed below the target; no filtering is needed.")
  outputs <- chromosome_stats <- input_stats <- list()
  for (i in seq_along(leaves)) {
    p <- .ofst_read_part(leaves[[i]]$positions)
    removed <- if (is.null(selection)) rep(FALSE, nrow(p)) else .ofst_removal_mask(p, selection)
    if (!is.null(selection)) chromosome_stats[[i]] <- .ofst_chromosome_removal_stats(p, removed)
    dt <- .ofst_read_part(leaves[[i]]$data)
    if (!is.null(selection)) {
      by_input <- dt[, .(rows_before = as.double(.N), score_before = sum(score)), by = c(source_col)]
      drop <- p[removed, .(seqnames, start)]
      dt <- dt[!drop, on = c("seqnames", "start")]
      after_input <- dt[, .(rows_after = as.double(.N), score_after = sum(score)), by = c(source_col)]
      stats <- merge(by_input, after_input, by = source_col, all.x = TRUE, sort = FALSE)
      for (col in c("rows_after", "score_after")) data.table::set(stats, which(is.na(stats[[col]])), col, 0)
      input_stats[[i]] <- stats
    }
    outputs[[i]] <- .ofst_leaf_output(dt, source_col, keys, lib_names, keep_all_scores)
    dt <- p <- NULL
  }
  message("OFST rescue: assembling retained alignment rows; the final returned table must still fit in RAM.")
  result <- if (length(outputs)) data.table::rbindlist(outputs, use.names = TRUE, fill = TRUE) else staged$template
  outputs <- NULL
  if (keep_all_scores) {
    for (col in setdiff(lib_names, names(result))) result[, (col) := NA_integer_]
    data.table::setcolorder(result, c(keys, unique(lib_names), "score"))
    .downcast_cols_if_safe(result, lib_names)
  }
  .downcast_score_if_safe(result)
  if (nrow(result) > controls$filter_target_rows)
    .ofst_abort("internal filtering error: retained rows still exceed the target; no result returned.")
  if (!is.null(selection)) {
    by_chromosome <- data.table::rbindlist(chromosome_stats)
    by_chromosome[, seqnames := as.character(seqnames)]
    by_chromosome <- by_chromosome[, lapply(.SD, sum), by = seqnames]
    by_chromosome[, `:=`(rows_after = rows_before - rows_removed,
                         positions_after = positions_before - positions_removed,
                         score_after = score_before - score_removed)]
    data.table::setorderv(by_chromosome, "seqnames")
    by_input <- data.table::rbindlist(input_stats)
    by_input <- by_input[, lapply(.SD, sum), by = c(source_col)]
    data.table::setnames(by_input, source_col, "input_index")
    by_input <- merge(data.table::data.table(input_index = seq_along(file_paths)), by_input,
                      by = "input_index", all.x = TRUE, sort = TRUE)
    for (col in setdiff(names(by_input), "input_index")) data.table::set(by_input, which(is.na(by_input[[col]])), col, 0)
    by_input[, `:=`(file = file_paths[input_index], library = lib_names[input_index],
                    raw_input_rows = row_counts[input_index], rows_removed = rows_before - rows_after,
                    score_removed = score_before - score_after)]
    summary <- structure(list(schema_version = 1L, applied = TRUE, method = "lowest_pooled_position_score",
                              input_kind = "files", alignment_keys = keys,
                              grouping = c("seqnames", "start"), ignores_strand_and_cigar = TRUE,
                              seed = controls$filter_seed, tie_method = "seeded_digest_53bit_priority",
                              target_rows = controls$filter_target_rows,
                              target_undershoot_rows = controls$filter_target_rows - nrow(result),
                              filter_chunk_rows = controls$filter_chunk_rows,
                              max_filter_score = controls$max_filter_score,
                              input_rows = sum(row_counts), rows_before = before,
                              rows_after = as.double(nrow(result)), rows_removed = before - nrow(result),
                              positions_before = positions_before,
                              positions_removed = sum(by_chromosome$positions_removed),
                              positions_after = positions_before - sum(by_chromosome$positions_removed),
                              score_before = score_before, score_after = sum(result$score),
                              score_removed = sum(by_chromosome$score_removed),
                              cutoff_score = selection$score, cutoff_priority = selection$priority,
                              cutoff_collision_positions = selection$boundary_ids,
                              by_chromosome = by_chromosome, by_input = by_input),
                         class = c("ofst_removal_summary", "list"))
    data.table::setattr(result, "removal_summary", summary)
    message("OFST filtering: removed ", .ofst_number(summary$positions_removed), " positions / ",
            .ofst_number(summary$rows_removed), " alignment rows / score ", .ofst_number(summary$score_removed),
            "; retained ", .ofst_number(summary$rows_after), " rows. Summary attached as attr(result, 'removal_summary').")
  }
  result
}

# Kept as an internal entry point for callers that already hold their tables.
# The public file entry point avoids materializing oversized inputs at all.
ofst_merge_internal <- function(dt_list, lib_names, keep_all_scores = TRUE,
                                keepCigar = TRUE, sort = TRUE, chunkified = FALSE,
                                allow_filtering = TRUE, max_filter_value = NULL,
                                filter_target_rows = 2^31 - 2, filter_seed = 1L,
                                max_filter_score = Inf, filter_chunk_rows = 5e6,
                                filter_tmpdir = tempdir()) {
  restore_rng <- .ofst_rng_restore()
  on.exit(restore_rng(), add = TRUE)
  if (!is.null(max_filter_value)) {
    if (!missing(max_filter_score)) .ofst_abort("use max_filter_score, not both max_filter_score and legacy max_filter_value.")
    warning("max_filter_value is deprecated; use max_filter_score (a ceiling on the pooled chromosome/start score).", call. = FALSE)
    max_filter_score <- max_filter_value
  }
  .ofst_filter_controls(allow_filtering, filter_target_rows, filter_seed,
                        max_filter_score, filter_chunk_rows, filter_tmpdir)
  for (arg in c("keep_all_scores", "keepCigar", "sort", "chunkified")) {
    x <- get(arg)
    if (!is.logical(x) || !is.null(dim(x)) || length(x) != 1L || is.na(x)) .ofst_abort(arg, " must be TRUE or FALSE.")
  }
  if (!is.list(dt_list) || !length(dt_list) ||
      !all(vapply(dt_list, data.table::is.data.table, logical(1))))
    .ofst_abort("dt_list must contain at least one data.table, with no non-table elements.")
  if (!is.character(lib_names) || anyNA(lib_names) || any(!nzchar(lib_names)) ||
      !is.null(dim(lib_names)) ||
      (!chunkified && length(lib_names) != length(dt_list)))
    .ofst_abort("lib_names must be non-empty character names, one per unchunked input table.")
  dt_list <- lapply(seq_along(dt_list), function(i) tryCatch(.normalize_dt(dt_list[[i]]),
    error = function(e) .ofst_abort("invalid input table ", i, ": ", conditionMessage(e))))
  keys <- setdiff(names(dt_list[[1L]]), c("score", if (chunkified && keep_all_scores) lib_names))
  for (i in seq_along(dt_list)) {
    input_keys <- setdiff(names(dt_list[[i]]), c("score", if (chunkified && keep_all_scores) lib_names))
    if (!setequal(keys, input_keys)) .ofst_abort("input table ", i, " has different alignment columns; identical keys are required.")
    if (chunkified && keep_all_scores) {
      present <- intersect(lib_names, names(dt_list[[i]]))
      if (!length(present)) .ofst_abort("chunk ", i, " has no per-library score columns.")
      for (label in present) {
        x <- dt_list[[i]][[label]]
        if (!is.numeric(x) || is.object(x) || !is.null(dim(x)))
          .ofst_abort("chunk ", i, ", library '", label, "' must have numeric scores.")
      }
    }
  }
  if (keep_all_scores && (anyDuplicated(lib_names) || any(lib_names %in% c(keys, "score"))))
    .ofst_abort("lib_names must be unique and must not collide with alignment columns or score.")
  rows <- vapply(dt_list, nrow, 0)
  if (sum(rows) <= filter_target_rows)
    return(.ofst_merge_lossless(dt_list, lib_names, keep_all_scores, keepCigar, sort, chunkified))
  if (!keepCigar) keys <- setdiff(keys, "cigar")
  scratch <- .ofst_scratch_dir(filter_tmpdir)
  on.exit(unlink(scratch, recursive = TRUE), add = TRUE)
  paths <- labels <- character()
  for (i in seq_along(dt_list)) {
    d <- dt_list[[i]]
    if (chunkified && keep_all_scores) {
      present <- intersect(lib_names, names(d))
      if (!length(present)) .ofst_abort("chunk ", i, " has no per-library score columns.")
      observed <- logical(nrow(d))
      for (label in present) {
        x <- d[[label]]
        observed <- observed | !is.na(x)
        index <- which(!is.na(x))
        part <- d[index, c(keys, label), with = FALSE]
        data.table::setnames(part, label, "score")
        paths <- c(paths, .ofst_spill(part, scratch, "input-"))
        labels <- c(labels, label)
      }
      if (any(!observed)) .ofst_abort("chunk ", i, " contains rows absent from every per-library score column; cannot reconstruct their library attribution.")
    } else {
      paths <- c(paths, .ofst_spill(d[, c(keys, "score"), with = FALSE], scratch, "input-"))
      labels <- c(labels, if (chunkified) paste0("chunk", i) else lib_names[i])
    }
  }
  controls <- list(allow_filtering = allow_filtering, filter_target_rows = filter_target_rows,
                   filter_seed = filter_seed, max_filter_score = max_filter_score,
                   filter_chunk_rows = filter_chunk_rows, filter_tmpdir = filter_tmpdir)
  result <- .ofst_merge_filtered(paths, labels,
    vapply(paths, function(p) as.double(fst::metadata_fst(p)$nrOfRows), 0), keys, keep_all_scores, controls)
  summary <- attr(result, "removal_summary", exact = TRUE)
  if (!is.null(summary)) {
    summary$input_kind <- if (chunkified && keep_all_scores) "chunk_library_contributions" else "in_memory_tables"
    summary$by_input[, file := NA_character_]
    data.table::setattr(result, "removal_summary", summary)
  }
  if (sort) .ofst_sort(result, setdiff(names(result), "score"))
  .ofst_finalize(result)
}

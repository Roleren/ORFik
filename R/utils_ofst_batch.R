# Bounded batch-first merging. Runs form a binary hash tree whose leaves fit
# the row budget. A binary carry merges equally sized input batches promptly;
# no full, source-labelled copy of the inputs is staged for pooled output.
.ofst_input_block <- function(path, from, to, keys) {
  dt <- tryCatch(.normalize_dt(.ofst_read_part(path, from = from, to = to)),
    error = function(e) {
      if (inherits(e, "ofst_scratch_error")) stop(e)
      .ofst_abort("invalid input '", path, "', rows ", from, "-", to, ": ", conditionMessage(e))
    })
  if (anyNA(dt$seqnames) || anyNA(dt$start))
    .ofst_abort("input '", path, "' contains missing chromosome/start keys. Correct these before positional filtering.")
  if (any(!is.finite(dt$score) & !is.na(dt$score)) || any(dt$score < 0, na.rm = TRUE))
    .ofst_abort("input '", path, "' has negative or infinite scores. Filtering requires non-negative finite scores; missing scores contribute zero.")
  data.table::set(dt, j = "score", value = as.double(dt$score))
  dt[is.na(score), score := 0]
  dt[, c(keys, "score"), with = FALSE]
}

.ofst_batch_budget <- function(paths, rows, requested) {
  # Rows remain the hard ceiling. The Linux estimate only REDUCES it, reserving
  # space for input, bind/group/sort workspaces, output and serialization.
  if (requested < 1e6 || !file.exists("/proc/meminfo")) return(requested)
  available <- tryCatch({
    line <- readLines("/proc/meminfo", warn = FALSE)
    as.double(sub("^MemAvailable:[[:space:]]*([0-9]+).*", "\\1", line[grepl("^MemAvailable:", line)])) * 1024
  }, error = function(e) NA_real_)
  if (length(available) != 1L || !is.finite(available) || available <= 0) return(requested)
  # A container can have a much smaller limit than /proc/meminfo reports.
  for (pair in list(c("/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current"),
                    c("/sys/fs/cgroup/memory/memory.limit_in_bytes", "/sys/fs/cgroup/memory/memory.usage_in_bytes"))) {
    if (all(file.exists(pair))) {
      value <- suppressWarnings(vapply(pair, function(p) as.double(readLines(p, n = 1L, warn = FALSE)), 0))
      if (all(is.finite(value))) available <- min(available, max(0, value[1L] - value[2L]))
    }
  }
  bytes <- 1
  for (i in which(rows > 0)) {
    sample <- .ofst_read_part(paths[i], from = 1, to = min(rows[i], 4096))
    bytes <- max(bytes, as.double(object.size(sample)) / max(1, nrow(sample)))
  }
  budget <- max(1, min(requested, floor(available / (8 * bytes))))
  if (budget < requested) message("OFST batch merge: reducing row budget from ", .ofst_number(requested),
    " to ", .ofst_number(budget), " using sampled uncompressed size and available RAM (8-fold workspace allowance).")
  budget
}

.ofst_run_leaves <- function(run) {
  if (is.null(run)) return(list())
  if (!is.null(run$path)) return(list(run))
  c(.ofst_run_leaves(run$children[[1L]]), .ofst_run_leaves(run$children[[2L]]))
}

.ofst_write_run <- function(dt, keys, budget, scratch, depth = 0L, prefix = 0) {
  if (!nrow(dt)) return(NULL)
  if (nrow(dt) <= budget) {
    .downcast_score_if_safe(dt)
    data.table::setorderv(dt, keys)
    return(list(path = .ofst_spill(dt, scratch, "batch-"), rows = as.double(nrow(dt)),
                depth = depth, prefix = prefix))
  }
  if (depth >= 32L)
    .ofst_abort("a chromosome/start group needs more than filter_chunk_rows=", .ofst_number(budget),
                " distinct rows. Increase filter_chunk_rows if RAM permits.")
  bucket <- .ofst_partition_id(dt, depth)
  occupied <- unique(bucket)
  if (length(occupied) == 1L) {
    # Do not duplicate a large, unsplittable position at every hash bit.
    bucket <- NULL
    children <- list(NULL, NULL)
    children[occupied + 1L] <- list(.ofst_write_run(dt, keys, budget, scratch,
      depth + 1L, prefix + occupied * 2^depth))
    return(list(children = children, depth = depth, prefix = prefix))
  }
  children <- lapply(0:1, function(i) {
    index <- which(bucket == i)
    .ofst_write_run(dt[index], keys, budget, scratch, depth + 1L, prefix + i * 2^depth)
  })
  list(children = children, depth = depth, prefix = prefix)
}

.ofst_split_run <- function(run, keys, budget, scratch) {
  if (is.null(run) || is.null(run$path)) return(run)
  dt <- .ofst_read_part(run$path)
  bucket <- .ofst_partition_id(dt, run$depth)
  children <- lapply(0:1, function(i) {
    index <- which(bucket == i)
    .ofst_write_run(dt[index], keys, budget, scratch, run$depth + 1L, run$prefix + i * 2^run$depth)
  })
  unlink(run$path)
  list(children = children, depth = run$depth, prefix = run$prefix)
}

.ofst_combine_runs <- function(a, b, keys, budget, scratch) {
  if (is.null(a)) return(b)
  if (is.null(b)) return(a)
  if (a$depth != b$depth || a$prefix != b$prefix)
    .ofst_abort("internal batch partition mismatch; no result returned.")
  if (!is.null(a$path) && !is.null(b$path)) {
    dt <- data.table::rbindlist(list(.ofst_read_part(a$path), .ofst_read_part(b$path)), use.names = TRUE)
    dt <- .ofst_reduce_scores(dt, keys)
    result <- .ofst_write_run(dt, keys, budget, scratch, a$depth, a$prefix)
    # Never remove consumed intermediates until their replacement is complete.
    unlink(c(a$path, b$path))
    return(result)
  }
  a <- .ofst_split_run(a, keys, budget, scratch)
  b <- .ofst_split_run(b, keys, budget, scratch)
  children <- lapply(1:2, function(i) .ofst_combine_runs(a$children[[i]], b$children[[i]], keys, budget, scratch))
  list(children = children, depth = a$depth, prefix = a$prefix)
}

.ofst_batch_runs <- function(paths, row_counts, keys, source_col, budget, scratch, verbose = TRUE) {
  levels <- pending <- list()
  pending_rows <- 0
  template <- NULL
  batch <- 0L
  flush <- function() {
    if (!length(pending)) return(invisible(NULL))
    combined <- data.table::rbindlist(pending, use.names = TRUE)
    pending <<- list()
    pending_rows <<- 0
    combined <- .ofst_reduce_scores(combined, c(keys, source_col))
    batch <<- batch + 1L
    if (verbose) message("OFST batch merge: batch ", batch, " collapsed to ", .ofst_number(nrow(combined)), " rows before writing.")
    run <- .ofst_write_run(combined, c(keys, source_col), budget, scratch)
    combined <- NULL
    level <- 1L
    while (level <= length(levels) && !is.null(levels[[level]])) {
      if (verbose) message("OFST batch merge: combining compacted runs at level ", level, ".")
      previous <- levels[[level]]
      levels[level] <<- list(NULL)
      run <- .ofst_combine_runs(previous, run, c(keys, source_col), budget, scratch)
      level <- level + 1L
    }
    levels[level] <<- list(run)
    if (budget >= 1e5) gc()
    invisible(NULL)
  }
  for (i in seq_along(paths)) {
    if (verbose) message("OFST rescue: reading input ", i, "/", length(paths), " (",
                          .ofst_number(row_counts[i]), " rows): ", paths[i])
    first <- 1
    while (first <= row_counts[i]) {
      last <- min(row_counts[i], first + budget - pending_rows - 1)
      dt <- .ofst_input_block(paths[i], first, last, keys)
      if (is.null(template)) template <- dt[0]
      if (length(source_col)) data.table::set(dt, j = source_col, value = rep.int(i, nrow(dt)))
      pending[[length(pending) + 1L]] <- dt
      pending_rows <- pending_rows + last - first + 1
      dt <- NULL
      first <- last + 1
      if (pending_rows >= budget) flush()
    }
  }
  flush()
  run <- NULL
  for (level in seq_along(levels)) if (!is.null(levels[[level]])) {
    run <- .ofst_combine_runs(run, levels[[level]], c(keys, source_col), budget, scratch)
    levels[level] <- list(NULL)
  }
  list(tree = run, template = template, batches = batch)
}

.ofst_batch_input_stats <- function(paths, row_counts, keys, budget, scratch, leaves, selection) {
  stats <- vector("list", length(paths))
  message("OFST removal summary: rereading inputs for requested per-study statistics; this adds disk I/O.")
  for (i in seq_along(paths)) {
    message("OFST removal summary: input ", i, "/", length(paths), ".")
    staged <- .ofst_batch_runs(paths[i], row_counts[i], keys, NULL, budget, scratch, verbose = FALSE)
    parts <- .ofst_run_leaves(staged$tree)
    total <- c(rows_before = 0, score_before = 0, rows_after = 0, score_after = 0)
    for (part in parts) {
      dt <- .ofst_read_part(part$path)
      total[1:2] <- total[1:2] + c(nrow(dt), sum(dt$score))
      # Hash-tree leaves are either disjoint or one contains the other.
      for (leaf in leaves) {
        depth <- min(part$depth, leaf$depth)
        if (part$prefix %% 2^depth != leaf$prefix %% 2^depth) next
        index <- if (leaf$depth <= part$depth) seq_len(nrow(dt)) else
          which(.ofst_partition_id(dt, 0L, 2^leaf$depth) == leaf$prefix)
        if (!length(index)) next
        p <- .ofst_read_part(leaf$positions)
        drop <- p[.ofst_removal_mask(p, selection), .(seqnames, start)]
        retained <- dt[index][!drop, on = c("seqnames", "start")]
        total[3:4] <- total[3:4] + c(nrow(retained), sum(retained$score))
      }
      dt <- NULL
      unlink(part$path)
    }
    stats[[i]] <- data.table::data.table(input_index = i, rows_before = total[1L], score_before = total[2L],
                                        rows_after = total[3L], score_after = total[4L])
  }
  data.table::rbindlist(stats)
}

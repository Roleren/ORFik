# Incremental OFST merge store.
#
# A merged master ("all.ofst") is the sum of many OFST inputs (libraries or study pools).
# The store keeps that sum together with every input's contribution, so inputs can be
# added, replaced or removed without re-merging everything. It is the wide
# "keys x libraries" table with the zeros left out, split by chromosome:
#
#   <store>/manifest.rds                     the single commit point (see .ofst_store_new_manifest)
#   <store>/manifest_libs.csv                human-readable copy of the library table (not read back)
#   <store>/<part>/keys.fst                  seq, <key columns>; append-only, row number = row id
#   <store>/<part>/order_g<gen>.fst          row ids in key order (rewritten per update)
#   <store>/<part>/total_g<gen>.fst          summed score per row id: the "all" column
#   <store>/<part>/contrib_<seg>.fst         row, score; append-only, one row range per input
#
# Row ids are int32 per partition (one chromosome, or a pool of small contigs), so they stay
# far below 2^31 (human chromosome 1: ~1.1e8 keys). The contribution rows of a partition can
# exceed 2^31 in total, so they are split into segments of at most `max_contrib_rows` rows;
# the manifest stores (segment, from, to) as doubles.
#
# Crash safety: appended rows and new generation files become visible only when manifest.rds
# is replaced (write + rename). Rows beyond the committed counts are cut off before the next
# append, and files of other generations are deleted after a commit.
# Writers take <store>/.lock (a directory); readers must not run during an update.

.ofst_store_version <- 1L
.ofst_store_key_types <- c("start", "end", "width", "strand", "size")

.ofst_store_abort <- function(...) stop("OFST store: ", ..., call. = FALSE)
.ofst_store_msg <- function(verbose, ...) if (verbose) message("OFST store: ", ...)

.ofst_store_workers <- function(workers = NULL) {
  if (!is.null(workers)) return(max(1L, as.integer(workers)))
  w <- getOption("ORFik.ofst_store_workers", NULL)
  if (!is.null(w)) return(max(1L, as.integer(w)))
  .ofst_merge_workers()
}

.ofst_store_mclapply <- function(x, f, workers) {
  res <- if (workers > 1L && length(x) > 1L && .Platform$OS.type == "unix") {
    parallel::mclapply(x, function(i) tryCatch(f(i), error = function(e) e),
                       mc.cores = min(workers, length(x)), mc.preschedule = FALSE)
  } else lapply(x, function(i) tryCatch(f(i), error = function(e) e))
  bad <- vapply(res, function(r) inherits(r, "error") || is.null(r), logical(1))
  if (any(bad)) {
    r <- res[[which(bad)[1L]]]
    .ofst_store_abort(if (inherits(r, "error")) conditionMessage(r) else
      "a worker returned no result (out of memory?)")
  }
  res
}

# Named list of columns without copying (as.list() dispatches through BiocGenerics and may copy).
.ofst_cols <- function(x, cols) stats::setNames(lapply(cols, function(k) x[[k]]), cols)

.ofst_store_path <- function(store) {
  if (!is.character(store) || length(store) != 1L || is.na(store) || !nzchar(store))
    .ofst_store_abort("store must be one directory path.")
  normalizePath(store, mustWork = FALSE)
}

.ofst_store_manifest_path <- function(store) file.path(store, "manifest.rds")

#' Read the manifest of an OFST store
#'
#' @param store character, the store directory.
#' @return a list: \code{key_columns}, \code{seqlevels}, \code{chr_part} (chromosome to
#' partition), \code{parts} (per partition: committed key rows and generation),
#' \code{segments} (contribution segments and committed rows), \code{libs} (one row per
#' input and partition: segment, row range, non-zero rows, score sum, source file, state)
#' and \code{generation}.
#' @family ofst_store
#' @export
ofst_store_manifest <- function(store) {
  store <- .ofst_store_path(store)
  p <- .ofst_store_manifest_path(store)
  if (!file.exists(p)) .ofst_store_abort("no store at '", store, "' (manifest.rds missing).")
  m <- readRDS(p)
  if (!identical(m$version, .ofst_store_version))
    .ofst_store_abort("unsupported store version in '", store, "'.")
  m
}

.ofst_store_commit <- function(store, m) {
  tmp <- paste0(.ofst_store_manifest_path(store), ".tmp")
  saveRDS(m, tmp, compress = FALSE)
  if (!file.rename(tmp, .ofst_store_manifest_path(store)))
    .ofst_store_abort("could not commit the manifest of '", store, "'.")
  csv <- file.path(store, "manifest_libs.csv")
  try(data.table::fwrite(m$libs, paste0(csv, ".tmp")), silent = TRUE)
  try(file.rename(paste0(csv, ".tmp"), csv), silent = TRUE)
  invisible(m)
}

.ofst_store_lock <- function(store) {
  lock <- file.path(store, ".lock")
  if (!dir.create(lock, showWarnings = FALSE))
    .ofst_store_abort("store '", store, "' is locked by another writer (", lock,
                      "). Remove the directory only if no update is running.")
  writeLines(c(Sys.info()[["nodename"]], as.character(Sys.getpid()), format(Sys.time(), tz = "UTC")),
             file.path(lock, "owner"))
  function() unlink(lock, recursive = TRUE)
}

.ofst_store_part_dir <- function(store, part) file.path(store, part)
.ofst_store_gen_file <- function(store, part, what, gen)
  file.path(store, part, sprintf("%s_g%06d.fst", what, gen))
.ofst_store_seg_file <- function(store, part, seg)
  file.path(store, part, sprintf("contrib_%04d.fst", seg))

# Remove generation files that the manifest does not reference (after a commit or a crash).
.ofst_store_gc <- function(store, m) {
  for (i in seq_len(nrow(m$parts))) {
    p <- m$parts$part[i]; g <- m$parts$gen[i]
    files <- list.files(.ofst_store_part_dir(store, p), pattern = "^(order|total)_g[0-9]+\\.fst$", full.names = TRUE)
    keep <- c(.ofst_store_gen_file(store, p, "order", g), .ofst_store_gen_file(store, p, "total", g))
    unlink(setdiff(files, keep))
  }
  live_segs <- if (nrow(m$segments)) .ofst_store_seg_file(store, m$segments$part, m$segments$seg) else character()
  segs <- list.files(store, pattern = "^contrib_[0-9]+\\.fst$", recursive = TRUE, full.names = TRUE)
  unlink(setdiff(normalizePath(segs, mustWork = FALSE), normalizePath(live_segs, mustWork = FALSE)))
  invisible(NULL)
}

.ofst_store_nrow <- function(path) if (file.exists(path)) as.numeric(fst::metadata_fst(path)$nrOfRows) else 0

# Cut a file back to its committed rows (rows written by an update that never committed).
.ofst_store_truncate <- function(path, n) {
  have <- .ofst_store_nrow(path)
  if (have <= n) return(invisible(FALSE))
  d <- if (n > 0) fst::read_fst(path, from = 1, to = n, as.data.table = TRUE) else
    fst::read_fst(path, from = 1, to = 1, as.data.table = TRUE)[0]
  tmp <- paste0(path, ".trunc")
  fst::write_fst(d, tmp)
  file.rename(tmp, path)
  invisible(TRUE)
}

.ofst_store_append <- function(dt, path, committed) {
  if (!nrow(dt)) return(invisible(NULL))
  .ofst_store_truncate(path, committed)
  if (committed == 0) {
    fst::write_fst(dt, path)
  } else if ("append_rows_fst" %in% getNamespaceExports("fst")) {
    fst::append_rows_fst(dt, path)
  } else {   # CRAN fst: rewrite (correct, slower)
    old <- fst::read_fst(path, as.data.table = TRUE)
    tmp <- paste0(path, ".tmp")
    fst::write_fst(data.table::rbindlist(list(old, dt)), tmp)
    file.rename(tmp, path)
  }
  invisible(NULL)
}

.ofst_store_write_atomic <- function(dt, path) {
  tmp <- paste0(path, ".tmp")
  fst::write_fst(dt, tmp)
  if (!file.rename(tmp, path)) .ofst_store_abort("could not write '", path, "'.")
}

# Row ranges of each chromosome in an input (several runs per chromosome are allowed).
.ofst_store_runs <- function(path) {
  s <- fst::read_fst(path, "seqnames")[[1L]]
  if (is.character(s)) s <- factor(s)
  if (!length(s)) return(data.table::data.table(chr = character(), from = numeric(), to = numeric()))
  if (anyNA(s)) .ofst_store_abort("missing seqnames in '", path, "'.")
  r <- ofst_code_runs(as.integer(s))
  data.table::data.table(chr = levels(s)[r$value], from = c(1, head(r$end, -1) + 1), to = r$end)
}

# Key columns of an input: everything except seqnames, score and an all-1 width.
.ofst_store_input_keys <- function(path) {
  cols <- fst::metadata_fst(path)$columnNames
  if (!all(c("seqnames", "score") %in% cols)) .ofst_store_abort("'", path, "' has no seqnames/score column.")
  setdiff(cols, c("seqnames", "score", "width"))
}

# Rows of one input for a set of chromosomes, as integer key columns (seq code, keys...) +
# double score, grouped (sorted, unique).
.ofst_store_read_input <- function(path, runs, chr_codes, key_columns) {
  runs <- runs[order(from)]
  has_width <- "width" %in% fst::metadata_fst(path)$columnNames && !"width" %in% key_columns
  cols <- c(key_columns, "score", if (has_width) "width")
  x <- .ofst_read_ranges(path, runs$from, runs$to, cols)
  if (has_width) {
    if (anyNA(x$width) || any(x$width != 1L))
      .ofst_store_abort("'", path, "' has a width column with values other than 1.")
    x[, width := NULL]
  }
  data.table::set(x, j = "seq", value = rep.int(chr_codes[runs$chr], runs$to - runs$from + 1))
  if ("strand" %in% key_columns) {
    s <- x$strand
    if (is.character(s)) s <- factor(s, levels = .ofst_strand_levels)
    if (any(!levels(s) %in% .ofst_strand_levels)) .ofst_store_abort("strand must be '+', '-' or '*'.")
    data.table::set(x, j = "strand", value = .ofst_strand_to_sortkey[match(levels(s), .ofst_strand_levels)][as.integer(s)])
  }
  for (k in setdiff(key_columns, "strand")) {
    v <- x[[k]]
    if (!is.integer(v)) {
      if (!is.numeric(v) || any(v != round(v), na.rm = TRUE) || any(abs(v) > .int32_max, na.rm = TRUE))
        .ofst_store_abort("key column '", k, "' of '", path, "' is not integer valued.")
      data.table::set(x, j = k, value = as.integer(v))
    }
  }
  sc <- as.numeric(x$score)
  sc[is.na(sc)] <- 0
  data.table::set(x, j = "score", value = sc)
  kc <- c("seq", key_columns)
  # pre-sorted inputs (per-library P-site files are): one linear check instead of a group-by
  if (ofst_keys_sorted_unique(.ofst_cols(x, kc))) {
    data.table::setcolorder(x, c(kc, "score"))
    g <- x
  } else {
    g <- x[, .(score = sum(score)), keyby = kc]
    data.table::setattr(g, "sorted", NULL)
  }
  g[score != 0]
}

.ofst_store_new_manifest <- function(key_columns, max_contrib_rows) {
  list(version = .ofst_store_version, key_columns = key_columns, seqlevels = character(),
       chr_part = data.table::data.table(chr = character(), part = character()),
       parts = data.table::data.table(part = character(), n_keys = numeric(), gen = integer()),
       segments = data.table::data.table(part = character(), seg = integer(), n_rows = numeric()),
       libs = data.table::data.table(lib = character(), part = character(), seg = integer(),
                                     from = numeric(), to = numeric(), nnz = numeric(),
                                     score_sum = numeric(), source = character(), source_size = numeric(),
                                     source_mtime = character(), state = character(), generation = integer()),
       generation = 0L, max_contrib_rows = max_contrib_rows)
}

#' Build an incremental OFST merge store
#'
#' Creates a store that holds the merged ("all") track of many OFST inputs together with
#' each input's contribution, so that inputs can later be added
#' (\code{ofst_store_add}), replaced (\code{ofst_store_replace}) or removed
#' (\code{ofst_store_remove}) without re-merging everything. The merged track is written
#' with \code{ofst_store_export_all}; the whole matrix can be loaded into RAM with
#' \code{ofst_store_read_matrix}.
#'
#' Layout: one partition per chromosome (chromosomes with fewer than
#' \code{partition_min_rows} input rows share one pooled partition), each with an
#' append-only key dictionary, a key order, the total ("all") column and the sparse
#' per-input contributions (row id, score). \code{manifest.rds} is the single commit point.
#' @param store character, directory of the new store (must not exist or be empty).
#' @param paths character, OFST files (P-site libraries or study pools). Every file needs
#' \code{seqnames} and \code{score}; all other columns except an all-1 \code{width} are keys
#' and must be integer valued (strand as \code{+ - *}). Inputs need not be sorted.
#' @param libs character, unique names of the inputs, default the file names without
#' extension.
#' @param batch_size integer, inputs read per batch (bounds memory), default 100.
#' @param partition_min_rows numeric, default 1e6: chromosomes with fewer input rows (in the
#' first batch) go to the pooled partition.
#' @param max_contrib_rows numeric, default 1e9: rows per contribution segment file.
#' @param workers integer, parallel partitions (forked); default
#' \code{getOption("ORFik.ofst_store_workers")} or half the data.table threads.
#' @param verbose logical, default TRUE.
#' @return invisible manifest (see \code{ofst_store_manifest}).
#' @family ofst_store
#' @export
#' @examples
#' lib <- function(seed) {
#'   set.seed(seed)
#'   data.frame(seqnames = factor(sample(c("1", "2"), 50, TRUE)), start = sample(100L, 50, TRUE),
#'              strand = factor(sample(c("+", "-"), 50, TRUE), levels = c("+", "-", "*")),
#'              size = 28L, score = sample(5L, 50, TRUE))
#' }
#' paths <- file.path(tempdir(), paste0("lib", 1:3, ".ofst"))
#' for (i in 1:3) fst::write_fst(lib(i), paths[i])
#' store <- tempfile("store_example")
#' ofst_store_build(store, paths[1:2])
#' ofst_store_add(store, paths[3])
#' merged <- ofst_store_export_all(store)   # = ofst_merge(paths, keep_all_scores = FALSE)
#' ofst_store_remove(store, "lib1")
#' ofst_store_read_matrix(store)
ofst_store_build <- function(store, paths, libs = sub("\\.ofst$", "", basename(paths)),
                             batch_size = 100L, partition_min_rows = 1e6, max_contrib_rows = 1e9,
                             workers = NULL, verbose = TRUE) {
  store <- .ofst_store_path(store)
  if (!is.character(paths) || !length(paths) || anyNA(paths)) .ofst_store_abort("paths must be file paths.")
  if (length(libs) != length(paths) || anyDuplicated(libs) || anyNA(libs) || any(!nzchar(libs)))
    .ofst_store_abort("libs must be unique non-empty names, one per path.")
  if (dir.exists(store) && length(list.files(store, all.files = TRUE, no.. = TRUE)))
    .ofst_store_abort("'", store, "' exists and is not empty.")
  if (!is.numeric(max_contrib_rows) || max_contrib_rows < 1 || max_contrib_rows > .int32_max)
    .ofst_store_abort("max_contrib_rows must be between 1 and 2^31 - 1.")
  keys <- lapply(paths, .ofst_store_input_keys)
  key_columns <- keys[[1L]]
  if (!all(vapply(keys, setequal, logical(1), key_columns)))
    .ofst_store_abort("all inputs need the same key columns (", paste(key_columns, collapse = ", "), ").")
  if (any(!key_columns %in% .ofst_store_key_types))
    .ofst_store_abort("unsupported key columns: ", paste(setdiff(key_columns, .ofst_store_key_types), collapse = ", "),
                      " (the store handles integer keys only, no cigar).")
  dir.create(store, recursive = TRUE, showWarnings = FALSE)
  m <- .ofst_store_new_manifest(key_columns, max_contrib_rows)
  m$partition_min_rows <- partition_min_rows
  .ofst_store_commit(store, m)
  batches <- split(seq_along(paths), ceiling(seq_along(paths) / batch_size))
  for (b in seq_along(batches)) {
    i <- batches[[b]]
    .ofst_store_msg(verbose, "build batch ", b, "/", length(batches), " (", length(i), " inputs)")
    m <- ofst_store_add(store, paths[i], libs[i], workers = workers, verbose = verbose)
  }
  invisible(ofst_store_manifest(store))
}

#' Add inputs to an OFST store
#'
#' @inheritParams ofst_store_build
#' @param store character, store directory (from \code{ofst_store_build}).
#' @param paths character, OFST files to add.
#' @param libs character, unique names of the new inputs (not live in the store).
#' @return invisible manifest after the update.
#' @family ofst_store
#' @export
ofst_store_add <- function(store, paths, libs = sub("\\.ofst$", "", basename(paths)),
                           workers = NULL, verbose = TRUE) {
  .ofst_store_update(store, add_paths = paths, add_libs = libs, remove_libs = character(),
                     workers = workers, verbose = verbose)
}

#' Remove inputs from an OFST store
#'
#' Subtracts the stored contributions; the original files are not needed (they may be
#' overwritten already, e.g. after a reshift). Removed rows stay in the contribution files
#' until \code{ofst_store_compact}.
#' @inheritParams ofst_store_add
#' @param libs character, names of live inputs to remove.
#' @return invisible manifest after the update.
#' @family ofst_store
#' @export
ofst_store_remove <- function(store, libs, workers = NULL, verbose = TRUE) {
  .ofst_store_update(store, add_paths = character(), add_libs = character(), remove_libs = libs,
                     workers = workers, verbose = verbose)
}

#' Replace inputs of an OFST store
#'
#' Removes the stored contribution of each named input and adds the new file under the same
#' name, in one transaction (one manifest commit).
#' @inheritParams ofst_store_add
#' @param libs character, names of live inputs, one per path.
#' @return invisible manifest after the update.
#' @family ofst_store
#' @export
ofst_store_replace <- function(store, paths, libs = sub("\\.ofst$", "", basename(paths)),
                               workers = NULL, verbose = TRUE) {
  .ofst_store_update(store, add_paths = paths, add_libs = libs, remove_libs = libs,
                     workers = workers, verbose = verbose)
}

.ofst_store_update <- function(store, add_paths, add_libs, remove_libs, workers = NULL, verbose = TRUE) {
  store <- .ofst_store_path(store)
  t0 <- Sys.time()
  m <- ofst_store_manifest(store)
  unlock <- .ofst_store_lock(store)
  on.exit(unlock(), add = TRUE)
  workers <- .ofst_store_workers(workers)
  live <- m$libs[state == "live"]
  if (length(add_paths) != length(add_libs) || anyDuplicated(add_libs) || anyNA(add_libs) || any(!nzchar(add_libs)))
    .ofst_store_abort("one unique, non-empty name per added path is required.")
  missing <- add_paths[!file.exists(add_paths)]
  if (length(missing)) .ofst_store_abort("inputs do not exist: ", paste(missing, collapse = ", "))
  not_live <- setdiff(remove_libs, live$lib)
  if (length(not_live)) .ofst_store_abort("not live in the store: ", paste(head(not_live, 10), collapse = ", "))
  clash <- intersect(setdiff(add_libs, remove_libs), live$lib)
  if (length(clash)) .ofst_store_abort("already in the store (use ofst_store_replace): ", paste(head(clash, 10), collapse = ", "))
  for (p in add_paths) if (!setequal(.ofst_store_input_keys(p), m$key_columns))
    .ofst_store_abort("'", p, "' does not have the store's key columns (", paste(m$key_columns, collapse = ", "), ").")

  # chromosome runs of the new inputs (one pass over each seqnames column)
  runs <- .ofst_store_mclapply(seq_along(add_paths), function(i) {
    r <- .ofst_store_runs(add_paths[i]); data.table::set(r, j = "input", value = rep.int(i, nrow(r))); r
  }, workers)
  runs <- data.table::rbindlist(runs)
  if (!nrow(runs)) runs <- data.table::data.table(chr = character(), from = numeric(), to = numeric(), input = integer())

  # new chromosomes: seq codes and partitions
  new_chr <- setdiff(unique(runs$chr), m$chr_part$chr)
  if (length(new_chr)) {
    rows <- runs[chr %in% new_chr, .(rows = sum(to - from + 1)), by = chr]
    min_rows <- if (is.null(m$partition_min_rows)) 1e6 else m$partition_min_rows
    rows[, part := ifelse(rows >= min_rows, "", "pooled")]
    n_old <- nrow(m$parts)
    own <- which(rows$part == "")
    if (length(own)) rows$part[own] <- sprintf("part_%04d", n_old + seq_along(own))
    m$chr_part <- data.table::rbindlist(list(m$chr_part, rows[, .(chr, part)]))
    m$seqlevels <- c(m$seqlevels, setdiff(new_chr, m$seqlevels))
    new_parts <- setdiff(unique(rows$part), m$parts$part)
    if (length(new_parts))
      m$parts <- data.table::rbindlist(list(m$parts, data.table::data.table(part = new_parts, n_keys = 0, gen = 0L)))
  }
  chr_codes <- stats::setNames(seq_along(m$seqlevels), m$seqlevels)
  runs[, part := m$chr_part$part[match(chr, m$chr_part$chr)]]
  rm_rows <- live[lib %in% remove_libs]
  todo <- sort(unique(c(runs$part, rm_rows$part)))
  new_gen <- m$generation + 1L
  .ofst_store_msg(verbose, "update generation ", new_gen, ": +", length(add_paths), " / -", length(remove_libs),
                  " inputs over ", length(todo), " partitions, ", min(workers, length(todo)), " worker(s)")
  m_ro <- m
  res <- .ofst_store_mclapply(todo, function(p) {
    .ofst_store_update_part(store, m_ro, p, new_gen, add_paths, add_libs, runs[part == p], chr_codes,
                            rm_rows[part == p])
  }, workers)
  names(res) <- todo

  # commit
  for (p in todo) {
    r <- res[[p]]
    data.table::set(m$parts, i = which(m$parts$part == p), j = c("n_keys", "gen"), value = list(r$n_keys, new_gen))
    m$segments <- data.table::rbindlist(list(m$segments[part != p], r$segments))
  }
  m$libs[lib %in% remove_libs & state == "live", state := "removed"]
  added <- data.table::rbindlist(lapply(res, `[[`, "libs"))
  if (nrow(added)) {
    info <- file.info(add_paths)
    i <- match(added$lib, add_libs)
    added[, `:=`(source = normalizePath(add_paths[i]), source_size = info$size[i],
                 source_mtime = format(info$mtime[i], "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
                 state = "live", generation = new_gen)]
    m$libs <- data.table::rbindlist(list(m$libs, added), use.names = TRUE)
  }
  # inputs without rows still get one manifest row, so they count as live
  empty <- setdiff(add_libs, added$lib)
  if (length(empty))
    m$libs <- data.table::rbindlist(list(m$libs, data.table::data.table(
      lib = empty, part = NA_character_, seg = NA_integer_, from = NA_real_, to = NA_real_, nnz = 0,
      score_sum = 0, source = normalizePath(add_paths[match(empty, add_libs)]),
      source_size = file.info(add_paths[match(empty, add_libs)])$size,
      source_mtime = format(file.info(add_paths[match(empty, add_libs)])$mtime, "%Y-%m-%dT%H:%M:%OS3Z", tz = "UTC"),
      state = "live", generation = new_gen)), use.names = TRUE)
  m$generation <- new_gen
  .ofst_store_commit(store, m)
  .ofst_store_gc(store, m)
  .ofst_store_msg(verbose, "generation ", new_gen, " committed in ",
                  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
  invisible(m)
}

# Load a partition's committed state: keys (list of int columns), order, total.
.ofst_store_load_part <- function(store, m, p, need_keys = TRUE) {
  pr <- m$parts[part == p]
  n <- if (nrow(pr)) pr$n_keys else 0
  kc <- c("seq", m$key_columns)
  if (n == 0) {
    keys <- if (need_keys) stats::setNames(rep(list(integer()), length(kc)), kc)
    return(list(n = 0, keys = keys, order = integer(), total = numeric()))
  }
  keys <- if (need_keys) .ofst_cols(fst::read_fst(file.path(store, p, "keys.fst"), from = 1, to = n), kc)
  ord <- fst::read_fst(.ofst_store_gen_file(store, p, "order", pr$gen))[[1L]]
  total <- as.numeric(fst::read_fst(.ofst_store_gen_file(store, p, "total", pr$gen))[[1L]])
  if (length(ord) != n || length(total) != n)
    .ofst_store_abort("partition '", p, "' is inconsistent with the manifest (", length(ord), "/",
                      length(total), " vs ", n, " keys).")
  list(n = n, keys = keys, order = ord, total = total)
}

# Contributions (row, score) of manifest rows h (one partition), concatenated in h order.
.ofst_store_read_contrib <- function(store, h) {
  h <- h[!is.na(seg) & nnz > 0]
  if (!nrow(h)) return(data.table::data.table(row = integer(), score = numeric()))
  pieces <- lapply(split(seq_len(nrow(h)), h$seg), function(i) {
    hh <- h[i]
    o <- order(hh$from)
    d <- .ofst_read_ranges(.ofst_store_seg_file(store, hh$part[1L], hh$seg[1L]), hh$from[o], hh$to[o], c("row", "score"))
    list(i = i[o], d = d)
  })
  idx <- unlist(lapply(pieces, `[[`, "i"))
  d <- data.table::rbindlist(lapply(pieces, `[[`, "d"))
  if (!identical(idx, seq_len(nrow(h))))   # back to h order (stable)
    d <- d[order(rep.int(idx, h$nnz[idx]), method = "radix")]
  data.table::set(d, j = "score", value = as.numeric(d$score))
  d
}

.ofst_store_update_part <- function(store, m, p, gen, add_paths, add_libs, runs, chr_codes, rm_rows) {
  dir.create(file.path(store, p), showWarnings = FALSE)
  st <- .ofst_store_load_part(store, m, p, need_keys = nrow(runs) > 0)
  total <- st$total
  # removals: subtract stored contributions
  if (nrow(rm_rows)) {
    d <- .ofst_store_read_contrib(store, rm_rows)
    total <- ofst_store_scatter_add(total, d$row, d$score, -1, TRUE)
    d <- NULL
  }
  segs <- m$segments[part == p]
  out_libs <- NULL
  n_keys <- st$n
  ord <- st$order
  if (nrow(runs)) {
    kc <- c("seq", m$key_columns)
    inputs <- sort(unique(runs$input))
    lib_tabs <- lapply(inputs, function(i)
      .ofst_store_read_input(add_paths[i], runs[input == i], chr_codes, m$key_columns))
    # match against the stored keys; unseen keys get new row ids (appended in key order)
    rows <- lapply(lib_tabs, function(g) ofst_store_match(st$keys, st$order, .ofst_cols(g, kc)))
    unseen <- data.table::rbindlist(Map(function(g, r) g[r == 0L, kc, with = FALSE], lib_tabs, rows))
    unseen <- unique(unseen)
    data.table::setorderv(unseen, kc)
    n_new <- nrow(unseen)
    if (st$n + n_new > .int32_max) .ofst_store_abort("partition '", p, "' would exceed 2^31 - 1 keys.")
    keys <- st$keys
    if (n_new) {
      new_ord <- seq_len(n_new)
      for (j in seq_along(lib_tabs)) {
        miss <- rows[[j]] == 0L
        if (any(miss)) rows[[j]][miss] <- as.integer(st$n) + ofst_store_match(.ofst_cols(unseen, kc), new_ord,
                                                                  .ofst_cols(lib_tabs[[j]][miss, kc, with = FALSE], kc))
      }
      .ofst_store_append(unseen, file.path(store, p, "keys.fst"), st$n)
      keys <- Map(c, keys, .ofst_cols(unseen, kc))
      ord <- ofst_store_merge_order(keys, ord, as.integer(st$n + new_ord))
      total <- c(total, numeric(n_new))
    }
    n_keys <- st$n + n_new
    # contributions, sorted by row id; segments capped at max_contrib_rows
    cur <- if (nrow(segs)) max(segs$seg) else 0L
    seg_rows <- if (nrow(segs)) segs$n_rows[segs$seg == cur] else Inf
    committed <- seg_rows
    batch <- list(); out_libs <- vector("list", length(lib_tabs))
    flush <- function() {
      if (length(batch)) .ofst_store_append(data.table::rbindlist(batch), .ofst_store_seg_file(store, p, cur), committed)
      batch <<- list()
    }
    for (j in seq_along(lib_tabs)) {
      r <- rows[[j]]; s <- lib_tabs[[j]]$score
      if (any(r == 0L)) .ofst_store_abort("internal: unmatched key in partition '", p, "'.")
      o <- order(r)
      cd <- data.table::data.table(row = as.integer(r[o]), score = as.numeric(s[o]))
      nn <- nrow(cd)
      if (nn > m$max_contrib_rows) .ofst_store_abort("input '", add_libs[inputs[j]], "' has more rows in partition '", p,
                                                     "' than max_contrib_rows.")
      if (seg_rows + nn > m$max_contrib_rows) {   # new segment
        flush()
        cur <- cur + 1L; seg_rows <- 0; committed <- 0
        segs <- data.table::rbindlist(list(segs, data.table::data.table(part = p, seg = cur, n_rows = 0)))
      }
      out_libs[[j]] <- data.table::data.table(lib = add_libs[inputs[j]], part = p, seg = cur, from = seg_rows + 1,
                                              to = seg_rows + nn, nnz = nn, score_sum = sum(s))
      batch[[length(batch) + 1L]] <- cd
      seg_rows <- seg_rows + nn
      data.table::set(segs, i = which(segs$seg == cur), j = "n_rows", value = seg_rows)
      total <- ofst_store_scatter_add(total, cd$row, cd$score, 1, FALSE)
    }
    flush()
    out_libs <- data.table::rbindlist(out_libs)
  }
  fst::write_fst(data.frame(row = as.integer(ord)), .ofst_store_gen_file(store, p, "order", gen))
  fst::write_fst(data.frame(score = total), .ofst_store_gen_file(store, p, "total", gen))
  list(n_keys = n_keys, segments = segs, libs = out_libs)
}

# Scores of one partition for all inputs (stored total) or a subset of live inputs.
.ofst_store_part_scores <- function(store, m, p, st, libs = NULL) {
  if (is.null(libs)) return(st$total)
  h <- m$libs[state == "live" & part == p & lib %in% libs]
  acc <- numeric(st$n)
  if (nrow(h)) {
    d <- .ofst_store_read_contrib(store, h)
    acc <- ofst_store_scatter_add(acc, d$row, d$score, 1, FALSE)
  }
  acc
}

# Sorted rows (score != 0) of one partition as a data.table(seq, keys..., score).
.ofst_store_part_table <- function(keys, ord, score) {
  keep <- ord[score[ord] != 0]
  d <- data.table::as.data.table(lapply(keys, `[`, keep))
  data.table::set(d, j = "score", value = score[keep])
  d
}

# Pieces (one per partition) -> the OFST table, rows and levels exactly as ofst_merge().
.ofst_store_assemble <- function(pieces, seqlevels, key_columns) {
  d <- data.table::rbindlist(pieces)
  pieces <- NULL
  if (!nrow(d)) {
    out <- data.table::data.table(seqnames = factor(character()))
    for (k in key_columns) data.table::set(out, j = k, value = integer())
    if ("strand" %in% key_columns) data.table::set(out, j = "strand", value = factor(character(), levels = .ofst_strand_levels))
    data.table::set(out, j = "score", value = integer())
    return(out[])
  }
  # chromosomes in character label order (as the sorted merge); each one is a contiguous block
  # within its partition, so ordering the blocks is enough
  runs <- ofst_code_runs(d$seq)
  chrs <- seqlevels[runs$value]
  chr_order <- order(chrs, method = "radix")
  if (is.unsorted(chr_order) || anyDuplicated(runs$value)) {
    if (anyDuplicated(runs$value)) .ofst_store_abort("internal: a chromosome is split in the export.")
    from <- c(1, head(runs$end, -1) + 1)
    idx <- unlist(lapply(chr_order, function(i) seq.int(from[i], runs$end[i])), use.names = FALSE)
    d <- d[idx]
    runs <- ofst_code_runs(d$seq)
    chrs <- seqlevels[runs$value]
  }
  data.table::set(d, j = "seq", value = structure(match(seqlevels[d$seq], chrs), levels = chrs, class = "factor"))
  data.table::setnames(d, "seq", "seqnames")
  if ("strand" %in% key_columns)
    data.table::set(d, j = "strand", value = structure(.ofst_sortkey_to_strand[d$strand],
                                                       levels = .ofst_strand_levels, class = "factor"))
  data.table::setcolorder(d, c("seqnames", key_columns, "score"))   # extra columns (wide) follow
  .downcast_score_if_safe(d)
  suppressMessages(.ofst_finalize(d))
}

#' Export the merged track of an OFST store
#'
#' Writes (or returns) the ordinary sorted OFST table: one row per key with a non-zero summed
#' score, identical to \code{ofst_merge(paths, keep_all_scores = FALSE)} of the live inputs.
#' @inheritParams ofst_store_add
#' @param store character, store directory, or an in-RAM matrix from
#' \code{ofst_store_read_matrix(as = "matrix")}.
#' @param out_file character or NULL (default). If given, the table is written there (via a
#' temporary file and rename, as \code{mergeLibs}) and invisibly returned.
#' @param libs character or NULL (default, all live inputs). A subset of live inputs gives a
#' subset track (e.g. inhibitor splits or study pools) without touching the rest.
#' @return a data.table (seqnames, key columns, score).
#' @family ofst_store
#' @export
ofst_store_export_all <- function(store, out_file = NULL, libs = NULL, workers = NULL, verbose = TRUE) {
  t0 <- Sys.time()
  if (inherits(store, "ofst_store_matrix")) {
    x <- store
    if (!is.null(libs) && length(setdiff(libs, x$libs))) .ofst_store_abort("unknown inputs: ", paste(head(setdiff(libs, x$libs)), collapse = ", "))
    pieces <- lapply(x$parts, function(pp) {
      sel <- if (is.null(libs)) rep(TRUE, length(pp$lib)) else x$libs[pp$lib] %in% libs
      score <- ofst_store_scatter_add(numeric(length(pp$keys[[1L]])), pp$row[sel], pp$score[sel], 1, FALSE)
      .ofst_store_part_table(pp$keys, ofst_store_key_order(pp$keys), score)
    })
    dt <- .ofst_store_assemble(pieces, x$seqlevels, x$key_columns)
  } else {
    store <- .ofst_store_path(store)
    m <- ofst_store_manifest(store)
    if (!is.null(libs)) {
      bad <- setdiff(libs, m$libs[state == "live"]$lib)
      if (length(bad)) .ofst_store_abort("not live in the store: ", paste(head(bad, 10), collapse = ", "))
    }
    parts <- m$parts[n_keys > 0]$part
    pieces <- .ofst_store_mclapply(parts, function(p) {
      st <- .ofst_store_load_part(store, m, p)
      .ofst_store_part_table(st$keys, st$order, .ofst_store_part_scores(store, m, p, st, libs))
    }, .ofst_store_workers(workers))
    dt <- .ofst_store_assemble(pieces, m$seqlevels, m$key_columns)
  }
  pieces <- NULL
  .ofst_store_msg(verbose, "exported ", .ofst_number(nrow(dt)), " rows in ",
                  round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), " s")
  if (!is.null(out_file)) {
    dir.create(dirname(out_file), recursive = TRUE, showWarnings = FALSE)
    suppressMessages(.ofst_save_merge(dt, out_file))
    return(invisible(dt))
  }
  dt
}

#' Verify an OFST store
#'
#' Re-sums the stored contributions of the live inputs per partition and compares them with
#' the stored totals and the manifest score sums. Run before publishing an export.
#' @inheritParams ofst_store_add
#' @return TRUE, or stops with the first mismatch.
#' @family ofst_store
#' @export
ofst_store_verify <- function(store, workers = NULL) {
  store <- .ofst_store_path(store)
  m <- ofst_store_manifest(store)
  ok <- .ofst_store_mclapply(m$parts$part, function(p) {
    st <- .ofst_store_load_part(store, m, p, need_keys = FALSE)
    acc <- .ofst_store_part_scores(store, m, p, st, libs = m$libs[state == "live"]$lib)
    h <- m$libs[state == "live" & part == p]
    if (!isTRUE(all.equal(acc, st$total, tolerance = 0)))
      .ofst_store_abort("partition '", p, "': totals differ from the summed contributions.")
    if (sum(st$total) != sum(h$score_sum))
      .ofst_store_abort("partition '", p, "': total score differs from the manifest score sums.")
    if (any(st$total < 0)) .ofst_store_abort("partition '", p, "': negative totals.")
    if (st$n && (anyDuplicated(st$order) || min(st$order) < 1 || max(st$order) > st$n))
      .ofst_store_abort("partition '", p, "': key order is not a permutation.")
    TRUE
  }, .ofst_store_workers(workers))
  all(unlist(ok))
}

#' Compact an OFST store
#'
#' Rewrites the contribution files of partitions where removed inputs hold more than
#' \code{min_dead_fraction} of the rows, keeping only live inputs (one new segment series per
#' partition). Keys and row ids are unchanged.
#' @inheritParams ofst_store_add
#' @param min_dead_fraction numeric, default 0.2.
#' @return invisible manifest after compaction.
#' @family ofst_store
#' @export
ofst_store_compact <- function(store, min_dead_fraction = 0.2, workers = NULL, verbose = TRUE) {
  store <- .ofst_store_path(store)
  m <- ofst_store_manifest(store)
  unlock <- .ofst_store_lock(store)
  on.exit(unlock(), add = TRUE)
  stats <- m$libs[!is.na(part), .(live = sum(nnz[state == "live"]), dead = sum(nnz[state != "live"])), by = part]
  todo <- stats[dead > 0 & dead / pmax(1, live + dead) > min_dead_fraction]$part
  if (!length(todo)) { .ofst_store_msg(verbose, "nothing to compact"); return(invisible(m)) }
  first_new <- max(c(0L, m$segments$seg)) + 1L   # new segment numbers never reuse a live file name
  res <- .ofst_store_mclapply(todo, function(p) {
    h <- m$libs[state == "live" & part == p & !is.na(seg)]
    cur <- first_new; rows <- 0; out <- vector("list", nrow(h)); segs <- list()
    for (i in seq_len(nrow(h))) {
      if (rows + h$nnz[i] > m$max_contrib_rows || i == 1L) {
        if (i > 1L) cur <- cur + 1L
        rows <- 0
      }
      d <- .ofst_store_read_contrib(store, h[i])
      .ofst_store_append(d, .ofst_store_seg_file(store, p, cur), rows)
      out[[i]] <- data.table::data.table(lib = h$lib[i], from = rows + 1, to = rows + h$nnz[i], seg = cur)
      rows <- rows + h$nnz[i]
      segs[[as.character(cur)]] <- data.table::data.table(part = p, seg = cur, n_rows = rows)
    }
    list(libs = data.table::rbindlist(out), segments = data.table::rbindlist(segs))
  }, .ofst_store_workers(workers))
  names(res) <- todo
  for (p in todo) {
    r <- res[[p]]
    m$segments <- data.table::rbindlist(list(m$segments[part != p], r$segments))
    i <- which(m$libs$state == "live" & m$libs$part %in% p & !is.na(m$libs$seg))
    j <- match(m$libs$lib[i], r$libs$lib)
    data.table::set(m$libs, i = i, j = c("seg", "from", "to"), value = list(r$libs$seg[j], r$libs$from[j], r$libs$to[j]))
    m$libs <- m$libs[!(part %in% p & state != "live")]
  }
  m$generation <- m$generation + 1L
  .ofst_store_commit(store, m)
  .ofst_store_gc(store, m)
  .ofst_store_msg(verbose, "compacted ", length(todo), " partition(s)")
  invisible(m)
}

# Key table in row-id order (no sorting), as OFST columns.
.ofst_store_keys_table <- function(keys, seqlevels, key_columns) {
  d <- data.table::as.data.table(keys)
  used <- sort(unique(d$seq))
  data.table::set(d, j = "seq", value = structure(match(d$seq, used), levels = seqlevels[used], class = "factor"))
  data.table::setnames(d, "seq", "seqnames")
  if ("strand" %in% key_columns)
    data.table::set(d, j = "strand", value = structure(.ofst_sortkey_to_strand[d$strand],
                                                       levels = .ofst_strand_levels, class = "factor"))
  d[]
}

# Key order of a list of integer key columns.
ofst_store_key_order <- function(keys) {
  if (!length(keys[[1L]])) return(integer())
  do.call(order, c(unname(keys), list(method = "radix")))
}

#' Read an OFST store into RAM as a sparse matrix
#'
#' The in-RAM mode: the whole store (or some inputs / chromosomes) as a keys x inputs sparse
#' matrix. With enough memory, updates and subset sums are plain vector operations on this
#' object; write it back with \code{ofst_store_write_matrix} and export it with
#' \code{ofst_store_export_all}.
#' @param store character, store directory.
#' @param libs character or NULL (default all live inputs).
#' @param seqnames character or NULL (default all chromosomes).
#' @param as character, "matrix" (default): an \code{ofst_store_matrix} list with
#' \code{key_columns}, \code{seqlevels}, \code{libs} and per partition \code{keys} (list of
#' integer columns: seq code, keys; strand as sort key 1 = "*", 2 = "+", 3 = "-") and the
#' triplets \code{row}, \code{lib} (index into \code{libs}), \code{score}.
#' "dgCMatrix": per partition a list(keys = OFST key table, m = Matrix::dgCMatrix with one
#' column per input). "wide": one OFST table with a score column per input plus the total
#' \code{score} (the \code{keep_all_scores = TRUE} layout of \code{ofst_merge}), rows sorted.
#' @param workers integer, parallel partitions.
#' @return see \code{as}.
#' @family ofst_store
#' @export
ofst_store_read_matrix <- function(store, libs = NULL, seqnames = NULL, as = c("matrix", "dgCMatrix", "wide"),
                                   workers = NULL) {
  as <- match.arg(as)
  store <- .ofst_store_path(store)
  m <- ofst_store_manifest(store)
  live <- m$libs[state == "live"]
  lib_names <- unique(live$lib)
  if (!is.null(libs)) {
    bad <- setdiff(libs, lib_names)
    if (length(bad)) .ofst_store_abort("not live in the store: ", paste(head(bad, 10), collapse = ", "))
    lib_names <- libs
  }
  parts <- m$parts[n_keys > 0]$part
  if (!is.null(seqnames)) parts <- intersect(parts, m$chr_part[chr %in% seqnames]$part)
  res <- .ofst_store_mclapply(parts, function(p) {
    st <- .ofst_store_load_part(store, m, p)
    h <- live[part == p & lib %in% lib_names & !is.na(seg)]
    d <- .ofst_store_read_contrib(store, h)
    keys <- st$keys
    li <- rep.int(match(h$lib, lib_names), h$nnz)
    row <- d$row; score <- d$score
    if (!is.null(seqnames) && length(unique(keys$seq)) > 1L) {   # pooled partition: keep the asked chromosomes
      want <- match(seqnames, m$seqlevels)
      keep_row <- keys$seq %in% want
      new_id <- cumsum(keep_row)
      sel <- keep_row[row]
      row <- new_id[row[sel]]; li <- li[sel]; score <- score[sel]
      keys <- lapply(keys, `[`, keep_row)
    }
    list(keys = keys, row = row, lib = li, score = score)
  }, .ofst_store_workers(workers))
  names(res) <- parts
  x <- structure(list(key_columns = m$key_columns, seqlevels = m$seqlevels, libs = lib_names, parts = res),
                 class = "ofst_store_matrix")
  if (as == "matrix") return(x)
  if (as == "dgCMatrix") {
    if (!requireNamespace("Matrix", quietly = TRUE)) .ofst_store_abort("as = 'dgCMatrix' needs the Matrix package.")
    return(lapply(x$parts, function(pp) {
      n <- length(pp$keys[[1L]])
      list(keys = .ofst_store_keys_table(pp$keys, x$seqlevels, x$key_columns),
           m = Matrix::sparseMatrix(i = pp$row, j = pp$lib, x = pp$score, dims = c(n, length(lib_names)),
                                    dimnames = list(NULL, lib_names)))
    }))
  }
  # wide: one OFST table with one column per input + total score
  if (any(lib_names %in% c("seq", "seqnames", m$key_columns, "score")))
    .ofst_store_abort("input names collide with OFST column names.")
  pieces <- lapply(x$parts, function(pp) {
    n <- length(pp$keys[[1L]])
    total <- ofst_store_scatter_add(numeric(n), pp$row, pp$score, 1, FALSE)
    ord <- ofst_store_key_order(pp$keys)
    keep <- ord[total[ord] != 0]
    pos <- integer(n); pos[keep] <- seq_along(keep)
    d <- data.table::as.data.table(lapply(pp$keys, `[`, keep))
    data.table::set(d, j = "score", value = total[keep])
    for (j in seq_along(lib_names)) {
      v <- numeric(length(keep)); s <- pp$lib == j
      v[pos[pp$row[s]]] <- pp$score[s]
      data.table::set(d, j = lib_names[j], value = v)
    }
    d
  })
  dt <- .ofst_store_assemble(pieces, x$seqlevels, x$key_columns)
  for (j in lib_names) { v <- dt[[j]]; if (.is_safe_int32(v)) data.table::set(dt, j = j, value = as.integer(v)) }
  dt[]
}

#' Write an in-RAM OFST store matrix as a new store
#'
#' The write half of the in-RAM mode: builds a fresh, compact store (keys sorted, one
#' contribution segment series per partition) from an \code{ofst_store_matrix}, e.g. after
#' editing it in RAM.
#' @param x an \code{ofst_store_matrix} from \code{ofst_store_read_matrix(as = "matrix")}.
#' @param store character, new store directory (must not exist or be empty).
#' @param max_contrib_rows numeric, rows per contribution segment, default 1e9.
#' @param partition_min_rows numeric, used for chromosomes added later, default 1e6.
#' @return invisible manifest.
#' @family ofst_store
#' @export
ofst_store_write_matrix <- function(x, store, max_contrib_rows = 1e9, partition_min_rows = 1e6) {
  store <- .ofst_store_path(store)
  if (!inherits(x, "ofst_store_matrix")) .ofst_store_abort("x must be an ofst_store_matrix.")
  if (dir.exists(store) && length(list.files(store, all.files = TRUE, no.. = TRUE)))
    .ofst_store_abort("'", store, "' exists and is not empty.")
  dir.create(store, recursive = TRUE, showWarnings = FALSE)
  m <- .ofst_store_new_manifest(x$key_columns, max_contrib_rows)
  m$partition_min_rows <- partition_min_rows
  m$seqlevels <- x$seqlevels
  m$generation <- 1L
  kc <- c("seq", x$key_columns)
  libs <- list(); segments <- list(); parts <- list(); chr_part <- list()
  for (p in names(x$parts)) {
    pp <- x$parts[[p]]
    dir.create(file.path(store, p), showWarnings = FALSE)
    n <- length(pp$keys[[1L]])
    o <- ofst_store_key_order(pp$keys)
    keys <- data.table::as.data.table(lapply(pp$keys, `[`, o))
    data.table::setnames(keys, kc)
    newid <- integer(n); newid[o] <- seq_len(n)          # rows renumbered in key order
    fst::write_fst(keys, file.path(store, p, "keys.fst"))
    fst::write_fst(data.frame(row = seq_len(n)), .ofst_store_gen_file(store, p, "order", 1L))
    row <- newid[pp$row]
    fst::write_fst(data.frame(score = ofst_store_scatter_add(numeric(n), row, pp$score, 1, FALSE)),
                   .ofst_store_gen_file(store, p, "total", 1L))
    lo <- order(pp$lib, row)
    cd <- data.table::data.table(row = row[lo], score = pp$score[lo], lib = pp$lib[lo])
    cnt <- cd[, .(nnz = .N, score_sum = sum(score)), by = lib]
    cur <- 0L; rows <- Inf; first <- 1
    for (i in seq_len(nrow(cnt))) {
      if (rows + cnt$nnz[i] > max_contrib_rows) {
        cur <- cur + 1L; rows <- 0
      }
      piece <- cd[first:(first + cnt$nnz[i] - 1), .(row, score)]
      .ofst_store_append(piece, .ofst_store_seg_file(store, p, cur), rows)
      libs[[length(libs) + 1L]] <- data.table::data.table(lib = x$libs[cnt$lib[i]], part = p, seg = cur, from = rows + 1,
                                                          to = rows + cnt$nnz[i], nnz = cnt$nnz[i], score_sum = cnt$score_sum[i])
      rows <- rows + cnt$nnz[i]; first <- first + cnt$nnz[i]
      segments[[paste(p, cur)]] <- data.table::data.table(part = p, seg = cur, n_rows = rows)
    }
    parts[[p]] <- data.table::data.table(part = p, n_keys = n, gen = 1L)
    chr_part[[p]] <- data.table::data.table(chr = x$seqlevels[unique(pp$keys[[1L]])], part = p)
  }
  m$parts <- data.table::rbindlist(c(list(m$parts), parts))
  m$chr_part <- data.table::rbindlist(c(list(m$chr_part), chr_part))
  m$segments <- data.table::rbindlist(c(list(m$segments), segments))
  l <- data.table::rbindlist(libs)
  if (nrow(l)) l[, `:=`(source = NA_character_, source_size = NA_real_, source_mtime = NA_character_,
                        state = "live", generation = 1L)]
  m$libs <- data.table::rbindlist(list(m$libs, l), use.names = TRUE)
  .ofst_store_commit(store, m)
  invisible(m)
}

#' @export
print.ofst_store_matrix <- function(x, ...) {
  nnz <- sum(vapply(x$parts, function(p) length(p$row), numeric(1)))
  keys <- sum(vapply(x$parts, function(p) length(p$keys[[1L]]), numeric(1)))
  cat("ofst_store_matrix:", .ofst_number(keys), "keys x", length(x$libs), "inputs,", .ofst_number(nnz),
      "non-zero cells in", length(x$parts), "partitions\n")
  invisible(x)
}

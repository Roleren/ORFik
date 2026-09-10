context("Global, reproducible OFST positional filtering")

filter_fixture <- function(n = 24L) {
  data.table::data.table(seqnames = rep(c("chr1", "chr2", "chr3"), length.out = n),
                         start = seq_len(n), strand = "+", cigar = "10M", score = 1L)
}

filter_plain <- function(x) {
  x <- data.table::copy(x)
  data.table::setattr(x, "removal_summary", NULL)
  for (column in names(x)) if (is.factor(x[[column]]))
    data.table::set(x, j = column, value = as.character(x[[column]]))
  data.table::setcolorder(x, sort(names(x)))
  data.table::setorderv(x, names(x))
  as.data.frame(x)
}

filter_files <- function(tables, ..., filter_input_summary = TRUE) {
  paths <- vapply(seq_along(tables), function(i) tempfile("filter-input-", fileext = ".ofst"), "")
  on.exit(unlink(paths))
  for (i in seq_along(tables)) fst::write_fst(tables[[i]], paths[i])
  suppressMessages(ORFik:::ofst_merge(paths, lib_names = paste0("lib", seq_along(paths)),
                                    filter_input_summary = filter_input_summary, ...))
}

test_that("ranking pools all libraries, strands and CIGARs before removal", {
  a <- filter_fixture(4)
  a[, score := c(1L, 4L, 2L, 20L)]
  b <- a[1:3]
  b[, `:=`(strand = "-", cigar = "5M100N5M", score = c(100L, 0L, 0L))]
  # Seven alignment keys: start 3 (score 2, two keys) must be removed before
  # start 2 (score 4), and start 1 has score 101 despite score 1 in file one.
  for (keep in c(FALSE, TRUE)) {
    result <- filter_files(list(a, b), keep_all_scores = keep,
                            filter_target_rows = 5, filter_chunk_rows = 4)
    expect_setequal(result$start, c(1L, 2L, 4L))
    expect_equal(nrow(result), 5L)
    expect_setequal(result$cigar, c("10M", "5M100N5M"))
    expect_setequal(as.character(result$strand), c("+", "-"))
    expect_equal(sum(result$score), 125)
    s <- attr(result, "removal_summary", exact = TRUE)
    expect_s3_class(s, "ofst_removal_summary")
    expect_equal(s$grouping, c("seqnames", "start"))
    expect_equal(c(s$rows_before, s$rows_after, s$rows_removed), c(7, 5, 2))
    expect_equal(c(s$positions_before, s$positions_removed), c(4, 1))
    expect_equal(c(s$score_before, s$score_after, s$score_removed), c(127, 125, 2))
    expect_equal(s$cutoff_score, 2)
    expect_equal(sum(s$by_chromosome$rows_removed), s$rows_removed)
    expect_equal(sum(s$by_chromosome$score_removed), s$score_removed)
    expect_equal(s$by_input$rows_removed, c(1, 1))
    expect_equal(s$by_input$score_removed, c(2, 0))
    if (keep) expect_equal(c(sum(result$lib1, na.rm = TRUE), sum(result$lib2, na.rm = TRUE)), c(25, 100))
  }
})

test_that("seeded ties cover the whole genome and do not depend on order or chunk size", {
  a <- filter_fixture(90L)
  b <- data.table::copy(a)
  b[, seqnames := factor(seqnames, levels = c("chr3", "unused", "chr2", "chr1"))]
  run <- function(tables, seed = 7L, chunk = 24L, ...) filter_files(tables,
    keep_all_scores = FALSE, filter_target_rows = 45, filter_seed = seed,
    filter_chunk_rows = chunk, ...)
  set.seed(531)
  rng <- .Random.seed
  first <- run(list(a, b))
  expect_identical(.Random.seed, rng)
  expect_equal(filter_plain(first), filter_plain(run(list(b[.N:1], a[.N:1]), chunk = 64)))
  expect_equal(filter_plain(first), filter_plain(run(list(a, b), dt_max_index_size = 5)))
  expect_false(identical(first$start, run(list(a, b), seed = 8L)$start))
  expect_equal(nrow(first), 45L)
  s <- attr(first, "removal_summary")
  expect_true(all(s$by_chromosome$rows_removed > 0))
  expect_true(all(s$by_chromosome$rows_after > 0))
  # Splitting the same inputs differently must not affect pooled selection.
  expect_equal(filter_plain(first), filter_plain(run(list(a[1:45], a[46:90], b))))
})

test_that("actual distinct rows, not input sum or chunk size, trigger filtering", {
  a <- filter_fixture(8)
  for (keep in c(FALSE, TRUE)) {
    x <- filter_files(list(a, a, a), keep_all_scores = keep,
                      filter_target_rows = 8, filter_chunk_rows = 8, allow_filtering = FALSE)
    expect_equal(nrow(x), 8L)
    expect_equal(sum(x$score), 24)
    expect_null(attr(x, "removal_summary", exact = TRUE))
    y <- filter_files(list(a, a), keep_all_scores = keep, dt_max_index_size = 9)
    expect_null(attr(y, "removal_summary", exact = TRUE))
  }
})

test_that("strict mode and real score ceilings give actionable errors", {
  a <- filter_fixture(6)
  a[, score := 1:6]
  expect_error(filter_files(list(a), filter_target_rows = 3, allow_filtering = FALSE),
               "6 distinct merged rows.*allow_filtering=FALSE")
  expect_error(filter_files(list(a), filter_target_rows = 2, max_filter_score = 3),
               "must remove at least 4.*permits only 3")
  x <- filter_files(list(a), filter_target_rows = 3, max_filter_score = 3)
  expect_equal(x$start, 4:6)
  expect_equal(attr(x, "removal_summary")$cutoff_score, 3)
})

test_that("whole positions can undershoot the target and zero target yields valid empty tables", {
  a <- filter_fixture(4)
  b <- a[1:2]
  b[, strand := "-"]
  a[, score := c(0L, 20L, 30L, 40L)]
  b[, score := c(0L, 20L)]
  for (keep in c(FALSE, TRUE)) {
    x <- filter_files(list(a, b), keep_all_scores = keep, filter_target_rows = 5)
    expect_equal(nrow(x), 4L)
    expect_false(1L %in% x$start)
    y <- filter_files(list(a, b, a[0]), keep_all_scores = keep, filter_target_rows = 0)
    expect_equal(nrow(y), 0L)
    expect_true(all(c("seqnames", "start", "strand", "cigar", "score") %in% names(y)))
    expect_true(is.factor(y$seqnames))
    s <- attr(y, "removal_summary")
    expect_equal(s$rows_removed, 6)
    expect_equal(s$score_after, 0)
    expect_equal(s$by_input$rows_after, c(0, 0, 0))
    expect_equal(s$by_input$raw_input_rows, c(4, 2, 0))
    if (keep) expect_true(all(paste0("lib", 1:3) %in% names(y)))
  }
})

test_that("NA scores are zero, fractional scores stay exact, and count overflow is avoided", {
  a <- filter_fixture(4)
  a[, score := c(NA_real_, 0.25, 0.5, 2^31)]
  x <- filter_files(list(a, a), filter_target_rows = 2)
  expect_equal(sort(x$score), c(1, 2^32))
  expect_type(x$score, "double")
  expect_equal(attr(x, "removal_summary")$score_removed, .5)
})

test_that("discarding CIGAR changes alignment counts but not positional pooling", {
  a <- filter_fixture(4)
  b <- data.table::copy(a)
  b[, cigar := "5M100N5M"]
  for (keep in c(FALSE, TRUE)) {
    x <- filter_files(list(a, b), keep_all_scores = keep, keepCigar = FALSE,
                      filter_target_rows = 3, filter_chunk_rows = 4)
    expect_false("cigar" %in% names(x))
    expect_equal(nrow(x), 3L)
    expect_equal(attr(x, "removal_summary")$rows_before, 4)
    expect_equal(sum(x$score), 6)
  }
})

test_that("invalid controls fail early and name the bad argument", {
  a <- filter_fixture(2)
  cases <- list(allow_filtering = list(NA, 1, c(TRUE, FALSE), matrix(TRUE)),
                 filter_target_rows = list(-1, 1.5, NA, Inf, 2^31 - 1, matrix(1)),
                 filter_seed = list(-1, 1.2, NA, 2^31),
                 max_filter_score = list(-1, NA, "3", matrix(1)),
                 filter_chunk_rows = list(0, 1.2, Inf, 2^31),
                 filter_tmpdir = list(NA_character_, "", 4, matrix(tempdir())))
  for (arg in names(cases)) for (value in cases[[arg]]) {
    args <- c(list(tables = list(a)), setNames(list(value), arg))
    expect_error(do.call(filter_files, args), arg)
  }
})

test_that("unsafe rescue inputs and impossible partition budgets fail clearly", {
  a <- filter_fixture(4)
  for (bad in c(-1, Inf, -Inf)) {
    x <- data.table::copy(a)
    x[, score := as.double(score)]
    x[1, score := bad]
    expect_error(filter_files(list(x), filter_target_rows = 2), "negative or infinite scores")
  }
  for (column in c("seqnames", "start")) {
    x <- data.table::copy(a)
    data.table::set(x, 1L, column, NA)
    expect_error(filter_files(list(x), filter_target_rows = 2), "missing chromosome/start")
  }
  paired <- data.table::copy(a)
  data.table::setnames(paired, c("start", "cigar"), c("start1", "cigar1"))
  paired[, `:=`(start2 = start1 + 20L, cigar2 = "10M")]
  expect_error(filter_files(list(paired), filter_target_rows = 2), "paired start1/start2")
  a[, `:=`(seqnames = "chr1", start = 1L, cigar = paste0(1:4, "M"))]
  expect_error(filter_files(list(a), filter_target_rows = 2, filter_chunk_rows = 2),
               "Increase filter_chunk_rows")
})

test_that("scratch cleanup is scoped on success and error and inputs are never modified", {
  parent <- tempfile("filter-scratch-test-")
  dir.create(parent)
  on.exit(unlink(parent, recursive = TRUE))
  path <- file.path(parent, "input.ofst")
  fst::write_fst(filter_fixture(8), path)
  checksum <- tools::md5sum(path)
  run <- function(...) suppressMessages(ORFik:::ofst_merge(path, filter_target_rows = 4,
    filter_tmpdir = parent, filter_chunk_rows = 4, ...))
  run()
  expect_identical(list.files(parent), "input.ofst")
  expect_error(run(allow_filtering = FALSE), "allow_filtering=FALSE")
  expect_identical(list.files(parent), "input.ofst")
  expect_identical(tools::md5sum(path), checksum)
  expect_error(filter_files(list(filter_fixture(4)), filter_target_rows = 2,
                             filter_tmpdir = file.path(parent, "missing")), "scratch directory.*does not exist")
})

test_that("internal entry point rescues original tables and chunkified library scores", {
  a <- ORFik:::.normalize_dt(filter_fixture(6))
  b <- data.table::copy(a)
  b[, strand := "-"]
  run <- function(...) suppressMessages(ORFik:::ofst_merge_internal(...,
    filter_target_rows = 6, filter_chunk_rows = 8))
  for (keep in c(FALSE, TRUE)) {
    actual <- run(list(data.table::copy(a), data.table::copy(b)), c("lib1", "lib2"), keep_all_scores = keep)
    expected <- filter_files(list(a, b), keep_all_scores = keep, filter_target_rows = 6, filter_chunk_rows = 8)
    expect_equal(filter_plain(actual), filter_plain(expected))
    expect_equal(attr(actual, "removal_summary")$input_kind, "in_memory_tables")
  }
  chunk1 <- data.table::copy(a)
  chunk1[, lib1 := score]
  chunk2 <- data.table::copy(b)
  chunk2[, lib2 := score]
  actual <- run(list(chunk1, chunk2), c("lib1", "lib2"), chunkified = TRUE)
  expected <- filter_files(list(a, b), filter_target_rows = 6, filter_chunk_rows = 8)
  expect_equal(filter_plain(actual), filter_plain(expected))
  expect_equal(attr(actual, "removal_summary")$input_kind, "chunk_library_contributions")
})

test_that("mergeLibs forwards controls and persists summaries without stale reports", {
  df <- ORFik::ORFik.template.experiment()
  df <- df[df$libtype == "RNA", ][1:2, ]
  out <- tempfile("filter-mergeLibs-")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE))
  paths <- file.path(out, c("a.ofst", "b.ofst"))
  for (p in paths) fst::write_fst(filter_fixture(12), p)
  run <- function(...) suppressMessages(ORFik::mergeLibs(df, out_dir = out, paths = paths,
    lib_names_full = c("lib1", "lib2"), filter_tmpdir = out, filter_chunk_rows = 8, ...))
  expect_error(run(filter_target_rows = 5, allow_filtering = FALSE), "allow_filtering=FALSE")
  expect_false(file.exists(file.path(out, "all.ofst")))
  run(filter_target_rows = 5, filter_seed = 43L, max_filter_score = 2)
  s <- readRDS(file.path(out, "all.ofst.removal_summary.rds"))
  expect_s3_class(s, "ofst_removal_summary")
  expect_equal(s$seed, 43L)
  expect_equal(s$rows_after, fst::metadata_fst(file.path(out, "all.ofst"))$nrOfRows)
  run()
  expect_null(readRDS(file.path(out, "all.ofst.removal_summary.rds")))
  expect_equal(fst::metadata_fst(file.path(out, "all.ofst"))$nrOfRows, 12)
})

test_that("randomized merges agree with an independent full-table weighted reference", {
  set.seed(1827)
  for (iteration in 1:10) {
    inputs <- lapply(1:3, function(i) data.table::data.table(
      seqnames = sample(c("chr1", "chr2", "chr3"), 35, TRUE),
      start = sample(1:12, 35, TRUE), strand = sample(c("+", "-"), 35, TRUE),
      cigar = sample(c("10M", "5M100N5M"), 35, TRUE), score = sample(0:5, 35, TRUE)))
    all <- data.table::rbindlist(inputs)
    merged <- all[, .(score = sum(score)), by = .(seqnames, start, strand, cigar)]
    positions <- merged[, .(pooled = sum(score), rows = .N), by = .(seqnames, start)]
    positions[, priority := ORFik:::.ofst_priority(seqnames, start, iteration)]
    data.table::setorderv(positions, c("pooled", "priority", "seqnames", "start"))
    target <- 30
    end <- which(cumsum(positions$rows) >= nrow(merged) - target)[1L]
    removed <- positions[seq_len(end), .(seqnames, start)]
    expected <- merged[!removed, on = c("seqnames", "start")]
    for (keep in c(FALSE, TRUE)) {
      actual <- filter_files(inputs, keep_all_scores = keep, filter_target_rows = target,
                              filter_seed = iteration, filter_chunk_rows = 16)
      scores <- actual$score
      summary <- attr(actual, "removal_summary")
      if (keep) {
        expect_equal(scores, rowSums(as.data.frame(actual[, paste0("lib", 1:3), with = FALSE]), na.rm = TRUE))
        actual[, (paste0("lib", 1:3)) := NULL]
      }
      expect_equal(filter_plain(actual), filter_plain(expected))
      expect_equal(summary$score_removed, sum(merged$score) - sum(expected$score))
      expect_equal(summary$rows_removed, nrow(merged) - nrow(expected))
      for (i in 1:3) {
        d <- inputs[[i]][!removed, on = c("seqnames", "start")]
        expect_equal(summary$by_input$score_after[i], sum(d$score))
      }
    }
  }
})

test_that("radix selection handles exact priority boundaries, duplicates and weights", {
  scratch <- tempfile("filter-radix-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE))
  values <- c(0, 1, 7, 8, 2^13 - 1, 2^13, 2^33, 2^43 - 1, 2^43, 2^53 - 1)
  h <- data.table::data.table(value = rep(values, 2), weight = rep(c(2, 1), each = length(values)))
  expected <- h[, .(weight = sum(weight)), by = value][order(value)]
  for (needed in c(1, 3, 4, 15, 29, 30)) for (budget in c(2, 100)) {
    paths <- vapply(split(h, rep(1:4, length.out = nrow(h))),
                     function(d) ORFik:::.ofst_spill(d, scratch), "")
    x <- ORFik:::.ofst_priority_cutoff(paths, needed, scratch, budget)
    index <- which(cumsum(expected$weight) >= needed)[1L]
    expect_equal(x$value, expected$value[index])
    expect_equal(x$below, sum(expected$weight[seq_len(index - 1L)]))
    unlink(list.files(scratch, full.names = TRUE))
  }
})

test_that("block reader uses double row offsets beyond the 32-bit boundary", {
  scratch <- tempfile("filter-offsets-")
  dir.create(scratch)
  on.exit(unlink(scratch, recursive = TRUE))
  seen <- list()
  sandbox <- new.env(parent = asNamespace("ORFik"))
  sandbox$.ofst_input_block <- function(path, from, to, keys) {
    args <- list(from = from, to = to)
    seen[[length(seen) + 1L]] <<- args
    filter_fixture(1)[, start := length(seen)]
  }
  stage <- ORFik:::.ofst_batch_runs
  environment(stage) <- sandbox
  suppressMessages(stage("mock-huge-input", 2^31 + 3, c("seqnames", "start", "strand", "cigar"),
                          ".source", 2^30, scratch))
  expect_equal(vapply(seen, `[[`, 0, "from"), c(1, 2^30 + 1, 2^31 + 1))
  expect_equal(vapply(seen, `[[`, 0, "to"), c(2^30, 2^31, 2^31 + 3))
  expect_false(anyNA(unlist(seen)))
})

test_that("extra alignment columns cannot shadow filtering helper variables", {
  a <- filter_fixture(12)
  for (column in c("keys", "merge_keys", "bucket", "index", "dt", "source_col", ".ofst_source", "drop"))
    a[, (column) := "metadata"]
  for (keep in c(FALSE, TRUE)) {
    plain <- filter_files(list(a, a), keep_all_scores = keep)
    expect_equal(nrow(plain), 12L)
    expect_equal(sum(plain$score), 24)
    filtered <- filter_files(list(a, a), keep_all_scores = keep,
                             filter_target_rows = 6, filter_chunk_rows = 4)
    expect_equal(nrow(filtered), 6L)
    expect_equal(sum(filtered$score), 12)
    expect_true(all(names(a) %in% names(filtered)))
    expect_true(all(filtered$bucket == "metadata"))
    expect_true(all(filtered$.ofst_source == "metadata"))
  }
})

test_that("messages distinguish rescue, no-op, filtering and scratch failures", {
  path <- tempfile("filter-messages-", fileext = ".ofst")
  on.exit(unlink(path))
  fst::write_fst(filter_fixture(6), path)
  messages <- character()
  x <- withCallingHandlers(ORFik:::ofst_merge(path, filter_target_rows = 3),
    message = function(m) { messages <<- c(messages, conditionMessage(m)); invokeRestart("muffleMessage") })
  expect_true(any(grepl("distinct merged rows = 6", messages, fixed = TRUE)))
  expect_true(any(grepl("Randomizing tied positions globally with seed 1", messages, fixed = TRUE)))
  expect_true(any(grepl("removed 3 positions / 3 alignment rows", messages, fixed = TRUE)))
  expect_true(any(grepl("Summary attached", messages, fixed = TRUE)))
  expect_equal(attr(x, "removal_summary")$positions_after, 3)
  expect_error(ORFik:::.ofst_read_part(paste0(path, "-missing")), "cannot read.*missing")
  expect_error(ORFik:::.ofst_spill(filter_fixture(1), paste0(path, "-missing")), "Check disk space and permissions")
  expect_error(ORFik:::.ofst_save_merge(filter_fixture(1), file.path(paste0(path, "-missing"), "out.ofst")),
               "existing output was not replaced")
})

test_that("filtering does not create an RNG state and scores cannot silently overflow", {
  path <- tempfile("filter-rng-", fileext = ".ofst")
  fst::write_fst(filter_fixture(4), path)
  on.exit(unlink(path), add = TRUE)
  existed <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  if (existed) saved <- get(".Random.seed", envir = .GlobalEnv)
  on.exit(if (existed) assign(".Random.seed", saved, envir = .GlobalEnv), add = TRUE)
  if (existed) rm(".Random.seed", envir = .GlobalEnv)
  suppressMessages(ORFik:::ofst_merge(path, filter_target_rows = 2))
  expect_false(exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
  a <- filter_fixture(4)
  data.table::set(a, j = "score", value = rep(1e308, nrow(a)))
  expect_error(filter_files(list(a, a), filter_target_rows = 2), "scores overflowed")
  expect_error(filter_files(list(a), filter_target_rows = 2), "aggregate input score overflowed")
})

test_that("native FST offset bridge reads real short blocks and dispatches long offsets", {
  path <- tempfile("filter-native-", fileext = ".ofst")
  on.exit(unlink(path))
  input <- filter_fixture(8)
  input[, seqnames := factor(seqnames)]
  fst::write_fst(input, path)
  actual <- ORFik:::.ofst_read_large_offset(path, 2, 5)
  expect_equal(filter_plain(actual), filter_plain(input[2:5]))
  actual <- ORFik:::.ofst_read_large_offset(path, 2, 5, columns = c("start", "cigar"))
  expect_equal(filter_plain(actual), filter_plain(input[2:5, .(start, cigar)]))
  expect_error(ORFik:::.ofst_read_large_offset(path, 2^31 + 1, 2^31 + 2),
               "Row selection.*range|row.*range", ignore.case = TRUE)
  sandbox <- new.env(parent = asNamespace("ORFik"))
  seen <- NULL
  sandbox$.ofst_read_large_offset <- function(path, from, to, columns) {
    seen <<- c(from, to)
    input[1]
  }
  read <- ORFik:::.ofst_read_part
  environment(read) <- sandbox
  read(path, from = 2^31 + 1, to = 2^31 + 5)
  expect_equal(seen, c(2^31 + 1, 2^31 + 5))
})

test_that("replicate and library groups match exact labels, not regular expressions", {
  out <- tempfile("filter-groups-")
  dir.create(out)
  on.exit(unlink(out, recursive = TRUE))
  paths <- file.path(out, paste0("input", 1:4, ".ofst"))
  for (i in 1:4) {
    d <- filter_fixture(4)
    d[, score := i]
    fst::write_fst(d, paths[i])
  }
  sandbox <- new.env(parent = asNamespace("ORFik"))
  sandbox$bamVarName <- function(...) c("RNA_A", "RNA_A+", "RNA_A", "RNA_A+")
  merge <- ORFik::mergeLibs
  environment(merge) <- sandbox
  df <- ORFik::ORFik.template.experiment()[1:4, ]
  for (mode in c("rep", "lib")) {
    suppressMessages(merge(df, out_dir = out, paths = paths, mode = mode,
      lib_names_full = paste0("lib", 1:4), filter_target_rows = 2))
    a <- readRDS(file.path(out, "RNA_A.ofst.removal_summary.rds"))
    b <- readRDS(file.path(out, "RNA_A+.ofst.removal_summary.rds"))
    expect_equal(a$by_input$library, c("lib1", "lib3"))
    expect_equal(b$by_input$library, c("lib2", "lib4"))
    expect_equal(a$score_after, 8)
    expect_equal(b$score_after, 12)
  }
})

test_that("internal validation and the deprecated ceiling provide useful diagnostics", {
  a <- ORFik:::.normalize_dt(filter_fixture(4))
  for (flag in c("keep_all_scores", "keepCigar", "sort", "chunkified")) {
    args <- c(list(dt_list = list(a), lib_names = "lib1"), setNames(list(NA), flag))
    expect_error(do.call(ORFik:::ofst_merge_internal, args), paste0(flag, " must"))
  }
  b <- data.table::copy(a)
  b[, lib1 := as.character(score)]
  expect_error(ORFik:::ofst_merge_internal(list(b), "lib1", chunkified = TRUE), "library 'lib1'.*numeric scores")
  expect_error(ORFik:::ofst_merge_internal(list(a), "lib1", max_filter_value = 1, max_filter_score = 1),
               "not both max_filter_score and legacy")
  expect_warning(suppressMessages(ORFik:::ofst_merge_internal(list(a), "lib1", filter_target_rows = 2,
                                                             max_filter_value = 1)), "deprecated")
})

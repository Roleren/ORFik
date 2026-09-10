context("Batch-first OFST merging")

batch_fixture <- function(n = 12L) data.table::data.table(
  seqnames = rep(c("chr1", "chr2"), length.out = n), start = seq_len(n),
  strand = "+", cigar = "10M", score = 1L)

batch_plain <- function(x) {
  x <- data.table::copy(x)
  data.table::setattr(x, "removal_summary", NULL)
  for (column in names(x)) if (is.factor(x[[column]]))
    data.table::set(x, j = column, value = as.character(x[[column]]))
  data.table::setcolorder(x, sort(names(x)))
  data.table::setorderv(x, names(x))
  as.data.frame(x)
}

batch_files <- function(tables, fun) {
  root <- tempfile("batch-test-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  paths <- file.path(root, paste0(seq_along(tables), ".ofst"))
  for (i in seq_along(tables)) fst::write_fst(tables[[i]], paths[i])
  fun(paths, root)
}

test_that("pooled batches collapse across inputs before writing and compact before the last input", {
  a <- batch_fixture()
  batch_files(rep(list(a), 12), function(paths, root) {
    spill <- ORFik:::.ofst_spill
    read <- ORFik:::.ofst_input_block
    seen <- 0L
    before_last <- 0L
    largest_spill <- 0
    peak_cached_rows <- 0
    testthat::local_mocked_bindings(
      .ofst_input_block = function(path, from, to, keys) {
        seen <<- max(seen, match(path, paths))
        read(path, from, to, keys)
      },
      .ofst_spill = function(dt, scratch, prefix = "part-") {
        if (prefix == "batch-") {
          expect_false(".ofst_source" %in% names(dt))
          largest_spill <<- max(largest_spill, nrow(dt))
          if (seen < length(paths)) before_last <<- before_last + 1L
        }
        path <- spill(dt, scratch, prefix)
        files <- list.files(scratch, pattern = "^batch-.*fst$", full.names = TRUE)
        peak_cached_rows <<- max(peak_cached_rows,
          sum(vapply(files, function(p) as.double(fst::metadata_fst(p)$nrOfRows), 0)))
        path
      }, .package = "ORFik")
    result <- suppressMessages(ORFik:::ofst_merge(paths, keep_all_scores = FALSE,
      filter_target_rows = 12, filter_chunk_rows = 36, filter_tmpdir = root))
    expect_equal(nrow(result), 12L)
    expect_equal(result$score, rep(12L, 12))
    expect_null(attr(result, "removal_summary"))
    expect_gt(before_last, 0L)
    expect_equal(largest_spill, 12)
    expect_lte(peak_cached_rows, 48)
    expect_setequal(list.files(root), basename(paths))
  })
})

test_that("adaptive partition trees agree with a direct merge and have no loss without filtering", {
  a <- batch_fixture(80)
  b <- data.table::copy(a)
  b[, `:=`(start = start + 30L, strand = "-", score = 3L)]
  batch_files(list(a, b, a[1:17], b[55:80]), function(paths, root) {
    expected <- data.table::rbindlist(list(a, b, a[1:17], b[55:80]))[
      , .(score = sum(score)), by = .(seqnames, start, strand, cigar)]
    for (budget in c(7, 29, 100)) {
      actual <- suppressMessages(ORFik:::ofst_merge(paths, keep_all_scores = FALSE,
        filter_target_rows = 160, allow_filtering = FALSE, filter_chunk_rows = budget, filter_tmpdir = root))
      expect_equal(batch_plain(actual), batch_plain(expected))
      expect_null(attr(actual, "removal_summary"))
    }
  })
})

test_that("optional detailed reports agree across batches and do not change global filtering", {
  a <- batch_fixture(28)
  a[, score := (start %% 5L) + 1L]
  b <- data.table::rbindlist(list(a, a[1:8], a[1:8]))
  batch_files(list(a, b, a[0]), function(paths, root) {
    run <- function(detail, budget) suppressMessages(ORFik:::ofst_merge(paths,
      keep_all_scores = FALSE, filter_target_rows = 13, filter_chunk_rows = budget,
      filter_input_summary = detail, filter_tmpdir = root))
    fast <- run(FALSE, 7)
    detailed <- run(TRUE, 11)
    expect_equal(batch_plain(fast), batch_plain(detailed))
    s <- attr(fast, "removal_summary")
    d <- attr(detailed, "removal_summary")
    expect_false(s$input_summary_computed)
    expect_true(all(is.na(s$by_input$rows_removed)))
    expect_true(d$input_summary_computed)
    expect_equal(d$by_input$rows_before, c(28, 28, 0))
    expect_equal(d$by_input$score_before, c(sum(a$score), sum(b$score), 0))
    expect_equal(d$by_input$score_after, c(sum(a[start %in% fast$start]$score),
      sum(b[start %in% fast$start]$score), 0))
    expect_equal(s$score_removed, d$score_removed)
    expect_equal(s$by_chromosome, d$by_chromosome)
    expect_identical(s$merge_strategy, "bounded_batch_first")
    expect_equal(s$schema_version, 2L)
  })
})

test_that("mergeLibs uses output storage immediately unless explicitly overridden", {
  batch_files(list(batch_fixture(), batch_fixture()), function(paths, root) {
    df <- ORFik::ORFik.template.experiment()[1:2, ]
    out <- file.path(root, "output")
    spill <- ORFik:::.ofst_spill
    testthat::local_mocked_bindings(.ofst_spill = function(dt, scratch, prefix = "part-") {
      expect_equal(dirname(scratch), normalizePath(out))
      spill(dt, scratch, prefix)
    }, .package = "ORFik")
    suppressMessages(ORFik::mergeLibs(df, out_dir = out, paths = paths,
      lib_names_full = c("a", "b"), keep_all_scores = FALSE,
      filter_target_rows = 6, filter_chunk_rows = 8))
    expect_setequal(list.files(out), c("all.ofst", "all.ofst.removal_summary.rds"))
    expect_false(readRDS(file.path(out, "all.ofst.removal_summary.rds"))$scratch_fallback_used)
  })
})

test_that("batch row offsets cover each original row once even at study boundaries", {
  batch_files(list(batch_fixture(9), batch_fixture(5), batch_fixture(17)), function(paths, root) {
    read <- ORFik:::.ofst_input_block
    seen <- vector("list", 3)
    testthat::local_mocked_bindings(.ofst_input_block = function(path, from, to, keys) {
      i <- match(path, paths)
      seen[[i]] <<- c(seen[[i]], seq(from, to))
      expect_lte(to - from + 1, 8)
      read(path, from, to, keys)
    }, .package = "ORFik")
    suppressMessages(ORFik:::ofst_merge(paths, keep_all_scores = FALSE,
      filter_target_rows = 17, filter_chunk_rows = 8, filter_tmpdir = root))
    expect_equal(seen, lapply(c(9, 5, 17), seq_len))
  })
})

test_that("invalid detail controls fail before touching input data", {
  for (bad in list(NA, 1, c(TRUE, FALSE), matrix(TRUE)))
    expect_error(ORFik:::.ofst_filter_controls(TRUE, 10, 1L, Inf, 10, tempdir(),
      filter_input_summary = bad), "filter_input_summary must be TRUE or FALSE")
})

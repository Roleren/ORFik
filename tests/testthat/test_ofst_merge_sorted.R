sorted_merge_fixture <- function(seed, chrs = c("1", "2", "10", "X", "MT"), n = 400, width = FALSE) {
  set.seed(seed)
  d <- data.table::data.table(seqnames = factor(sample(chrs, n, TRUE), levels = chrs),
                              start = sample.int(300L, n, TRUE),
                              strand = factor(sample(c("+", "-", "*"), n, TRUE), levels = c("+", "-", "*")),
                              size = sample(26:30, n, TRUE), score = sample.int(5L, n, TRUE))
  # same layout as merged study pools: chromosomes contiguous in factor-code order, then keys
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand, size)]
  data.table::setkey(d, NULL)
  if (width) { d[, width := 1L]; data.table::setcolorder(d, c("seqnames", "start", "width")) }
  d
}

sorted_merge_write <- function(tables) {
  paths <- vapply(seq_along(tables), function(i) tempfile(fileext = ".ofst"), "")
  for (i in seq_along(tables)) fst::write_fst(tables[[i]], paths[i])
  paths
}

general_merge <- function(paths) {
  d <- .ofst_merge_lossless(.read_fst_list(paths), lib_names = paste0("l", seq_along(paths)),
                            keep_all_scores = FALSE)
  .ofst_finalize(d)
}

test_that("sorted chromosome-unit merge equals the general lossless merge", {
  paths <- sorted_merge_write(lapply(1:6, sorted_merge_fixture))
  expected <- general_merge(paths)
  keys <- c("seqnames", "start", "strand", "size")
  for (unit_rows in c(1e9, 150, 1)) for (workers in c(1L, 2L)) {
    got <- .ofst_finalize(.ofst_merge_sorted(paths, keys, unit_rows = unit_rows, workers = workers, verbose = FALSE))
    expect_identical(got[], expected[])
  }
  # public entry point takes the sorted path for pooled output
  expect_message(public <- ofst_merge(paths, keep_all_scores = FALSE), "sorted merge")
  expect_identical(public[], expected[])
  # chromosome labels sort as characters ("1" < "10" < "2"), as .ofst_sort
  expect_identical(unique(as.character(public$seqnames)), c("1", "10", "2", "MT", "X"))
})

test_that("missing scores count as zero and differing chromosome sets merge", {
  a <- sorted_merge_fixture(1, chrs = c("1", "2"))
  b <- sorted_merge_fixture(2, chrs = c("2", "3", "Y"))
  b[c(1, 5), score := NA_integer_]
  paths <- sorted_merge_write(list(a, b))
  got <- ofst_merge(paths, keep_all_scores = FALSE)
  expect_identical(got[], general_merge(paths)[])
  expect_equal(sum(got$score), sum(a$score) + sum(b$score, na.rm = TRUE))
})

test_that("an all-1 width column is dropped when only some inputs have it", {
  a <- sorted_merge_fixture(3)
  b <- sorted_merge_fixture(4, width = TRUE)
  paths <- sorted_merge_write(list(a, b))
  b_plain <- data.table::copy(b)[, width := NULL]
  expected <- general_merge(sorted_merge_write(list(a, b_plain)))
  got <- ofst_merge(paths, keep_all_scores = FALSE)
  expect_identical(got[], expected[])
  expect_false("width" %in% names(got))
  # width other than 1 cannot be dropped
  b[1, width := 2L]
  expect_error(ofst_merge(sorted_merge_write(list(a, b)), keep_all_scores = FALSE), "width column")
  # strict schema rules stay for per-library output
  expect_error(ofst_merge(paths, keep_all_scores = TRUE), "identical columns")
})

test_that("inputs with a chromosome split over several blocks use the general merge", {
  a <- sorted_merge_fixture(5)
  b <- sorted_merge_fixture(6)
  b <- b[c(seq(2, .N, 2), seq(1, .N, 2))]          # chromosomes interleaved
  paths <- sorted_merge_write(list(a, b))
  expect_null(.ofst_merge_sorted(paths, c("seqnames", "start", "strand", "size"), verbose = FALSE))
  expect_message(got <- ofst_merge(paths, keep_all_scores = FALSE), "general merge")
  expect_identical(got[], general_merge(paths)[])
})

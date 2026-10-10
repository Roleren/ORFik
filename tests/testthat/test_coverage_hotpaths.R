library(ORFik)

fst_hotpath_fixture <- function() {
  root <- tempfile("fst-hotpath-")
  dir.create(root)
  index <- data.table::rbindlist(lapply(c("1", "2"), function(chr) {
    data.table::rbindlist(lapply(c(1L, 11L), function(begin) {
      paths <- file.path(root, paste0(chr, "_", begin, "_", c("plus", "minus"), ".fst"))
      offset <- if (chr == "1") 0L else 100L
      fst::write_fst(data.frame(A = offset + begin:(begin + 9L), B = -begin:-(begin + 9L)), paths[1])
      fst::write_fst(data.frame(A = offset + 1000L + begin:(begin + 9L), B = 2L * begin:(begin + 9L)), paths[2])
      data.table::data.table(chr = chr, start = begin, end = begin + 9L,
        file_forward = paths[1], file_reverse = paths[2])
    }))
  }))
  path <- file.path(root, "coverage_index.fst")
  fst::write_fst(index, path)
  list(root = root, index = path)
}

test_that("FST page batching preserves exon, transcript and requested column order", {
  f <- fst_hotpath_fixture(); on.exit(unlink(f$root, recursive = TRUE))
  grl <- GRangesList(plus = GRanges("1", IRanges(c(12, 3), c(14, 5)), "+"),
    minus = GRanges("2", IRanges(c(15, 2), c(17, 4)), "-"),
    crossing = GRanges("1", IRanges(8, 13), "+"),
    duplicate = GRanges("1", IRanges(c(3, 3), c(5, 5)), "+"))
  ans <- coverageByTranscriptFST(grl, f$index)
  expect_identical(names(ans), names(grl))
  expect_equal(ans$plus$A, c(12:14, 3:5))
  expect_equal(ans$minus$A, c(1117:1115, 1104:1102))
  expect_equal(ans$crossing$A, 8:13)
  expect_equal(ans$duplicate$A, rep(3:5, 2L))
  cols <- coverageByTranscriptFST(grl, f$index, c("B", "A"))
  expect_identical(names(cols$plus), c("B", "A"))
  expect_equal(cols$plus$B, -ans$plus$A)
  expect_equal(coverageByTranscriptFST(grl, f$index, "B")$minus$B,
               2L * c(17:15, 4:2))
  # Identical slices must remain independent data.tables for callers using :=.
  same <- coverageByTranscriptFST(rep(grl[1], 2L), f$index)
  same[[1]][, A := 0L]
  expect_equal(same[[2]]$A, c(12:14, 3:5))
  empty <- coverageByTranscriptFST(GRangesList(empty = GRanges()), f$index)
  expect_identical(names(empty), "empty")
  expect_equal(nrow(empty[[1]]), 0L)
  # The chromosome key prevents matching a range to another chromosome's page.
  mixed <- GRangesList(tx = GRanges(c("1", "2"), IRanges(c(3, 3), c(5, 5)), "+"))
  expect_equal(coverageByTranscriptFST(mixed, f$index)[[1]]$A, c(3:5, 103:105))
})

test_that("FST indexes are reread and relocated page paths are supported", {
  f <- fst_hotpath_fixture(); on.exit(unlink(f$root, recursive = TRUE))
  grl <- GRangesList(tx = GRanges("1", IRanges(3, 5), "+"))
  idx <- fst::read_fst(f$index)
  idx$file_forward <- file.path("/old/location", basename(idx$file_forward))
  idx$file_reverse <- file.path("/old/location", basename(idx$file_reverse))
  fst::write_fst(idx, f$index)
  expect_equal(coverageByTranscriptFST(grl, f$index)[[1]]$A, 3:5)
  page <- file.path(f$root, basename(idx$file_forward[1]))
  d <- fst::read_fst(page); d$A <- d$A + 40L
  fst::write_fst(d, page)
  expect_equal(coverageByTranscriptFST(grl, f$index)[[1]]$A, 43:45)
  expect_error(coverageByTranscriptFST(grl, paste0(f$index, "missing")), "index")
  expect_error(coverageByTranscriptFST(grl, f$index, "missing"))
  expect_error(coverageByTranscriptFST(GRangesList(GRanges("3", IRanges(1, 2))), f$index))
})

test_that("coverage tables preserve positions, integer types, frames and gene factors", {
  cov <- RleList(z = Rle(c(0L, 2L, 2L, 0L, 4L)), a = Rle(c(3L, 0L)),
                 empty = Rle(integer()))
  dt <- ORFik:::coverage_to_dt(cov, drop.zero.dt = TRUE, withFrames = TRUE, fraction = 28L)
  expect_identical(dt$count, c(2L, 2L, 4L, 3L))
  expect_identical(dt$position, c(2L, 3L, 5L, 1L))
  expect_identical(dt$frame, c(1, 2, 1, 0))
  expect_identical(dt$genes, factor(c("z", "z", "z", "a")))
  expect_identical(dt$fraction, rep(28L, 4L))
  num <- ORFik:::coverage_to_dt(cov, keep.names = FALSE, drop.zero.dt = TRUE)
  expect_identical(num$genes, c(1L, 1L, 1L, 2L))
  full <- ORFik:::coverage_to_dt(cov, withFrames = TRUE)
  expect_identical(full$position, c(1:5, 1:2))
  expect_identical(full$count, c(0L, 2L, 2L, 0L, 4L, 3L, 0L))
  fractional <- cov / 3
  expect_equal(ORFik:::coverage_to_dt(fractional, drop.zero.dt = TRUE)$count, c(2, 2, 4, 3) / 3)
  expect_warning(zero <- ORFik:::coverage_to_dt(RleList(tx = Rle(0L, 20L)),
    drop.zero.dt = TRUE, withFrames = TRUE), "No coverage")
  expect_identical(names(zero), c("count", "genes", "position", "frame"))
  expect_equal(nrow(zero), 0L)
})

test_that("dense and sparse coverage tables agree with independently expanded profiles", {
  set.seed(119L)
  for (i in seq_len(15L)) {
    values <- lapply(seq_len(5L), function(j) sample(c(0L, 0L, 2L, 4L),
      sample.int(50L, 1L), replace = TRUE))
    cov <- RleList(lapply(values, Rle)); names(cov) <- paste0("tx", seq_along(values))
    full <- ORFik:::coverage_to_dt(cov, withFrames = TRUE)
    expect_identical(full$count, unlist(values, use.names = FALSE))
    expect_identical(full$position, unlist(lapply(values, seq_along), use.names = FALSE))
    sparse <- ORFik:::coverage_to_dt(cov, drop.zero.dt = TRUE, withFrames = TRUE)
    positive <- full$count > 0L
    expect_identical(sparse$count, full$count[positive])
    expect_identical(sparse$position, full$position[positive])
    expect_identical(sparse$frame, full$frame[positive])
  }
})

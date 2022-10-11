context("covRLE")
library(ORFik)

seqlengths <- as.integer(c(200, 300))
names(seqlengths) <- c("chr1", "chr2")
gr <- GRanges(seqnames = c("chr1", "chr1", "chr2", "chr2"),
              ranges = IRanges(start = c(10, 50, 100, 150), end = c(40, 80, 129, 179)),
              strand = c("+", "+", "-", "-"), seqlengths = seqlengths)
ga1 <- GAlignments(seqnames = c("chr1"), pos = as.integer(1), cigar = "1M1S",
                  strand = factor(c("+"), levels = c("+", "-", "*")))
ga2 <- GAlignments(seqnames = c("chr2"), pos = as.integer(170), cigar = "1M1S",
                  strand = factor(c("-"), levels = c("+", "-", "*")))
ga <- suppressWarnings(c(ga1, ga2))
seqlengths(ga) <- seqlengths

tx <- resize(gr, 9)
tx <- sortPerGroup(groupGRangesBy(tx, paste("tx", rep(c(1,2), each = 2))))


covRle_both <- covRleFromGR(gr)
covRle_pluss <- covRleFromGR(gr, ignore.strand = TRUE)
covRle_empty <- covRleFromGR(GRanges(seqlengths = seqlengths), ignore.strand = TRUE)
covRle_ga <- covRleFromGR(ga)
covRleList_both <- covRleList(list(covRle_both, covRle_both), fraction = c(28,29))
covRleList_pluss <- covRleList(list(covRle_pluss, covRle_pluss), fraction = c(28,29))

test_that("covRLE is created properly", {
  expect_is(covRle_pluss, "covRle")
  expect_is(covRle_both, "covRle")
  expect_equal(lengths(covRle_both), seqlengths)
  expect_equal(lengths(covRle_pluss), seqlengths)
  expect_equal(length(covRle_both), 2)
  expect_equal(length(covRle_pluss), 2)
  expect_equal(seqlevels(covRle_both), seqlevels(gr))
  expect_equal(seqlevels(covRle_pluss), seqlevels(gr))
  expect_equal(strandMode(covRle_both), 1)
  expect_equal(strandMode(covRle_pluss), 0)
})

test_that("covRLE fails when it should", {
  gr_bad <- gr
  seqlengths(gr_bad) <- NA
  expect_error(covRleFromGR(gr_bad))
  expect_error(covRle(f(covRle_both), r(covRle_both)[[1]]))
})

test_that("covRLE works in coveragePerTiling", {
  dt <- coveragePerTiling(grl = tx, covRle_both, as.data.table = TRUE)
  expect_is(dt, "data.table")
  expect_equal(nrow(dt), 36)
  expect_equal(sum(dt$count), 36)
  dt_pluss <- coveragePerTiling(grl = tx, covRle_pluss, as.data.table = TRUE)
  expect_equal(nrow(dt_pluss), 36)
  expect_equal(sum(dt_pluss$count), 36)
})

test_that("covRLEList is created properly", {
  expect_is(covRleList_pluss, "covRleList")
  expect_is(covRleList_both, "covRleList")
  expect_equal(lengths(covRleList_both), seqlengths)
  expect_equal(lengths(covRleList_pluss), seqlengths)
  expect_equal(length(covRleList_both), 2)
  expect_equal(length(covRleList_pluss), 2)
  expect_equal(seqlevels(covRleList_both), seqlevels(gr))
  expect_equal(seqlevels(covRleList_pluss), seqlevels(gr))
  expect_equal(strandMode(covRleList_both), 1)
  expect_equal(strandMode(covRleList_pluss), 0)
})

test_that("covRleList fails when it should", {
  expect_error(covRleList(list(RleList())))
  expect_error(covRleList(list(covRle_both, RleList())))
})

test_that("covRLE works in windowPerReadLength", {
  dt <- windowPerReadLength(tx, reads = covRleList_both, upstream = -1, downstream = 6)
  expect_is(dt, "data.table")
  expect_equal(nrow(dt), 12)
  expect_equal(sum(dt$score), 4)
})

test_that("covRLE works in regionPerReadLength", {
  dt <- regionPerReadLength(tx, reads = covRleList_both)
  expect_is(dt, "data.table")
  expect_equal(nrow(dt), 36)
  expect_equal(sum(dt$score), 4)
})

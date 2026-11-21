context("coveragePerTiling keep names (issue #100)")
library(ORFik)

test_that("as.data.table keeps names as factor with keep.names = 'factor'", {
  grl <- GRangesList(
    txB = GRanges("1", IRanges(c(1, 10), c(5, 15)), "+"),
    txA = GRanges("1", IRanges(30, 33), "+"),
    txC = GRanges("1", IRanges(c(40, 50), c(44, 55)), "-"))
  reads <- GRanges("1", IRanges(c(2, 3, 12, 31, 41, 52), width = 1),
                   c("+", "+", "+", "+", "-", "-"), seqlengths = c("1" = 60))
  # default (TRUE) keeps the integer gene index when zeros are kept
  dt_int <- coveragePerTiling(grl, reads, as.data.table = TRUE)
  expect_true(is.integer(dt_int$genes))
  # Opt-in: names as factor, in input order, positions unchanged
  dt <- coveragePerTiling(grl, reads, as.data.table = TRUE,
                          keep.names = "factor", withFrames = TRUE)
  expect_true(is.factor(dt$genes))
  expect_equal(levels(dt$genes), names(grl))
  expect_equal(as.integer(dt$genes), dt_int$genes)
  expect_equal(dt$position, dt_int$position)
  expect_equal(dt$count, dt_int$count)
  expect_equal(as.vector(table(dt$genes)[names(grl)]),
               as.vector(widthPerGroup(grl, FALSE)))
  expect_equal(dt[genes == "txA", sum(count)], 1)
  # same windows size path
  grl_eq <- grl[c(1, 3)]
  dt_eq <- coveragePerTiling(grl_eq, reads, as.data.table = TRUE,
                             keep.names = "factor")
  expect_equal(levels(dt_eq$genes), names(grl_eq))
  expect_equal(dt_eq$position, rep(1:11, 2))
  # drop.zero.dt path also returns names
  dt_sparse <- coveragePerTiling(grl, reads, as.data.table = TRUE,
                                 keep.names = "factor", drop.zero.dt = TRUE)
  expect_true(is.factor(dt_sparse$genes))
  expect_equal(sort(unique(as.character(dt_sparse$genes))), sort(names(grl)))
  # Rle output keeps names, FALSE still drops them
  expect_equal(names(coveragePerTiling(grl, reads, keep.names = "factor")),
               names(grl))
  expect_null(names(coveragePerTiling(grl, reads, keep.names = FALSE)))
  # duplicated names: positions still per group
  grl_dup <- grl; names(grl_dup) <- c("g", "g", "h")
  dt_dup <- coveragePerTiling(grl_dup, reads, as.data.table = TRUE,
                              keep.names = "factor")
  expect_equal(levels(dt_dup$genes), c("g", "h"))
  expect_equal(dt_dup$position, dt_int$position)
  expect_error(coveragePerTiling(grl, reads, keep.names = "yes"))
})

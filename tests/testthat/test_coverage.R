context("Coverage Helpers")
library(ORFik)

ORF1 <- GRanges("1", IRanges(21, 49), "+")
ORF2 <- GRanges("1", IRanges(30, 49), "+")
grl <- GRangesList("tx1_1" = ORF1, "tx1_2" = ORF2)
tx <- resize(resize(grl[1], width = 50), width = 70, fix = "end")
names(tx) <- "tx1"
footprintsGood <- GRanges("1", IRanges(seq.int(21, 49, 3), width = 1), "+")
footprintsGood$size <- 29
footprintsBad <- GRanges()
footprintsMiss <- GRanges("1", IRanges(500, width = 1), "+")
footprintsMiss$size <- 29

ga1 <- GAlignments(seqnames = c("1"), pos = as.integer(22), cigar = "1M1S",
                   strand = factor(c("+"), levels = c("+", "-", "*")))
ga2 <- GAlignments(seqnames = c("2"), pos = as.integer(170), cigar = "1M1S",
                   strand = factor(c("-"), levels = c("+", "-", "*")))
ga <- suppressWarnings(c(ga1, ga2))

test_that("coveragePerTiling: Input checks", {
  rle <- coveragePerTiling(grl, footprintsGood)
  expect_equal(as.integer(sum(runValue(rle))), c(10, 7))
  rle <- coveragePerTiling(grl, ga)
  expect_equal(as.integer(sum(runValue(rle))), c(1, 0))
})

test_that("coverageScorings works as intended", {
  # no reads hit
  coverage <- coveragePerTiling(grl, footprintsGood, is.sorted = TRUE,
                                as.data.table = TRUE)
  expect_is(coverage, "data.table")
  coverage2 <- coveragePerTiling(grl, footprintsGood, is.sorted = TRUE,
                                as.data.table = TRUE, drop.zero.dt = TRUE,
                                fraction = 28)
  expect_equal(nrow(coverage[count > 0,]), nrow(coverage2))
  expect_equal(unique(coverage2$fraction), 28)

  # fracPos
  dt <- coverageScorings(coverage, "fracPos")
  expect_equal(sum(dt$score), length(grl))

  # Sum
  dt <- coverageScorings(coverage, "sum")
  expect_equal(sum(dt$score), sum(countOverlaps(grl, footprintsGood)))

  # log2 sum
  dt <- coverageScorings(coverage, "log2sum")
  expect_equal(sum(dt$score[is.finite(dt$score)]),
               7)

  # Mean
  dt <- coverageScorings(coverage, "mean")
  expect_equal(sum(dt$score), 10)

  # Zscore
  dt <- coverageScorings(coverage, "zscore")
  expect_equal(round(sum(dt$score), 2), -0.11)

  # Transcript Normalized
  dt <- coverageScorings(coverage, "transcriptNormalized")
  expect_equal(round(sum(dt$score), 2), length(grl))

})

test_that("windowPerReadLength works as intended", {
  # per group coverage
  grltest <- windowPerReadLength(grl, tx, footprintsGood,
                                 scoring = "fracPos")
  expect_is(grltest, "data.table")
  expect_equal(nrow(grltest), 52)
  expect_equal(grltest$fraction[1], 29)
  expect_equal(c(min(grltest$position), max(grltest$position)), c(-5, 20))
  expect_equal(round(grltest$score[6], 3) , 0.143)
  # meta coverage
  grltest <- windowPerReadLength(grl, tx, footprintsGood)
  expect_is(grltest, "data.table")
  expect_equal(nrow(grltest), 26)
  expect_equal(round(grltest$score[6], 3) , 0.268)

  # - strand
  strand(grl) <- "-"
  strand(tx) <- "-"
  strand(footprintsGood) <- "-"
  grltest <- windowPerReadLength(grl, tx, footprintsGood,
                                 scoring = "fracPos")
  # Optimization works
  expect_equal(windowPerReadLength(grl, tx, footprintsGood,
                                   scoring = "transcriptNormalized"),
               windowPerReadLength(grl, tx, footprintsGood,
                                   scoring = "transcriptNormalized",
                                   drop.zero.dt = TRUE, append.zeroes = TRUE))
  expect_equal(windowPerReadLength(grl, tx, footprintsBad,
                                   scoring = "transcriptNormalized"),
               windowPerReadLength(grl, tx, footprintsBad,
                                   scoring = "transcriptNormalized",
                                   drop.zero.dt = TRUE, append.zeroes = TRUE))
  expect_equal(windowPerReadLength(grl, tx, footprintsMiss,
                                   scoring = "transcriptNormalized"),
               suppressWarnings(windowPerReadLength(grl, tx, footprintsMiss,
                                   scoring = "transcriptNormalized",
                                   drop.zero.dt = TRUE, append.zeroes = TRUE)))

})

test_that("windowPerReadLength works as intended strange cases", {
  # no reads hit
  grltest <- windowPerReadLength(grl, tx, footprintsBad,
                                 scoring = "fracPos")
  expect_is(grltest, "data.table")
  expect_equal(nrow(grltest), 0)
  # no grl
  grltest <- windowPerReadLength(GRangesList(), tx, footprintsGood)
  expect_is(grltest, "data.table")
  expect_equal(nrow(grltest), 0)
})

test_that("regionPerReadLength works as intended", {
  # Per frame
  grltest <- regionPerReadLength(grl, footprintsGood, scoring = "frameSumPerLG",
                                 drop.zero.dt = FALSE)
  expect_is(grltest, "data.table")
  expect_equal(nrow(grltest), 6)
  expect_equal(grltest$score[1], 10)
  expect_equal(grltest$score[4], 7)
  grltest <- regionPerReadLength(grl, footprintsGood, scoring = "frameSumPerL",
                                 drop.zero.dt = TRUE)
  expect_equal(nrow(grltest), 1)
  expect_equal(grltest$score[1], 17)
  grltest <- regionPerReadLength(grl, footprintsBad, scoring = "frameSumPerL",
                                 drop.zero.dt = TRUE)
  suppressWarnings(grltest <- regionPerReadLength(grl, footprintsMiss, scoring = "frameSumPerL",
                                 drop.zero.dt = TRUE))
  expect_is(grltest, "data.table")
  expect_equal(nrow(grltest), 0)
  expect_equal(ncol(grltest), 3)
})

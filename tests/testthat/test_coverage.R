context("GRanges Helpers")
library(ORFik)

ORF1 <- GRanges("1", IRanges(21, 49), "+")
ORF2 <- GRanges("1", IRanges(30, 49), "+")
grl <- GRangesList("tx1_1" = ORF1, "tx1_2" = ORF2)
tx <- resize(resize(grl[1], width = 50), width = 70, fix = "end")
names(tx) <- "tx1"
footprintsGood <- GRanges("1", IRanges(seq.int(21, 49, 3), width = 1), "+")
footprintsGood$score <- 29
footprintsBad <- GRanges()

test_that("windowPerReadLength works as intended", {
  # per group coverage
  grltest <- windowPerReadLength(grl, tx, footprintsGood,
                                         scoring = "fracPos")
  expect_is(grltest, "data.table")
  expect_equal(nrow(grltest), 52)
  expect_equal(grltest$fraction[1], 29)
  expect_equal(c(min(grltest$position), max(grltest$position)), c(-5,20))
  expect_equal(round(grltest$score[6], 3) , 0.143)
  # meta coverage
  grltest <- windowPerReadLength(grl, tx, footprintsGood)
  expect_is(grltest, "data.table")
  expect_equal(nrow(grltest), 26)
  expect_equal(round(grltest$score[6], 3) , 0.268)
})

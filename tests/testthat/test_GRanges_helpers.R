library(ORFik)
context("GRanges Helpers")



ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                     strand = Rle(strand(rep("+", 3))))

ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(20, 30, 40), end = c(25, 35, 45)),
                      strand = Rle(strand(rep("+", 3))))
names(ORFranges) = rep("tx1_1",3)
names(ORFranges2) = rep("tx1_2",3)
grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2)
gr <- unlist(grl, use.names = F)

test_that("groupGRangesBy works as intended", {

  grltest <- groupGRangesBy(gr)
  expect_is(grltest,"GRangesList")
  expect_equal(length(grltest), 2)
  expect_equal(length(unlist(grl[1], use.names = F)), 3)
})

test_that("tile1 works as intended", {

  tilex <- tile1(grl)
  expect_is(tilex,"GRangesList")
  expect_equal(as.integer(unlist(start(tilex[1]))),
    c(1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15, 20, 21, 22, 23, 24, 25))
  expect_equal(as.integer(unlist(end(tilex[2]))),
    c(20, 21, 22, 23, 24, 25, 30, 31, 32, 33, 34, 35, 40, 41, 42, 43, 44, 45))
})

test_that("widthPerGroup works as intended", {

  widths <- widthPerGroup(grl, F)
  expect_is(widths,"integer")
  expect_equal(widths, c(17,18))
})

test_that("firstExonPerGroup works as intended", {

  firstExons <- firstExonPerGroup(grl)
  expect_is(firstExons,"GRangesList")
  expect_equal(length(firstExons), 2)
  expect_equal(as.integer(unlist(start(firstExons))), c(1,20))
})






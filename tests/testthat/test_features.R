library(ORFik)
context("features")

RFP <- GRanges(seqnames = Rle(rep("1", 5)),
                            ranges = IRanges(start = c(1, 10, 20, 30, 40),
                                             end = c(5, 15, 25, 35, 45)),
                            strand = Rle(strand(rep("+", 5))))
RFP2 <- GRanges(seqnames = Rle(rep("1", 2)),
               ranges = IRanges(start = c(1000, 1010),
                                end = c(1005, 1015)),
               strand = Rle(strand(rep("+", 2))))

RFP3 <- GRanges(seqnames = Rle(rep("1", 3)),
                ranges = IRanges(start = c(1, 10, 20),
                                 end = c(5, 15, 25)),
                strand = Rle(strand(rep("+", 3))))

RFP4 <- GRanges(seqnames = Rle(rep("1", 3)),
                ranges = IRanges(start = c(1, 4, 12),
                                 end = c(1, 4, 12)),
                strand = Rle(strand(rep("+", 3))))

ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                     strand = Rle(strand(rep("+", 3))))

ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(20, 30, 40), end = c(25, 35, 45)),
                      strand = Rle(strand(rep("+", 3))))

ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(30, 40, 50), end = c(35, 45, 55)),
                      strand = Rle(strand(rep("+", 3))))
names(ORFranges) = rep("tx1_1" ,3)
names(ORFranges2) = rep("tx1_2", 3)
names(ORFranges3) = rep("tx1_3", 3)
grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2, tx1_3 = ORFranges3)

test_that("entropy works as intended", {

  entropy <- entropy(grl, RFP)
  expect_is(entropy,"numeric")
  expect_equal(round(entropy, 2), c(1.00, 1.00, 0.93))

  entropy <- entropy(grl, RFP2)
  expect_is(entropy,"numeric")
  expect_equal(entropy, c(0, 0, 0))

  entropy <- entropy(grl, RFP3)
  expect_is(entropy,"numeric")
  expect_equal(round(entropy, 2), c(1.00, 0.55, 0.00))
})

test_that("ORFScore works as intended", {

  scores <- ORFScores(grl, RFP)
  expect_is(scores,"list")
  expect_equal(as.numeric(scores$ORFscore), c(0, 0, 0))

  scores <- ORFScores(grl, RFP4)
  expect_is(scores,"list")
  expect_equal(as.numeric(round(scores$ORFscore, 2)), c(0.58, 0.00, 0.00))

})



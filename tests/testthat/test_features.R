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

RFP5 <- GRanges(seqnames = Rle(rep("1", 6)),
                ranges = IRanges(start = c(1, 4, 30, 60, 80, 90),
                                 end = c(30, 33, 63, 90, 110, 120)),
                strand = Rle(strand(rep("+", 6))))

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

tx_len <- c(100, 200, 300)
names(tx_len) <- c(ORFik:::OrfToTxNames(grl, unique = T)[1], "tx2", "tx3")

test_that("fractionLength works as intended", {

  scores <- ORFik:::fractionLength(grl, tx_len)
  expect_is(scores,"numeric")
  expect_equal(scores, c(0.17, 0.18, 0.18))


  tx_len[2] <- 300
  scores <- ORFik:::fractionLength(grl, tx_len)
  expect_is(scores,"numeric")
  expect_equal(scores, c(0.17, 0.18, 0.18))

  tx_len[1] <- 200
  scores <- ORFik:::fractionLength(grl, tx_len)
  expect_is(scores,"numeric")
  expect_equal(scores, c(0.085, 0.090, 0.090))
})

tx1 <- GRanges(seqnames = Rle(rep("1", 6)),
               ranges = IRanges(start = c(1, 10, 20, 30, 40, 50),
                                end = c(5, 15, 25, 35, 45, 55)),
               strand = Rle(strand(rep("+", 6))))
tx2 <- GRanges(seqnames = Rle(rep("1", 1)),
               ranges = IRanges(start = c(100),
                                end = c(300)),
               strand = Rle(strand(rep("+", 1))))
tx3 <- GRanges(seqnames = Rle(rep("1", 1)),
               ranges = IRanges(start = c(400),
                                end = c(700)),
               strand = Rle(strand(rep("+", 1))))
tx <- GRangesList(tx1 = tx1, tx2 = tx2, tx3 = tx3)

test_that("disengagementScore works as intended", {

  scores <- ORFik:::disengagementScore(grl, RFP, GtfOrTx = tx)
  expect_is(scores,"numeric")
  expect_equal(scores, c(1, 2, 3))

  scores <- ORFik:::disengagementScore(grl, RFP2, GtfOrTx = tx)
  expect_is(scores,"numeric")
  expect_equal(scores, c(1, 1, 1))
})

threeUTRs <- GRangesList(tx1 = ORFranges3)
names(threeUTRs) <- "tx1"
test_that("RibosomeReleaseScore works as intended", {
  grl <- grl[1]
  scores <- ORFik:::RibosomeReleaseScore(grl, RFP, threeUTRs)
  expect_is(scores,"numeric")
  expect_equal(round(scores,2), 1.41)

  scores <- ORFik:::RibosomeReleaseScore(grl, RFP2, threeUTRs)
  expect_is(scores,"numeric")
  expect_equal(round(scores,2), 1.06)
})

cds_1 <- GRanges(seqnames = Rle(rep("1", 1)),
                 ranges = IRanges(start = c(40),
                                  end = c(1020)),
                 strand = Rle(strand(rep("+", 1))))
cds <- GRangesList(cds_1)
names(cds) <- "tx1"
test_that("floss works as intended", {

  scores <- ORFik:::floss(grl, RFP5, cds, 26, 34)
  expect_is(scores,"numeric")
  expect_equal(round(scores,2), c(0.25, 0.08, 0.08))

  scores <- ORFik:::floss(grl, RFP2, cds, 26, 34)
  expect_is(scores,"numeric")
  expect_equal(scores, c(0, 0, 0))
})

RNA  <- GRanges(seqnames = Rle(rep("1", 3)),
                ranges = IRanges(start = c(1,100,400),
                                 end = c(55,300,700)),
                strand = Rle(strand(rep("+", 3))))
test_that("te works as intended", {
  tx_len <- width(RNA)
  names(tx_len) <- names(tx)
  scores <- ORFik:::te(grl, RNA, RFP, tx_len)
  expect_is(scores,"numeric")
  expect_equal(round(scores,2), c(5.82, 5.50, 3.67))

  scores <- ORFik:::te(grl, RNA, RFP2, tx_len)
  expect_is(scores,"numeric")
  expect_equal(scores, c(0, 0, 0))
})


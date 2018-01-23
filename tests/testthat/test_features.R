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

tx_len <- c(100, 200, 300)
names(tx_len) <- names(tx)

fiveUTRs1 <- GRanges(seqnames = Rle(rep("1", 4)),
                     ranges = IRanges(start = c(1, 10, 20, 30),
                                      end = c(5, 15, 25, 35)),
                     strand = Rle(strand(rep("+", 4))))
fiveUTRs2 <- GRanges(seqnames = Rle(rep("1", 2)),
                     ranges = IRanges(start = c(150, 180),
                                      end = c(170, 190)),
                     strand = Rle(strand(rep("+", 2))))

fiveUTRs <- GRangesList(tx1 = fiveUTRs1, tx2 = fiveUTRs2)

cds_1 <- GRanges(seqnames = Rle(rep("1", 1)),
                 ranges = IRanges(start = c(40),
                                  end = c(1020)),
                 strand = Rle(strand(rep("+", 1))))
cds <- GRangesList(tx1 = cds_1)


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



test_that("disengagementScore works as intended", {

  scores <- disengagementScore(grl, RFP, GtfOrTx = tx)
  expect_is(scores,"numeric")
  expect_equal(scores, c(1, 2, 3))

  scores <- disengagementScore(grl, RFP2, GtfOrTx = tx)
  expect_is(scores,"numeric")
  expect_equal(scores, c(1, 1, 1))
})

threeUTRs <- GRangesList(tx1 = ORFranges3)
names(threeUTRs) <- "tx1"
test_that("RibosomeReleaseScore works as intended", {
  grl <- grl[1]
  scores <- RibosomeReleaseScore(grl, RFP, threeUTRs)
  expect_is(scores,"numeric")
  expect_equal(round(scores,2), 1.41)

  scores <- RibosomeReleaseScore(grl, RFP2, threeUTRs)
  expect_is(scores,"numeric")
  expect_equal(round(scores,2), 1.06)
})

test_that("floss works as intended", {

  scores <- floss(grl, RFP5, cds, 26, 34)
  expect_is(scores,"numeric")
  expect_equal(round(scores,2), c(0.25, 0.08, 0.08))

  scores <- floss(grl, RFP2, cds, 26, 34)
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
  scores <- te(grl, RNA, RFP, tx_len)
  expect_is(scores,"numeric")
  expect_equal(round(scores,2), c(5.82, 5.50, 3.67))

  scores <- te(grl, RNA, RFP2, tx_len)
  expect_is(scores,"numeric")
  expect_equal(scores, c(0, 0, 0))
})

test_that("insideOutsideORF works as intended", {

  scores <- insideOutsideORF(grl, RFP, tx)
  expect_is(scores,"numeric")
  expect_equal(scores, c(0.8, 0.8, 0.6))

  scores <- insideOutsideORF(grl, RFP2, tx)
  expect_is(scores,"numeric")
  expect_equal(scores, c(1, 1, 1))
})

cds_1 <- GRanges(seqnames = Rle(rep("1", 1)),
                 ranges = IRanges(start = c(56),
                                  end = c(200)),
                 strand = Rle(strand(rep("+", 1))))
cds_2 <- GRanges(seqnames = Rle(rep("1", 1)),
                 ranges = IRanges(start = c(151),
                                  end = c(191)),
                 strand = Rle(strand(rep("+", 1))))
cds <- GRangesList(tx1 = cds_1, tx2 = cds_2)


fiveUTRs1 <- GRanges(seqnames = Rle(rep("1", 6)),
               ranges = IRanges(start = c(1, 10, 20, 30, 40, 50),
                                end = c(5, 15, 25, 35, 45, 55)),
               strand = Rle(strand(rep("+", 6))))
fiveUTRs2 <- GRanges(seqnames = Rle(rep("1", 2)),
                    ranges = IRanges(start = c(150, 180),
                                     end = c(170, 190)),
                    strand = Rle(strand(rep("+", 2))))

fiveUTRsOri <- GRangesList(tx1 = fiveUTRs1, tx2 = fiveUTRs2)


test_that("distOrfToCds works as intended", {

  scores <- distOrfToCds(grl, fiveUTRsOri, extension = 0)
  expect_is(scores,"numeric")
  expect_equal(scores, c(19, 7, 1))

  scores <- distOrfToCds(grl, fiveUTRsOri, cds, 0)
  expect_is(scores,"numeric")
  expect_equal(scores, c(19, 7, 1))

  scores <- distOrfToCds(grl, fiveUTRsOri, cds, 5)
  expect_is(scores,"numeric")
  expect_equal(scores, c(19, 7, 1))
})


scores <- distOrfToCds(grl, fiveUTRsOri, extension = 0)
test_that("inFrameWithCDS works as intended", {

  inFrame <- inFrameWithCDS(scores)
  expect_is(inFrame,"numeric")
  expect_equal(inFrame, c(1, 1, 1))
})

test_that("isOverlappingCds works as intended", {

  overlaps <- isOverlappingCds(scores)
  expect_is(overlaps,"logical")
  expect_equal(overlaps, c(F, F, F))
})

test_that("OrfRankOrder works as intended", {

  ranks <- OrfRankOrder(grl)
  expect_is(ranks,"integer")
  expect_equal(ranks, c(1, 2, 3))
})


test_that("allFeatures works as intended", {
  dt <- ORFik:::allFeatures(grl = grl,orfFeatures = T, RFP = RFP5, RNA = RNA,
                            tx = tx, fiveUTRs = fiveUTRsOri, cds = cds,
                            threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                            extension = 5)
  expect_is(dt,"data.table")
  expect_equal(ncol(dt), 10)
  expect_equal(nrow(dt), 3)

})


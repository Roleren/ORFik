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
RFP6 <- GRanges(seqnames = Rle(rep("1", 6)),
                ranges = IRanges(start = c(1, 4, 30, 60, 80, 90),
                                 end = c(30, 33, 63, 90, 110, 120)),
                strand = Rle(strand(rep("-", 6))))
RFP7 <- c(RFP5, RFP6)

RFP5GAlign <- as(RFP5, "GAlignments")



ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                     strand = Rle(strand(rep("+", 3))))

ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(20, 30, 40), end = c(25, 35, 45)),
                      strand = Rle(strand(rep("+", 3))))

ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(30, 40, 50), end = c(35, 45, 55)),
                      strand = Rle(strand(rep("+", 3))))
ORFranges4 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(30, 40, 50), end = c(35, 45, 55)),
                      strand = Rle(strand(rep("-", 3))))
ORFranges4 <- sort(ORFranges4, decreasing = TRUE)
names(ORFranges) <- rep("tx1_1" ,3)
names(ORFranges2) <- rep("tx1_2", 3)
names(ORFranges3) <- rep("tx1_3", 3)
names(ORFranges4) <- rep("tx4_1", 3)
grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2,
                   tx1_3 = ORFranges3, tx4_1 = ORFranges4)



tx1 <- GRanges(seqnames = Rle(rep("1", 6)),
               ranges = IRanges(start = c(1, 10, 20, 30, 40, 50),
                                end = c(5, 15, 25, 35, 45, 200)),
               strand = Rle(strand(rep("+", 6))))
tx2 <- GRanges(seqnames = Rle(rep("1", 1)),
               ranges = IRanges(start = c(250, 280),
                                end = c(270, 310)),
               strand = Rle(strand(rep("+", 1))))
tx3 <- GRanges(seqnames = Rle(rep("1", 1)),
               ranges = IRanges(start = c(400),
                                end = c(700)),
               strand = Rle(strand(rep("+", 1))))
tx4 <- GRanges(seqnames = Rle(rep("1", 6)),
               ranges = IRanges(start = c(1, 150, 160, 170, 180, 190),
                                end = c(145, 155, 165, 175, 185, 195)),
               strand = Rle(strand(rep("-", 6))))
tx4 <- sort(tx4, decreasing = TRUE)
tx <- GRangesList(tx1 = tx1, tx2 = tx2, tx3 = tx3, tx4 = tx4)

tx_len <- ORFik:::widthPerGroup(tx, TRUE)

RNA  <- GRanges(seqnames = Rle(rep("1", 3)),
                ranges = IRanges(start = c(1,250,400),
                                 end = c(200,310,700)),
                strand = Rle(strand(rep("+", 3))))
RNA2 <- c(RNA, GRanges(seqnames = Rle(rep("1", 1)),
                       ranges = IRanges(start = c(1),
                                        end = c(195)),
                       strand = Rle(strand(rep("-", 1)))))


RNAGAlign   <- as(RNA, "GAlignments")


fiveUTRs1 <- GRanges(seqnames = Rle(rep("1", 6)),
                     ranges = IRanges(start = c(1, 10, 20, 30, 40, 50),
                                      end = c(5, 15, 25, 35, 45, 55)),
                     strand = Rle(strand(rep("+", 6))))
fiveUTRs2 <- GRanges(seqnames = Rle(rep("1", 2)),
                     ranges = IRanges(start = c(250, 280),
                                      end = c(270, 290)),
                     strand = Rle(strand(rep("+", 2))))
fiveUTRs4 <- GRanges(seqnames = Rle(rep("1", 6)),
                     ranges = IRanges(start = c(140, 150, 160, 170, 180, 190),
                                      end = c(145, 155, 165, 175, 185, 195)),
                     strand = Rle(strand(rep("-", 6))))
fiveUTRs4 <- sort(fiveUTRs4, decreasing = TRUE)
fiveUTRs <- GRangesList(tx1 = fiveUTRs1, tx2 = fiveUTRs2, tx4 = fiveUTRs4)


cds1 <- GRanges(seqnames = Rle(rep("1", 1)),
                 ranges = IRanges(start = c(56),
                                  end = c(150)),
                 strand = Rle(strand(rep("+", 1))))
cds2 <- GRanges(seqnames = Rle(rep("1", 1)),
                 ranges = IRanges(start = c(291),
                                  end = c(302)),
                 strand = Rle(strand(rep("+", 1))))
cds4 <- GRanges(seqnames = Rle(rep("1", 1)),
                 ranges = IRanges(start = c(50),
                                  end = c(139)),
                 strand = Rle(strand(rep("-", 1))))
cds4 <- sort(cds4, decreasing = TRUE)
cds <- GRangesList(tx1 = cds1, tx2 = cds2, tx4 = cds4)
threeUTRs1 <- GRanges(seqnames = Rle(rep("1", 1)),
                     ranges = IRanges(start = c(151),
                                      end = c(200)),
                     strand = Rle(strand(rep("+", 1))))
threeUTRs2 <- GRanges(seqnames = Rle(rep("1", 1)),
                      ranges = IRanges(start = c(303),
                                       end = c(310)),
                      strand = Rle(strand(rep("+", 1))))
threeUTRs4 <- GRanges(seqnames = Rle(rep("1", 1)),
                      ranges = IRanges(start = c(1),
                                       end = c(49)),
                      strand = Rle(strand(rep("-", 1))))
threeUTRs4 <- sort(threeUTRs4, decreasing = TRUE)
threeUTRs <- GRangesList(tx1 = threeUTRs1, tx2 = threeUTRs2, tx4 = threeUTRs4)


test_that("entropy works as intended", {

  entropy <- entropy(grl, RFP)
  expect_is(entropy, "numeric")
  expect_equal(round(entropy, 2), c(1.00, 1.00, 0.93, 0.00))

  entropy <- entropy(grl, RFP2)
  expect_is(entropy, "numeric")
  expect_equal(entropy, c(0, 0, 0, 0))

  entropy <- entropy(grl, RFP3)
  expect_is(entropy, "numeric")
  expect_equal(round(entropy, 2), c(1.00, 0.55, 0.00, 0.00))

  entropy <- entropy(grl, RFP7)
  expect_is(entropy, "numeric")
  expect_equal(round(entropy, 2), c(0.93, 0.99, 0.97, 0.94))
})

test_that("orfScore works as intended", {

  scores <- orfScore(grl, RFP)
  expect_is(scores, "list")
  expect_equal(as.numeric(scores$ORFScores), c(0, 0, 0, 0))

  scores <- orfScore(grl, RFP4)
  expect_is(scores, "list")
  expect_equal(as.numeric(round(scores$ORFScores, 2)), c(0.58, 0.00, 0.00, 0.00))

  scores <- orfScore(grl, RFP7)
  expect_is(scores, "list")
  expect_equal(as.numeric(round(scores$ORFScores, 2)), c(0.00, 0.00, 0.36, 0.00))

})

test_that("fractionLength works as intended", {

  scores <- ORFik:::fractionLength(grl, tx_len)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(0.09, 0.10, 0.10, 0.10))

  changedtx_len <- tx_len
  changedtx_len[2] <- 300
  scores <- ORFik:::fractionLength(grl, changedtx_len)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(0.09, 0.10, 0.10, 0.10))

  changedtx_len[1] <- 200
  scores <- ORFik:::fractionLength(grl, changedtx_len)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(0.08, 0.09, 0.09, 0.10))
})



test_that("disengagementScore works as intended", {

  scores <- disengagementScore(grl, RFP, GtfOrTx = tx)
  expect_is(scores, "numeric")
  expect_equal(scores, c(1, 2, 3, 1))

  scores <- disengagementScore(grl, RFP2, GtfOrTx = tx)
  expect_is(scores, "numeric")
  expect_equal(scores, c(1, 1, 1, 1))

  scores <- disengagementScore(grl, RFP7, GtfOrTx = tx)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(0.43, 0.80, 0.80, 1.00))
})


test_that("ribosomeReleaseScore works as intended", {
  scores <- ribosomeReleaseScore(grl, RFP, threeUTRs)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(11.76, 11.11, 8.33, 2.72))

  scores <- ribosomeReleaseScore(grl, RFP2, threeUTRs)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(2.94, 2.78, 2.78, 2.72))

  scores <- ribosomeReleaseScore(grl, RFP7, threeUTRs)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(8.82, 11.11, 11.11, 2.72))
})

test_that("ribosomeStallingScore works as intended", {
  scores <- ribosomeStallingScore(grl, RFP)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(2.83, 3.00, 2.00, 6.00))

  scores <- ribosomeStallingScore(grl, RFP2)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(5.67, 6.00, 6.00, 6.00))

  scores <- ribosomeStallingScore(grl, RFP7)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(5.67, 3.00, 3.00, 6.00))
})



test_that("floss works as intended", {

  scores <- floss(grl, RFP5, cds, 26, 34)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(0.25, 0.08, 0.08, 0.00))

  scores <- floss(grl, RFP2, cds, 26, 34)
  expect_is(scores, "numeric")
  expect_equal(scores, c(0, 0, 0, 0))

  scores <- floss(grl, RFP7, cds, 26, 34)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(0.25, 0.08, 0.08, 0.08))
})



test_that("translationalEff works as intended", {

  scores <- translationalEff(grl, RNA, RFP, tx)
  expect_is(scores, "numeric")
  expect_equal(round(scores,2), c(19.06, 18.00, 12.00, NaN))

  scores <- translationalEff(grl, RNA, RFP2, tx)
  expect_is(scores, "numeric")
  expect_equal(scores, c(0, 0, 0, NaN))

  scores <- translationalEff(grl, RNA2, RFP7, tx)
  expect_is(scores, "numeric")
  expect_equal(round(scores,2), c(7.06, 10.00, 10.00, 9.72))
})

test_that("insideOutsideORF works as intended", {

  scores <- insideOutsideORF(grl, RFP, tx)
  expect_is(scores, "numeric")
  expect_equal(scores, c(0.8, 0.8, 0.6, 1.0))

  scores <- insideOutsideORF(grl, RFP2, tx)
  expect_is(scores, "numeric")
  expect_equal(scores, c(1, 1, 1, 1))

  scores <- insideOutsideORF(grl, RFP7, tx)
  expect_is(scores, "numeric")
  expect_equal(round(scores, 2), c(0.43, 0.57, 0.57, 0.57))
})

test_that("distToCds works as intended", {

  scores <- distToCds(grl, fiveUTRs, extension = 0)
  expect_is(scores, "numeric")
  expect_equal(scores, c(19, 7, 1, 38))

  scores <- distToCds(grl, fiveUTRs, cds, 0)
  expect_is(scores, "numeric")
  expect_equal(scores, c(19, 7, 1, 38))
  # TODO: Decide if this is correct behavior ?
  # My idea now is to keep it.
  scores <- distToCds(grl, fiveUTRs, cds, 5)
  expect_is(scores, "numeric")
  expect_equal(scores, c(19, 7, 1, 43))
})


scores <- distToCds(grl, fiveUTRs, extension = 0)
test_that("isInFrame works as intended", {

  inFrame <- isInFrame(scores)
  expect_is(inFrame, "numeric")
  expect_equal(inFrame, c(0, 0, 0, 1))
})

test_that("isOverlapping works as intended", {

  overlaps <- isOverlapping(scores)
  expect_is(overlaps, "logical")
  expect_equal(overlaps, c(F, F, F, F))
})

test_that("rankOrder works as intended", {

  ranks <- rankOrder(grl)
  expect_is(ranks, "integer")
  expect_equal(ranks, c(1, 2, 3, 1))
})


test_that("computeFeatures works as intended", {
  # without RNA in input
  dt <- computeFeaturesSpecial(grl = grl,orfFeatures = T, RFP = RFP5,
                            tx = tx, fiveUTRs = fiveUTRs, cds = cds,
                            threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                            extension = 0)
  expect_is(dt, "data.table")
  expect_equal(ncol(dt), 13)
  expect_equal(nrow(dt), 4)

  # normal inputs
  dt <- computeFeaturesSpecial(grl = grl,orfFeatures = T, RFP = RFP5, RNA = RNAGAlign,
                            tx = tx, fiveUTRs = fiveUTRs, cds = cds,
                            threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                            extension = 0)
  expect_is(dt, "data.table")
  expect_equal(ncol(dt), 15)
  expect_equal(nrow(dt), 4)

  dt <- computeFeaturesSpecial(grl = grl,orfFeatures = T, RFP = RFP5, RNA = RNA,
                            tx = tx, fiveUTRs = fiveUTRs, cds = cds,
                            threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                            extension = 5, cageFiveUTRs = fiveUTRs)
  expect_is(dt, "data.table")
  expect_equal(ncol(dt), 15)
  expect_equal(nrow(dt), 4)

  dt <- computeFeaturesSpecial(grl = grl,orfFeatures = T, RFP = RFP5GAlign, RNA = RNAGAlign,
                            tx = tx, fiveUTRs = fiveUTRs, cds = cds,
                            threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                            extension = 5, cageFiveUTRs = fiveUTRs)
  expect_is(dt, "data.table")
  expect_equal(ncol(dt), 15)
  expect_equal(nrow(dt), 4)

  # only nonvarying by Ribo-seq
  dt <- computeFeaturesSpecial(grl = grl,orfFeatures = T, RFP = RFP5GAlign, RNA = RNAGAlign,
                            tx = tx, fiveUTRs = fiveUTRs, cds = cds,
                            threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                            extension = 5, cageFiveUTRs = fiveUTRs, includeNonVarying = F)
  expect_is(dt, "data.table")
  expect_equal(ncol(dt), 10)
  expect_equal(nrow(dt), 4)

  # test from example table in orfik
  dt <- computeFeaturesSpecial(grl = grl,orfFeatures = T, RFP = RFP7, RNA = RNAGAlign,
                            tx = tx, fiveUTRs = fiveUTRs, cds = cds,
                            threeUTRs = threeUTRs, riboStart = 26, riboStop = 34,
                            extension = 0)
  expect_is(dt, "data.table")
  expect_equal(ncol(dt), 15)
  expect_equal(nrow(dt), 4)
  # load file
  load(system.file("extdata", "featureTable.rdata",
                         package = "ORFik")) # loaded as featureExamples
  expect_equal(dt, featureExamples) # should be equal to saved version
})


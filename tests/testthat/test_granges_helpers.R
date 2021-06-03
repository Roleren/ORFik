context("GRanges Helpers")
library(ORFik)

ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = c(1, 10, 20),
                                      end = c(5, 15, 25)),
                     strand = Rle(strand(rep("+", 3))))

ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(20, 30, 40),
                                       end = c(25, 35, 45)),
                      strand = Rle(strand(rep("+", 3))))
names(ORFranges) = rep("tx1_1",3)
names(ORFranges2) = rep("tx1_2",3)
grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2)
gr <- unlist(grl, use.names = FALSE)

cds1 <- ORFranges <- GRanges(seqnames = Rle(rep("1", 2)),
                             ranges = IRanges(start = c(100, 110),
                                              end = c(105, 115)),
                             strand = Rle(strand(rep("+", 2))))
cds2 <- ORFranges <- GRanges(seqnames = Rle(rep("1", 2)),
                             ranges = IRanges(start = c(200, 210),
                                              end = c(205, 215)),
                             strand = Rle(strand(rep("+", 2))))
cds <- GRangesList( tx1 = cds1, tx2 = cds2)
fiveUTRs1 <- GRanges(seqnames = Rle(rep("1", 6)),
                     ranges = IRanges(start = c(1, 10, 20, 30, 40, 50),
                                      end = c(5, 15, 25, 35, 45, 55)),
                     strand = Rle(strand(rep("+", 6))))
fiveUTRs2 <- GRanges(seqnames = Rle(rep("1", 2)),
                     ranges = IRanges(start = c(150, 180),
                                      end = c(170, 190)),
                     strand = Rle(strand(rep("+", 2))))
fiveUTRs <- GRangesList(tx1 = fiveUTRs1, tx2 = fiveUTRs2)

tx1 <-  GRanges(seqnames = Rle(rep("1", 5)),
                ranges = IRanges(start = c(1, 10, 20, 30, 40),
                                 end = c(5, 15, 25, 35, 45)),
                strand = Rle(strand(rep("+", 5))))

tx2 <- GRanges(seqnames = Rle(rep("1", 5)),
               ranges = IRanges(start = c(20, 30, 40, 50, 60),
                                end = c(25, 35, 45, 55, 65)),
               strand = Rle(strand(rep("+", 5))))
txl <- GRangesList(tx1 = tx1, tx2 = tx2)

txTst <- GRangesList(GRanges("1", IRanges(c(1,3,5), width = 1), "+"),
                     GRanges("1", IRanges(c(5,3,1), width = 1), "-"))

test_that("groupGRangesBy works as intended", {

  grltest <- groupGRangesBy(gr)
  expect_is(grltest,"GRangesList")
  expect_equal(length(grltest), 2)
  expect_equal(length(unlist(grl[1], use.names = FALSE)), 3)

  ggg <- GRanges(seqnames = "chrI",
                IRanges(start = c(10, 50, 100, 200), end = c(20,60,110,210)),
                strand = factor("+", levels = c("+", "-", "*")))
  names(ggg) <- c("a", "a", "b", "a")
  res <- groupGRangesBy(ggg)
  expect_equal(length(res), 2)
})

test_that("tile1 works as intended", {

  tilex <- tile1(grl)
  expect_is(tilex,"GRangesList")
  expect_equal(as.integer(unlist(start(tilex[1]))),
               c(1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15,
                 20, 21, 22, 23, 24, 25))
  expect_equal(as.integer(unlist(end(tilex[2]))),
               c(20, 21, 22, 23, 24, 25, 30, 31, 32, 33,
                 34, 35, 40, 41, 42, 43, 44, 45))
})

test_that("widthPerGroup works as intended", {

  widths <- widthPerGroup(grl, FALSE)
  expect_is(widths,"integer")
  expect_equal(widths, c(17,18))
})

test_that("firstExonPerGroup works as intended", {

  firstExons <- firstExonPerGroup(grl)
  expect_is(firstExons,"GRangesList")
  expect_equal(length(firstExons), 2)
  expect_equal(as.integer(unlist(start(firstExons))), c(1, 20))
})


test_that("lastExonPerGroup works as intended", {

  lastExons <- lastExonPerGroup(grl)
  expect_is(lastExons,"GRangesList")
  expect_equal(length(lastExons), 2)
  expect_equal(as.integer(unlist(start(lastExons))), c(20, 40))
})

test_that("strandBool works as intended", {

  strandLogical <- strandBool(grl)
  expect_is(strandLogical,"logical")
  expect_equal(length(strandLogical), 2)
  expect_equal(sum(strandLogical), 2)
})

test_that("assignFirstExonsStartSite works as intended", {
  newStarts <- as.integer(c(2, 21))
  reassigned <- assignFirstExonsStartSite(grl, newStarts)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(firstStartPerGroup(reassigned, FALSE), newStarts)
})

test_that("assignLastExonsStopSite works as intended", {
  newStops<- as.integer(c(26, 46))
  reassigned <- assignLastExonsStopSite(grl, newStops)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(lastExonEndPerGroup(reassigned, FALSE), newStops)
})


test_that("downstreamFromPerGroup works as intended", {
  downstreamFrom <- as.integer(c(3, 3))
  reassigned <- downstreamFromPerGroup(txTst, downstreamFrom)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(firstStartPerGroup(reassigned, FALSE), downstreamFrom)
})

test_that("downstreamOfPerGroup works as intended", {
  downstreamOf <- as.integer(c(11, 31))
  reassigned <- downstreamOfPerGroup(tx = grl, downstreamOf)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(firstStartPerGroup(reassigned, FALSE), downstreamOf + 1)
})

test_that("upstreamFromPerGroup works as intended", {
  upstreamFrom <- as.integer(c(10, 20, 41))
  reassigned <- upstreamFromPerGroup(tx = grl[c(1,1,2)], upstreamFrom)
  expect_is(reassigned,"GRangesList")
  expect_equal(stopSites(reassigned, is.sorted = TRUE), upstreamFrom)
})

test_that("upstreamOfPerGroup works as intended", {
  upstreamOf <- as.integer(c(12, 32))
  reassigned <- upstreamOfPerGroup(tx = grl, upstreamOf)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(lastExonEndPerGroup(reassigned, FALSE), upstreamOf)
  upstreamOf <- 20
  reassigned <- upstreamOfPerGroup(txl[1], 20)
  expect_equal(stopSites(reassigned), 15)
  reassigned2 <- upstreamOfPerGroup(txl[1], 20, FALSE)
  expect_equal(stopSites(reassigned2), 15)
})

test_that("asTX works as intended", {
  reassigned <- asTX(grl, txl)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(lastExonEndPerGroup(reassigned, FALSE), as.integer(c(17, 29)))
})

test_that("extendLeaders works as intended", {
  reassigned <- extendLeaders(fiveUTRs, 5, cds)
  expect_is(reassigned, "GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(firstStartPerGroup(reassigned, FALSE), as.integer(c(1, 145)))
  expect_equal(lastExonEndPerGroup(reassigned, FALSE), as.integer(c(115, 215)))

  circular_fives <- fiveUTRs
  isCircular(circular_fives) <- rep(TRUE, length(isCircular(circular_fives)))
  reassigned <- extendLeaders(circular_fives, 5, cds)
  expect_equal(firstStartPerGroup(reassigned, FALSE), as.integer(c(-4, 145)))
  expect_equal(lastExonEndPerGroup(reassigned, FALSE), as.integer(c(115, 215)))

  reassigned <- extendLeaders(fiveUTRs, fiveUTRs)
  expect_equal(startSites(fiveUTRs, is.sorted = TRUE),
               startSites(reassigned, is.sorted = TRUE))
})

test_that("extendTrailers works as intended", {
  reassigned <- extendTrailers(txl, 5)
  expect_is(reassigned, "GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(firstStartPerGroup(reassigned, FALSE), as.integer(c(1, 20)))
  expect_equal(lastExonEndPerGroup(reassigned, FALSE), as.integer(c(50, 70)))
})

test_that("matchNaming works as intended", {

  ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                       ranges = IRanges(start = c(1, 2, 3),
                                        end = c(1, 2, 3)),
                       strand = Rle(strand(rep("+", 3))))

  ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(4, 5, 7),
                                         end = c(4, 5, 7)),
                        strand = Rle(strand(rep("+", 3))))

  names(ORFranges) = rep("tx1_1",3)
  names(ORFranges2) = rep("tx1_2",3)
  grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2)
  gr <- unlist(grl, use.names = FALSE)
  test_result <- ORFik:::matchNaming(gr, grl)
  # should stay 0 meta columns
  expect_equal(ncol(elementMetadata(unlist(test_result))), 0)
  # create some example meta columns
  gr2 <- gr
  df <- DataFrame(matrix(NA, ncol = 3, nrow = length(gr2)))

  colnames(df) <- c("orf_id", "orf_name", "exon_id")
  class(df[,1]) <- "integer"
  class(df[,2]) <- "character"
  class(df[,3]) <- "integer"

  elementMetadata(gr2) <- df
  # should now loose all meta
  test_result <- ORFik:::matchNaming(gr2, grl)
  expect_equal(ncol(elementMetadata(unlist(test_result))), 0)


  grl2 <- groupGRangesBy(gr2)
  # should now get all meta data
  test_result <- ORFik:::matchNaming(gr, grl2)
  expect_equal(ncol(elementMetadata(unlist(test_result))), 3)
})

test_that("reduceKeepAttr works as intended", {

  ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                       ranges = IRanges(start = c(1, 2, 3),
                                        end = c(1, 2, 3)),
                       strand = Rle(strand(rep("+", 3))))

  ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                        ranges = IRanges(start = c(4, 5, 7),
                                         end = c(4, 5, 7)),
                        strand = Rle(strand(rep("+", 3))))

  names(ORFranges) = rep("tx1_1",3)
  names(ORFranges2) = rep("tx1_2",3)
  grl <- GRangesList(tx1_1 = ORFranges, tx1_2 = ORFranges2)
  reassigned <- reduceKeepAttr(grl, keep.names = TRUE)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  unlreassigned <- unlist(reassigned, use.names = FALSE)
  expect_equal(as.integer(start(unlreassigned)), c(1,4,7))
  expect_equal(names(unlreassigned), c("tx1_1", "tx1_2", "tx1_2"))
})


test_that("windowPerGroup works as intended", {

  gr <- GRanges(seqnames = "1", ranges = IRanges(start = c(40),end = c(40)),
                strand = "+")

  txgr <- GRanges(seqnames = "1", ranges = IRanges(start = c(20, 45, 100),
                                                   end = c(40, 70, 100)),
                  strand = "+")

  names(gr) = rep("tx1",1)
  names(txgr) = c(rep("tx1",2), "tx2")
  tx <- groupGRangesBy(txgr)

  test_result <- windowPerGroup(gr, tx, 20, 20)

  expect_equal(as.integer(unlist(start(test_result), use.names = FALSE)),
               c(20,45))
  expect_equal(as.integer(unlist(end(test_result), use.names = FALSE)),
               c(40, 64))

})

test_that("readWidths works as intended", {
  ga <- GAlignments(seqnames = "1", pos = as.integer(1), cigar = "1M1S",
                    strand = factor("+", levels = c("+", "-", "*")))

  expect_equal(readWidths(ga), 1) # With soft-clip
  expect_equal(readWidths(ga, after.softclips = FALSE), 2) # Without soft-clip
})

test_that("convertToOneBasedRanges works as intended", {
  # Soft clipping should not matter
  ga <- GAlignments(seqnames = "1", pos = as.integer(5), cigar = "22S6M",
                    strand = factor("+", levels = c("+", "-", "*")))
  ga2 <- GAlignments(seqnames = "1", pos = as.integer(5), cigar = "3S6M",
                    strand = factor("+", levels = c("+", "-", "*")))
  ga3 <- GAlignments(seqnames = "1", pos = as.integer(5), cigar = "4S6M",
                     strand = factor("+", levels = c("+", "-", "*")))

  ga <- c(rep(ga, 2), rep(ga2, 2), ga3)

  res <- convertToOneBasedRanges(ga, addScoreColumn = TRUE,
                                 addSizeColumn = TRUE)
  expect_equal(readWidths(res), 6)

  res <- convertToOneBasedRanges(ga, addScoreColumn = FALSE,
                                 addSizeColumn = TRUE)
  expect_equal(readWidths(res), rep(6, 5))

  res <- convertToOneBasedRanges(ga, addScoreColumn = TRUE,
                                 addSizeColumn = FALSE)
  expect_equal(score(res), 5)
  # Introns gaps matter
  ga <- GAlignments(seqnames = "1", pos = as.integer(5), cigar = "6M6N6M",
                    strand = factor("+", levels = c("+", "-", "*")))
  ga2 <- GAlignments(seqnames = "1", pos = as.integer(5), cigar = "6M7N6M",
                     strand = factor("+", levels = c("+", "-", "*")))
  ga3 <- GAlignments(seqnames = "1", pos = as.integer(5), cigar = "6M8N7M",
                     strand = factor("+", levels = c("+", "-", "*")))

  ga <- c(rep(ga, 2), rep(ga2, 2), ga3)
  res <- convertToOneBasedRanges(ga, addScoreColumn = TRUE,
                                 addSizeColumn = TRUE, method = "3prime")
  expect_equal(start(res), c(22, 23, 25))
})

test_that("pmapToTranscriptF works as intended", {
  res <- pmapToTranscriptF(grl, grl)
  expect_equal(startSites(res), c(1,1))
  expect_equal(stopSites(res), c(17,18))
  expect_equal(names(res), c("tx1_1", "tx1_2"))

  res <- pmapToTranscriptF(ranges(grl), grl)
  expect_equal(unlist(end(res), use.names = FALSE), c(17,18))
  res <- pmapToTranscriptF(stopSites(grl, asGR = TRUE), grl)
  expect_equal(end(res), c(17,18))
  res <- pmapToTranscriptF(ranges(stopSites(grl, asGR = TRUE)), grl)
  expect_equal(end(res), c(17,18))
})

test_that("pmapFromTranscriptF works as intended", {
  temp <- pmapToTranscriptF(grl, grl)
  x <- ranges(unlist(temp, use.names = TRUE))
  names(x) <- c(1,2)
  res <- pmapFromTranscriptF(x, grl)
  expect_equal(ranges(res), ranges(grl))
})

library(ORFik)
context("GRanges Helpers")



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
    c(1, 2, 3, 4, 5, 10, 11, 12, 13, 14, 15,
        20, 21, 22, 23, 24, 25))
  expect_equal(as.integer(unlist(end(tilex[2]))),
    c(20, 21, 22, 23, 24, 25, 30, 31, 32, 33,
        34, 35, 40, 41, 42, 43, 44, 45))
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
  expect_equal(firstStartPerGroup(reassigned, F), newStarts)
})

test_that("assignLastExonsStopSite works as intended", {
  newStops<- as.integer(c(26, 46))
  reassigned <- assignLastExonsStopSite(grl, newStops)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(lastExonEndPerGroup(reassigned, F), newStops)
})

test_that("downstreamOfPerGroup works as intended", {
  downstreamOf <- as.integer(c(11, 31))
  reassigned <- downstreamOfPerGroup(tx = grl, downstreamOf)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(firstStartPerGroup(reassigned, F), downstreamOf)
})

test_that("upstreamOfPerGroup works as intended", {
  upstreamOf <- as.integer(c(12, 32))
  reassigned <- upstreamOfPerGroup(tx = grl, upstreamOf)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(lastExonEndPerGroup(reassigned, F), upstreamOf)
})

tx1 <-  GRanges(seqnames = Rle(rep("1", 5)),
                                 ranges = IRanges(start = c(1, 10, 20, 30, 40),
                                                  end = c(5, 15, 25, 35, 45)),
                                 strand = Rle(strand(rep("+", 5))))

tx2 <- GRanges(seqnames = Rle(rep("1", 5)),
                      ranges = IRanges(start = c(20, 30, 40, 50, 60),
                                       end = c(25, 35, 45, 55, 65)),
                      strand = Rle(strand(rep("+", 5))))
txl <- GRangesList(tx1 = tx1, tx2 = tx2)

test_that("asTX works as intended", {
  reassigned <- asTX(grl, txl)
  expect_is(reassigned,"GRangesList")
  expect_equal(length(reassigned), 2)
  expect_equal(lastExonEndPerGroup(reassigned, F), as.integer(c(17,29)))
})




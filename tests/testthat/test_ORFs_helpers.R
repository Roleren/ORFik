library(ORFik)
context("ORF helpers")

transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
                            ranges = IRanges(start = c(1, 10, 20, 30, 40), end = c(5, 15, 25, 35, 45)),
                            strand = Rle(strand(rep("+", 5))))
ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
                     strand = Rle(strand(rep("+", 3))))
ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(10, 20, 30), end = c(15, 25, 35)),
                      strand = Rle(strand(rep("+", 3))))
ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = c(20, 30, 40), end = c(25, 35, 45)),
                      strand = Rle(strand(rep("+", 3))))

test_that("define_trailer works as intended for plus strand", {

  #at the start
  trailer <- define_trailer(ORFranges, transcriptRanges)
  expect_is(trailer, "GRanges")
  expect_equal(start(trailer), c(30, 40))
  expect_equal(end(trailer), c(35, 45))

  #middle
  trailer2 <- define_trailer(ORFranges2, transcriptRanges)
  expect_equal(start(trailer2), 40)
  expect_equal(end(trailer2), 45)

  #at the end
  trailer3 <- define_trailer(ORFranges3, transcriptRanges)
  expect_is(trailer3, "GRanges")
  expect_equal(length(trailer3), 0)

  #trailer size 3
  trailer4 <- define_trailer(ORFranges2, transcriptRanges, 3)
  expect_equal(start(trailer4), 40)
  expect_equal(end(trailer4), 42)
})


transcriptRanges <- GRanges(seqnames = Rle(rep("1", 5)),
                            ranges = IRanges(start = rev(c(1, 10, 20, 30, 40)), end = rev(c(5, 15, 25, 35, 45))),
                            strand = Rle(strand(rep("-", 5))))
ORFranges <- GRanges(seqnames = Rle(rep("1", 3)),
                     ranges = IRanges(start = rev(c(1, 10, 20)), end = rev(c(5, 15, 25))),
                     strand = Rle(strand(rep("-", 3))))
ORFranges2 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = rev(c(10, 20, 30)), end = rev(c(15, 25, 35))),
                      strand = Rle(strand(rep("-", 3))))
ORFranges3 <- GRanges(seqnames = Rle(rep("1", 3)),
                      ranges = IRanges(start = rev(c(20, 30, 40)), end = rev(c(25, 35, 45))),
                      strand = Rle(strand(rep("-", 3))))

test_that("define_trailer works as intended for minus strand", {

  #at the end
  trailer <- define_trailer(ORFranges, transcriptRanges)
  expect_is(trailer, "GRanges")
  expect_is(trailer, "GRanges")
  expect_equal(length(trailer), 0)

  #middle
  trailer2 <- define_trailer(ORFranges2, transcriptRanges)
  expect_equal(start(trailer2), 1)
  expect_equal(end(trailer2), 5)

  #at the start
  trailer3 <- define_trailer(ORFranges3, transcriptRanges)
  expect_equal(start(trailer3), c(1, 10))
  expect_equal(end(trailer3), c(5, 15))

  #trailer size 3
  trailer4 <- define_trailer(ORFranges2, transcriptRanges, 3)
  expect_equal(start(trailer4), 3)
  expect_equal(end(trailer4), 5)
})

transcriptRanges <- GRanges(seqnames = Rle(rep("1", 4)),
                            ranges = IRanges(start = rev(c(10, 20, 30, 40)), end = rev(c(15, 25, 35, 45))),
                            strand = Rle(strand(rep("-", 4))))

test_that("map_to_GRanges works as intended for minus strand", {
  mapback <- map_to_GRanges(IRanges(start = c(1, 10, 20), width = c(17, 5, 4)), transcriptRanges)
  #first ORF
  expect_equal(sum(width(mapback[mapback$names == "_1"])), 17)
  expect_equal(start(mapback[mapback$names == "_1"]), c(40, 30, 21))
  expect_equal(end(mapback[mapback$names == "_1"]), c(45, 35, 25))
  #second ORF
  expect_equal(sum(width(mapback[mapback$names == "_2"])), 5)
  expect_equal(start(mapback[mapback$names == "_2"]), c(30, 24))
  expect_equal(end(mapback[mapback$names == "_2"]), c(32, 25))
  #third ORF
  expect_equal(sum(width(mapback[mapback$names == "_3"])), 4)
  expect_equal(start(mapback[mapback$names == "_3"]), c(11))
  expect_equal(end(mapback[mapback$names == "_3"]), c(14))
})

strand(transcriptRanges) <- "+"
transcriptRanges <- sort(transcriptRanges)
test_that("map_to_GRanges works as intended for plus strand", {
  mapback <- map_to_GRanges(IRanges(start = c(1, 10, 20), width = c(17, 5, 4)), transcriptRanges)
  #first ORF
  expect_equal(sum(width(mapback[mapback$names == "_1"])), 17)
  expect_equal(start(mapback[mapback$names == "_1"]), c(10, 20, 30))
  expect_equal(end(mapback[mapback$names == "_1"]), c(15, 25, 34))
  #second ORF
  expect_equal(sum(width(mapback[mapback$names == "_2"])), 5)
  expect_equal(start(mapback[mapback$names == "_2"]), c(23, 30))
  expect_equal(end(mapback[mapback$names == "_2"]), c(25, 31))
  #third ORF
  expect_equal(sum(width(mapback[mapback$names == "_3"])), 4)
  expect_equal(start(mapback[mapback$names == "_3"]), c(41))
  expect_equal(end(mapback[mapback$names == "_3"]), c(44))
})

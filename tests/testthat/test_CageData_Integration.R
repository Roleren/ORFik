library(ORFik)
context("CageData Integration")

library(GenomicFeatures)
samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite", package = "GenomicFeatures")
txdb <- loadDb(samplefile)
fiveUTRs <- fiveUTRsByTranscript(txdb) # <- extract only 5' leaders
cds <- cdsBy(txdb)

cage <- GRanges(seqnames = as.character(seqnames(fiveUTRs)[1:2]), ranges =  IRanges(as.integer(start(fiveUTRs)[1:2]-500) ,
                as.integer(start(fiveUTRs)[1:2])), strand = as.character(strand(fiveUTRs)[1:2]), score = c(5,10))

cageEqualStart <- rep(cage[1],3)
cageEqualStart$score[3] <- 50
start(cageEqualStart)[3] <- start(cageEqualStart)[3] - 1


test_that("reassignTSSbyCage picks best one max peak of several", {
  # second granges have higher score, and both are within frame, then max one should be picked
  test_result <- reassignTSSbyCage(fiveUTRs[1], cds = cds, cageAsGR = cageEqualStart)

  expect_is(test_result, "GRangesList")
  expect_is(strand(test_result),"CompressedRleList")
  expect_is(seqnames(test_result),"CompressedRleList")
  expect_equal(length(test_result), 1)
  expect_equal(as.integer(unlist(start(test_result))), 32670735)
  expect_equal(as.integer(unlist(end(test_result))), 32671324)
  })


test_that("reassignTSSbyCage picks all needed peaks, no more, no less", {

  test_result <- reassignTSSbyCage(fiveUTRs[1:2], cds = cds, cageAsGR = cage)

  expect_is(test_result, "GRangesList")
  expect_is(strand(test_result),"CompressedRleList")
  expect_is(seqnames(test_result),"CompressedRleList")
  expect_equal(length(test_result), 2)
  expect_equal(as.integer(unlist(start(test_result))), c(32670736, 32670736))
  expect_equal(as.integer(unlist(end(test_result))), c(32671324, 32671564))

})



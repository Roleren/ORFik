context("CageData Integration")
library(ORFik)

library(GenomicFeatures)
samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                          package = "GenomicFeatures")
txdb <- loadDb(samplefile)
fiveUTRs <- fiveUTRsByTranscript(txdb) # <- extract only 5' leaders
cds <- cdsBy(txdb,"tx",use.names = TRUE)[1:length(fiveUTRs)]
names(cds) <- names(fiveUTRs)
rm(txdb)

cage <- GRanges(seqnames = as.character(seqnames(fiveUTRs)[1:2]),
                ranges =  IRanges(as.integer(start(fiveUTRs)[1 : 2] - 500) ,
                                  as.integer(start(fiveUTRs)[1 : 2])),
                strand = as.character(strand(fiveUTRs)[1 : 2]), score = c(5, 10))

cageEqualStart <- rep(cage[1], 3)
cageEqualStart$score[3] <- 50
start(cageEqualStart)[3] <- start(cageEqualStart)[3] - 1


test_that("reassignTSSbyCage picks best one max peak of several", {
  # second granges have higher score, and both are within frame, then max one should be picked
  test_result <- suppressWarnings(reassignTSSbyCage(fiveUTRs[1], cage = cageEqualStart, cds = cds))

  expect_is(test_result, "GRangesList")
  expect_is(strand(test_result),"CompressedRleList")
  expect_is(seqnames(test_result),"CompressedRleList")
  expect_equal(length(test_result), 1)
  expect_equal(as.integer(unlist(start(test_result))), 32670735)
  expect_equal(as.integer(unlist(end(test_result))), 32671324)
})


test_that("reassignTSSbyCage picks all needed peaks, no more, no less", {

  test_result <- reassignTSSbyCage(fiveUTRs[1:2], cage = cage, cds = cds )

  expect_is(test_result, "GRangesList")
  expect_is(strand(test_result), "CompressedRleList")
  expect_is(seqnames(test_result), "CompressedRleList")
  expect_equal(length(test_result), 2)
  expect_equal(as.integer(unlist(start(test_result))), c(32670736, 32670736))
  expect_equal(as.integer(unlist(end(test_result))), c(32671324, 32671564))

})

fiveAsGR <- unlist(fiveUTRs, use.names = TRUE)
cage <- GRanges(seqnames = c("1","1","2","2","3","3"),
                ranges =  IRanges(as.integer(start(fiveAsGR)[1 : 6] -500) ,
                                  as.integer(start(fiveAsGR)[1 : 6])),
                strand = as.character(strand(fiveAsGR)[1:6]),
                score = c(5, 10, 1, 2, 1, 2))

test_that("matchSeqlevels fixes seqlevel discrepancies correctly", {

  test_result <- matchSeqlevels(cage, fiveUTRs)

  expect_is(test_result, "GRanges")
  expect_equal(length(test_result), 6)
  expect_is(strand(test_result), "Rle")
  expect_is(seqnames(test_result), "Rle")
  expect_equal(seqlevels(test_result), c("chr1","chr2","chr3"))

})


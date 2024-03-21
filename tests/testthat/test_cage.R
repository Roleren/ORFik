context("CageData Integration")
library(ORFik)
library(GenomicFeatures)
samplefile <- system.file("extdata", "hg19_knownGene_sample.sqlite",
                          package = "GenomicFeatures")
txdb <- loadTxdb(samplefile)
loadRegions(txdb, c("leaders", "cds"))
cds <- cds[names(leaders)] # Subset to cds with leaders

cage <- GRanges(seqnames = as.character(seqnames(leaders)[1:2]),
                ranges =  IRanges(as.integer(start(leaders)[1 : 2] - 500) ,
                                  as.integer(start(leaders)[1 : 2])),
                strand = as.character(strand(leaders)[1 : 2]), score = c(5, 10))

cageEqualStart <- rep(cage[1], 3)
cageEqualStart$score[3] <- 50
start(cageEqualStart)[3] <- start(cageEqualStart)[3] - 1


test_that("reassignTSSbyCage picks best one max peak of several", {
  # second granges have higher score, and both are within frame,
  # then max one should be picked
  test_result <- suppressWarnings(reassignTSSbyCage(leaders[1],
                                                    cage = cageEqualStart))

  # expect_is(test_result, "GRangesList")
  expect_is(strand(test_result),"CompressedRleList")
  expect_is(seqnames(test_result),"CompressedRleList")
  expect_equal(length(test_result), 1)
  expect_equal(as.integer(unlist(start(test_result))), 32670735)
  expect_equal(as.integer(unlist(end(test_result))), 32671282)
})


test_that("reassignTSSbyCage picks all needed peaks, no more, no less", {

  test_result <- reassignTSSbyCage(leaders[1:2], cage = cage)

  expect_is(test_result, "GRangesList")
  expect_is(strand(test_result), "CompressedRleList")
  expect_is(seqnames(test_result), "CompressedRleList")
  expect_equal(length(test_result), 2)
  expect_equal(as.integer(unlist(start(test_result))), c(32670736, 32670736))
  expect_equal(as.integer(unlist(end(test_result))), c(32671282, 32671282))
})

fiveAsGR <- unlist(leaders, use.names = TRUE)
cage <- GRanges(seqnames = c("1","1","2","2","3","3"),
                ranges =  IRanges(as.integer(start(fiveAsGR)[1 : 6] -500) ,
                                  as.integer(start(fiveAsGR)[1 : 6])),
                strand = as.character(strand(fiveAsGR)[1:6]),
                score = c(5, 10, 1, 2, 1, 2))
extended <- GRanges("1", IRanges(c(32671236, 32671236-4, 32671236-6),
                                    width = 1), "+", score = c(1, 5, 5))
cage <- c(cage, extended)
if (seqlevels(leaders)[1] == "chr1")
  seqlevels(cage) <- paste0("chr", seqlevels(cage))

test_that("filterCage filters correctly", {

  test_result <- filterCage(cage, 1,leaders)
  expect_is(test_result, "GRanges")
  expect_equal(length(test_result), 5)
})

test_that("restrictTSSByUpstreamLeader filters correctly", {
  f <- leaders[1]
  start(f) <- IntegerList(start(f) - 800)
  end(f) <- IntegerList(end(f) - 800)
  ff <- c(leaders[1], f)
  shifted <- extendsTSSexons(ff, 1000)
  test_result <- restrictTSSByUpstreamLeader(ff, shifted)
  expect_is(test_result, "GRangesList")
  expect_equal(length(test_result), 2)

  expect_equal(as.integer(start(test_result[1])), as.integer(end(ff[2])))
})


test_that("addNewTSSOnLeaders assigns new TSS correctly", {
  test_result <- downstreamFromPerGroup(leaders[3], 32671324)

  expect_is(test_result, "GRangesList")
  test_result <- startSites(test_result, is.sorted = TRUE)
  expect_equal(test_result, 32671324)
})

test_that("reassignTSSbyCage removes correct leaders", {

  test_result <- reassignTSSbyCage(leaders[c(1:2,8)], cage = cage,
                                   removeUnused = TRUE)
  expect_equal(length(test_result), 2)
})

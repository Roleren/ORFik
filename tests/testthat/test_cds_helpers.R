library(ORFik)
context("cds helpers")

###########################NB########################
# Until AnnotationDbi version 1.41.4 is released
# a warning will appear, see:
# https://github.com/Bioconductor/AnnotationHub/issues/1
# remove this text when that happens.

gene1 <- GRanges(Rle(c("1"), c(4)),
                 IRanges(c(925942, 930155, 939040, 939275), width=c(72, 182, 90, 17)),
                 Rle(strand(c("+", "+", "+", "+"))))
gene2 <- GRanges(Rle(c("1"), c(4)),
                 IRanges(c(245929870, 245927931, 245863799, 245858562), width=c(32, 103, 88, 109)),
                 Rle(strand(c("-", "-", "-", "-"))))

cds <- GRangesList("gene_plus_strand" = gene1, "gene_minus_strand" = gene2)

test_that("startSites works as intended", {

  cds_starts <- startSites(cds, asGR = T)

  expect_is(cds_starts, "GRanges")
  expect_equal(start(cds_starts)[1], 925942)
  expect_equal(start(cds_starts)[2], 245929901)
})

test_that("stopSites works as intended", {

  cds_stops <- stopSites(cds, asGR = T)

  expect_is(cds_stops, "GRanges")
  expect_equal(start(cds_stops)[1], 939291)
  expect_equal(start(cds_stops)[2], 245858562)
})

test_that("windowResize works as intended", {

  cds_starts <- startSites(cds, asGR = T)
  resized <- windowResize(cds_starts, window_size = 30)

  expect_is(resized, "GRanges")
  expect_equal(start(resized)[1], 925912)
  expect_equal(end(resized)[1], 925972)
  expect_equal(start(resized)[2], 245929871)
  expect_equal(end(resized)[2], 245929931)
})

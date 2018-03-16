library(ORFik)
context("shift footprints")

gtf <- system.file("extdata", "example.gtf",
        package = "ORFik") ## location of the gtf file
suppressWarnings(txdb <-
 GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf"))
bam <- system.file("extdata", "example.bam",
        package = "ORFik") ## location of the bam file

footprints <- readFootprints(bam)


test_that("detectRibosomeShifts works as intended", {

  shifts <- detectRibosomeShifts(footprints, txdb, start=TRUE,
                                 stop=FALSE, offset_plots=FALSE, top_tx=10)
  expect_is(shifts, "data.frame")
  expect_equal(shifts$fragment_length, c(28,29,30))
  expect_equal(shifts$offsets_start, c(-11, -12, -13))
})


test_that("shiftFootprints works as intended", {
  shifts <- detectRibosomeShifts(footprints, txdb, start=TRUE,
                                 stop=FALSE, offset_plots=FALSE, top_tx=10)

  shiftedReads <- shiftFootprints(footprints, shifts$fragment_length,
                                  shifts$offsets_start)

  expect_is(shiftedReads, "GRanges")
  expect_equal(length(shiftedReads), length(footprints))
  expect_equal(start(shiftedReads[1]), 24066285)

})

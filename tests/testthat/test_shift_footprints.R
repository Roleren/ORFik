context("shift footprints")
library(ORFik)


gtf <- system.file("extdata", "annotations.gtf",
        package = "ORFik") ## location of the gtf file
suppressWarnings(txdb <-
 GenomicFeatures::makeTxDbFromGFF(gtf, format = "gtf"))
riboSeq_file <- system.file("extdata", "ribo-seq.bam",
                            package = "ORFik")

footprints <- GenomicAlignments::readGAlignments(
  riboSeq_file, param = ScanBamParam(flag = scanBamFlag(
    isDuplicate = FALSE, isSecondaryAlignment = FALSE)))

test_that("detectRibosomeShifts works as intended", {

  shifts <- detectRibosomeShifts(footprints, txdb)
  expect_is(shifts, "data.frame")
  expect_equal(shifts$fragment_length, c(28, 29, 30))
  expect_equal(shifts$offsets_start, c(-11, -12, -13))
})


test_that("shiftFootprints works as intended", {
  shifts <- detectRibosomeShifts(footprints, txdb)
  shiftedReads <- shiftFootprints(footprints, shifts$fragment_length,
                                  shifts$offsets_start)

  expect_is(shiftedReads, "GRanges")
  expect_equal(length(shiftedReads), length(footprints))
  expect_equal(start(shiftedReads[1]), 24066285)
})

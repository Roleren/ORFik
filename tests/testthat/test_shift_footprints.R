context("shift footprints")
library(ORFik)


gtf <- system.file("extdata", "annotations.gtf",
        package = "ORFik") ## location of the gtf file
suppressWarnings(txdb <-
 GenomicFeatures::makeTxDbFromGFF(gtf))
riboSeq_file <- system.file("extdata", "ribo-seq.bam",
                            package = "ORFik")

footprints <- GenomicAlignments::readGAlignments(
  riboSeq_file, param = ScanBamParam(flag = scanBamFlag(
    isDuplicate = FALSE, isSecondaryAlignment = FALSE)))

test_that("ribosome shifting works as intended", {

  shifts <- detectRibosomeShifts(footprints, txdb)
  expect_is(shifts, "data.frame")
  expect_equal(shifts$fraction, c(28, 29, 30))
  expect_equal(shifts$offsets_start, c(-11, -12, -13))

  # shiftFootprints
  shiftedReads <- shiftFootprints(footprints, shifts)

  # Shifted positive
  expect_is(shiftedReads, "GRanges")
  expect_equal(length(shiftedReads), length(footprints))
  expect_equal(start(shiftedReads[1]), 24066285)
  expect_equal(as.character(strand(shiftedReads[1])), "+")
  expect_equal(shiftedReads[1]$size, 28)

  # Shifted negative
  expect_equal(start(shiftedReads[16553]), 22711511)
  expect_equal(as.character(strand(shiftedReads[16553])), "-")
  expect_equal(shiftedReads[16553]$size, 30)
})

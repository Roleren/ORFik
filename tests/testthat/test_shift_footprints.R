context("shift footprints")
library(ORFik)
df <- ORFik.template.experiment.zf()
gtf <- df@txdb
suppressWarnings(txdb <-
 GenomicFeatures::makeTxDbFromGFF(gtf))
riboSeq_file <- filepath(df, "default")

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

test_that("ribosome shifting works as intended for gapped reads", {
  footprints <- GAlignments(seqnames = Rle(factor(c("1", "2"))),
                      pos = c(161524L, 336440L),
              cigar = c("20M469N9M", "11M207N18M"),
              strand = Rle(factor(c("+", "-"), levels = c("+", "-", "*"))))

  shift12 <- shiftFootprints(footprints,
                             data.frame(fraction = 29, offsets_start = -12))
  shift5 <- shiftFootprints(footprints,
                            data.frame(fraction = 29, offsets_start = -5))
  expect_equal(start(shift12), c(161536, 336663))
  expect_equal(start(shift5), c(161529, 336670))

  shift21 <- shiftFootprints(footprints,
                             data.frame(fraction = 29, offsets_start = -21))
  shift24 <- shiftFootprints(footprints,
                             data.frame(fraction = 29, offsets_start = -24))
  expect_equal(start(shift21), c(162014, 336447))
  expect_equal(start(shift24), c(162017, 336444))
})


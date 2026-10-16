context("shift footprints")
library(ORFik)
df <- ORFik.template.experiment.zf()
gtf <- df@txdb
suppressWarnings(txdb <- loadTxdb(gtf))
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

test_that("ribosome shifting works as intended for all types", {
  footprints <- GAlignments(seqnames = Rle(factor(rep(c("1", "2"), 5))),
                    pos = c(rep(3L, 4), rep(1L, 6)),
                    cigar = c("2S20M", "20M2S", "2S18M2S", "2S18M2S", "22M", "22M",
                              "10M10N10M", "10M10N10M", "10M2D10M","10M2D10M"),
                    strand = Rle(factor(rep(c("+", "-"), 5), levels = c("+", "-", "*"))))

  shift12 <- shiftFootprints(footprints,
    data.frame(fraction = c(18, 20, 22), offsets_start = -c(12,12,12)), FALSE)
  shift5 <- shiftFootprints(footprints,
    data.frame(fraction = c(18, 20, 22), offsets_start = -c(5,5,5)), FALSE)
  expect_equal(start(shift12), c(15,10,15,8,13,10,23,8,15,8))
  expect_equal(start(shift5), c(8,17,8,15,6,17,6,25,6,17))
})





test_that("shiftFootprintsByExperiment detects files missing from a partial shift.list", {
  # Regression: a shift.list covering only some files set txdb/tx/txNames to NULL for all files,
  # so the remaining files failed auto-detection ("asTX: not all references are present").
  d <- ORFik.template.experiment.zf()
  src_bam <- filepath(d, "default"); src_ofst <- filepath(d, "ofst")
  tmp <- file.path(tempdir(), "partial_shiftlist"); dir.create(file.path(tmp, "ofst"), recursive = TRUE, showWarnings = FALSE)
  for (n in c("ribo-seq_a", "ribo-seq_b")) {
    file.copy(src_bam, file.path(tmp, paste0(n, ".bam")), overwrite = TRUE)
    file.copy(src_ofst, file.path(tmp, "ofst", paste0(n, ".ofst")), overwrite = TRUE)
  }
  d2 <- d[c(1, 1), ]
  d2@listData$filepath <- file.path(tmp, c("ribo-seq_a.bam", "ribo-seq_b.bam"))
  d2@listData$rep <- c(1, 2)
  files <- filepath(d2, "ofst")
  known <- list(data.frame(fraction = c(28, 29, 30), offsets_start = c(-11, -12, -13)))
  names(known) <- files[1]
  out <- file.path(tmp, "pshifted/")
  expect_silent(suppressMessages(suppressWarnings(
    shiftFootprintsByExperiment(d2, out.dir = out, shift.list = known, output_format = "ofst",
                                accepted.lengths = 28:30, BPPARAM = BiocParallel::SerialParam()))))
  expect_true(all(file.exists(file.path(out, c("ribo-seq_a_pshifted.ofst", "ribo-seq_b_pshifted.ofst")))))
})

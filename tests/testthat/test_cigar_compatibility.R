test_that("readWidths preserves query and reference CIGAR semantics", {
  reads <- GAlignments(
    seqnames = Rle("chr1", 3), pos = c(1L, 100L, 300L),
    cigar = c("3S10M2I5M20N4D7M", "5M100N5M", "2H6M2S"),
    strand = Rle(strand(c("+", "-", "+")))
  )
  expect_identical(readWidths(reads), c(24L, 10L, 6L))
  expect_identical(readWidths(reads, after.softclips = FALSE), c(27L, 10L, 8L))
  expect_identical(readWidths(reads, along.reference = TRUE), c(26L, 10L, 6L))
})

test_that("reference CIGAR compatibility preserves intron and flag options", {
  widths <- ORFik:::cigarWidthAlongReferenceSpace_compat
  expect_identical(widths(c("10M20N5M", "3M2D7M")), c(35L, 12L))
  expect_identical(widths(c("10M20N5M", "3M2D7M"),
                          N.regions.removed = TRUE), c(15L, 12L))
  expect_identical(widths(c("10M", "10M"), flag = c(0L, 4L)), c(10L, NA_integer_))
  expect_identical(widths(character()), integer())
})

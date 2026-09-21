context("Export helpers")
library(ORFik)

test_that("export.bed12 writes expected grouped BED fields", {
  grl <- GRangesList(
    tx1 = GRanges("chr1", IRanges(c(5, 20), c(10, 25)), "+"),
    tx2 = GRanges("chr2", IRanges(c(100, 120), c(105, 130)), "-")
  )
  out <- file.path(tempdir(), "orfik_test_export.bed12")

  export.bed12(grl, out, rgb = c(0L, 255L))

  lines <- readLines(out)
  expect_length(lines, 2)

  fields1 <- strsplit(lines[1], "\t", fixed = TRUE)[[1]]
  fields2 <- strsplit(lines[2], "\t", fixed = TRUE)[[1]]

  expect_equal(fields1[1:6], c("chr1", "4", "25", "tx1", "12", "+"))
  expect_equal(fields1[7:12], c("4", "25", "0", "2", "6,6", "0,15"))

  expect_equal(fields2[1:6], c("chr2", "99", "130", "tx2", "17", "-"))
  expect_equal(fields2[7:12], c("99", "130", "255", "2", "6,11", "0,20"))
})

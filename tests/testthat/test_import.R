context("shift footprints")
library(ORFik)
df <- ORFik.template.experiment()
df_z <- ORFik.template.experiment.zf()

test_that("import ofst works", {
  reads <- fimport(filepath(df[1,], "default"))
  reads2 <- import.ofst(filepath(df[1,], "default"))
  expect_is(reads, "GRanges")
  expect_equal(length(reads), 102)
  expect_identical(reads, reads2)
})

test_that("import bam works", {
  reads <- fimport(filepath(df_z[1,], "default"))
  reads2 <- readBam(filepath(df_z[1,], "default"))
  expect_is(reads, "GAlignments")
  expect_equal(length(reads), 16649)
  expect_identical(reads, reads2)
})


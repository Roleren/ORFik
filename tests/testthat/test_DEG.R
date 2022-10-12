context("Differential_analysis")
library(ORFik)

df <- ORFik.template.experiment()
df.rna <- df[df$libtype == "RNA",]
df.rfp <- df[df$libtype == "RFP",]

test_that("DEG analysis works", {
  dt <- suppressWarnings(DEG.analysis(df.rna, output.dir = NULL))
  expect_is(dt, "data.table")
  expect_equal(nrow(dt), 6)
  expect_equal(unique(as.character(dt$Regulation)), "No change")
})

test_that("te.table works", {
  dt <- suppressWarnings(te.table(df.rfp, df.rna))
  expect_is(dt, "data.table")
  expect_equal(nrow(dt), 24)
})


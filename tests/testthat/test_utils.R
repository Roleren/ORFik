context("Utils")
library(ORFik)


test_that("NGS pair detection works as intended", {
  # As data.table
  paths <- c("/aligned_GRCz10/wigs/Bud_fwd.wig", "/aligned_GRCz10/wigs/Bud_MZDicer_fwd.wig",
             "/aligned_GRCz10/wigs/Bud_MZDicer_rev.wig", "/aligned_GRCz10/wigs/Bud_rev.wig")
  res <- findNGSPairs(paths)
  expect_equal(res$match[1], TRUE)
  # Compression
  paths <- paste0(paths, ".gz")
  res <- findNGSPairs(paths)
  expect_equal(res$match[1], TRUE)
  # No match
  paths[4] <- "/aligned_GRCz10/wigs/Bud1_rev.wig"
  res <- findNGSPairs(paths)
  expect_equal(res[4], paths[4])

  # Bed files
  paths <- c("1KCell_fwd.bed.gz", "2KCell_fwd.bed.gz",
             "2KCell_rev.bed.gz","1KCell_rev.bed.gz")
  res <- findNGSPairs(paths, format = "bed")
  expect_equal(res$match[1], TRUE)
})

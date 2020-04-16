context("Utils")
library(ORFik)


test_that("wig pair detection wokrs as intended", {
  # As data.table
  paths <- c("/aligned_GRCz10/wigs/Bud_fwd.wig", "/aligned_GRCz10/wigs/Bud_MZDicer_fwd.wig",
             "/aligned_GRCz10/wigs/Bud_MZDicer_rev.wig", "/aligned_GRCz10/wigs/Bud_rev.wig")
  res <- findWigPairs(paths)
  expect_equal(res$match[1], TRUE)
  # No match
  paths[4] <- "/aligned_GRCz10/wigs/Bud1_rev.wig"
  res <- findWigPairs(paths)
  expect_equal(res[4], paths[4])
})

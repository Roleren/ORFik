context("cell line naming patterns")
library(ORFik)

test_that("U-251MG is not labelled U2OS and maps to brain", {
  cl <- ORFik:::cellLineNames()
  x <- c("U-251MG ribo-seq rep1", "U251 cells", "U2OS_ribo", "U2-OS cells", "U2 OS", "rep U2", "HeLa U2 snRNA")
  got <- ORFik:::findFromPath(x, cl, "auto")
  expect_equal(got[1:2], c("U251", "U251"))
  expect_equal(got[3:6], rep("U2OS", 4))
  expect_equal(ORFik:::findFromPath("U251", ORFik:::cellLineNames(TRUE), "auto"), "brain")
  expect_equal(ORFik:::findFromPath("U2OS", ORFik:::cellLineNames(TRUE), "auto"), "bone")
})

test_that("cell line and tissue tables stay aligned", {
  expect_equal(nrow(ORFik:::cellLineNames()), nrow(ORFik:::cellLineNames(TRUE)))
})

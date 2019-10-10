context("Report Integration")
library(ORFik)

# Make test data
template <- create.experiment(dir = system.file("extdata", "", package = "ORFik"),
                  exper = "ORFik", txdb = system.file("extdata",
                                                      "annotations.gtf",
                                                      package = "ORFik"),
                  viewTemplate = FALSE)
template$X5[6] <- "heart"


test_that("Experiment class works as intended", {
  # test from example table in orfik
  df <- read.experiment(template)
  expect_equal(ncol(df), 6)
  # load file
  outputLibs(df)
  expect_equal(exists("ORFik_cage"), TRUE)
})


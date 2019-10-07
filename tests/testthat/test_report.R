context("Report Integration")
library(ORFik)

# Make test data
# create.experiment(dir = system.file("extdata", "", package = "ORFik"),
#                   exper = "ORFik", txdb = system.file("extdata",
#                                                       "annotations.gtf",
#                                                       package = "ORFik"),
#                   saveDir = "~/Desktop")

test_that("Experiment class works as intended", {
  # test from example table in orfik
  df <- read.experiment(system.file("extdata", "ORFik.csv",
                                            package = "ORFik"))
  expect_equal(ncol(df), 6)
  # load file
  outputLibs(df)
})


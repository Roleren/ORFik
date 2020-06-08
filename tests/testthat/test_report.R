context("Report Integration")
library(ORFik)

# Make test data
template <- create.experiment(dir = system.file("extdata", "", package = "ORFik"),
                  exper = "ORFik", txdb = system.file("extdata",
                                                      "annotations.gtf",
                                                      package = "ORFik"),
                  viewTemplate = FALSE)

df <- read.experiment(template)

test_that("Experiment class created as intended", {
  # test from example table in orfik
  expect_equal(ncol(df), 6)

})

test_that("Experiment class loaded as intended", {
  # load file
  outputLibs(df)
  expect_equal(exists("ORFik_CAGE_heart"), TRUE)
})

test_that("Experiment class correct naming", {
  # load file
  names <- bamVarName(df)
  expect_equal(names, c("ORFik_CAGE_heart", "ORFik_RFP_heart",
                        "ORFik_RFP","ORFik_RNA_heart"))
  names <- bamVarName(df, skip.experiment = TRUE)
  expect_equal(names, c("CAGE_heart", "RFP_heart",
                        "RFP","RNA_heart"))

  names <- bamVarName(df, skip.experiment = TRUE, skip.fraction = TRUE)
  expect_equal(names, c("CAGE_heart", "RFP_heart",
                        "RFP", "RNA_heart"))

  names <- bamVarName(df, skip.experiment = TRUE, skip.stage = TRUE)
  expect_equal(names, c("CAGE", "RFP",
                        "RFP","RNA"))
})

test_that("Experiment class correct renaming", {
  names <- c("run1", "_r2_", "rep3", "Rep5")
  dt <- repNames()
  res <- mainNames(names, dt)
  expect_equal(res, c("1", "2", "3", "5"))

  names <- c("rna-seq", "ribo-seq")
  dt <- libNames()
  res <- mainNames(names, dt)
  expect_equal(res, c("RNA", "RFP"))
})

# Count tables
test_that("count tables created as intended", {
  # Summairzed Experiment load
  expect_warning(SE <- makeSummarizedExperimentFromBam(df, region = "mrna"))
  expect_equal(assay(SE)[1,3], 8670)
  expect_warning(collapsed <- scoreSummarizedExperiment(SE, score = "count",
                                                        collapse = TRUE))
  expect_equal(assay(collapsed)[1,2], 8670)
  expect_warning(collapsed <- scoreSummarizedExperiment(SE, score = "count",
                                                        collapse = "all"))
  expect_equal(assay(collapsed)[3,1], 1971)
})

test_that("count tables loaded as intended", {
  # Summairzed Experiment load
  expect_warning(table <- countTable(df, "mrna"))
  expect_equal(table[1,3], data.table("ORFik_RFP" = 8670))
})

test_that("filepath work as intended", {
  # Summairzed Experiment load
  res <- filepath(df, "default")
  expect_equal(length(res), 4)

  # dff <- ORFik:::experiment(listData = list(libtype = c("RFP", "RFP"), rep = c(1,2),
  #                                           filepath = c("file1",
  #                                                        "file2"),
  #                                           reverse = c("", "paired-end")))

})

test_that("transcriptWindow plots correctly", {
  df <- df[3,]
  expect_warning(loadRegions(df))
  transcriptWindow(leaders, get("cds", mode = "S4"), trailers, df)
})


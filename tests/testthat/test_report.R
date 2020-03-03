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
  expect_equal(exists("ORFik_cage_heart"), TRUE)
})

test_that("Experiment class correct naming", {
  # load file
  names <- bamVarName(df)
  expect_equal(names, c("ORFik_cage_heart", "ORFik_ribo-seq_heart",
                        "ORFik_ribo-seq","ORFik_rna-seq_heart"))
  names <- bamVarName(df, skip.experiment = TRUE)
  expect_equal(names, c("cage_heart", "ribo-seq_heart",
                        "ribo-seq","rna-seq_heart"))

  names <- bamVarName(df, skip.experiment = TRUE, skip.fraction = TRUE)
  expect_equal(names, c("cage_heart", "ribo-seq_heart",
                        "ribo-seq","rna-seq_heart"))

  names <- bamVarName(df, skip.experiment = TRUE, skip.stage = TRUE)
  expect_equal(names, c("cage", "ribo-seq",
                        "ribo-seq","rna-seq"))
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
  expect_equal(table[1,3], data.table("ORFik_ribo-seq" = 8670))
})

test_that("transcriptWindow plots correctly", {
  df <- df[3,]
  expect_warning(loadRegions(df))
  transcriptWindow(leaders, get("cds", mode = "S4"), trailers, df)
})


context("Report Integration")
library(ORFik)

# Make test data
df <- ORFik.template.experiment()

# Count tables
test_that("count tables created as intended", {
  # Summairzed Experiment load
  SE <- makeSummarizedExperimentFromBam(df, region = "mrna")
  expect_gt(assay(SE)[1,3], 0)
  collapsed <- scoreSummarizedExperiment(SE, score = "count", collapse = TRUE)
  expect_gt(assay(collapsed)[1,2], 0)
  expect_warning(collapsed <- scoreSummarizedExperiment(SE, score = "count",
                                                        collapse = "all"))
  expect_gt(assay(collapsed)[3,1], 0)
})

test_that("count tables loaded as intended", {
  # Summairzed Experiment load
  #countTable_regions(df, out.dir = "~/Desktop/ORFik/extdata/Homo_sapiens_sample/QC_STATS/")
  table <- countTable(df, "mrna")
  expect_equal(table[1,3], data.table("CAGE_WT_r1" = 69))
  table <- countTable(df[2:3,], "mrna")
  expect_equal(colnames(table), c("CAGE_Mutant_r2", "CAGE_WT_r1"))
})

test_that("transcriptWindow plots correctly", {
  loadRegions(df)
  expect(exists("leaders", mode = "S4"), "Leaders does not exist in global env")
  expect(exists("cds", mode = "S4"), "CDS does not exist in global env")
  expect(exists("trailers", mode = "S4"), "trailers does not exist in global env")
  expect_warning(transcriptWindow(leaders, get("cds", mode = "S4"), trailers, df[3,],
                                  BPPARAM = BiocParallel::SerialParam()))

})

test_that("QCreport work as intended", {
  QCreport(df[9,], out.dir = tempdir(), BPPARAM = BiocParallel::SerialParam())
  expect(file.exists(file.path(tempdir(), "QC_STATS", "STATS.csv")),
         "QC_STATS/STATS.csv does not exist in tempdir")
})

context("Run ids in per-library output tables")
library(ORFik)

df <- ORFik.template.experiment()
df_run <- df
df_run@listData$Run <- sprintf("SRR%07d", seq_len(nrow(df)))

test_that("runIDsForOutput: NULL without (non-empty) Run, NA for empty entries", {
  expect_null(ORFik:::runIDsForOutput(df))
  expect_null(ORFik:::runIDsForOutput(NULL))
  expect_equal(ORFik:::runIDsForOutput(df_run), df_run$Run)
  expect_null(ORFik:::runIDsForOutput(df_run, n = nrow(df_run) - 1))
  df_empty <- df; df_empty@listData$Run <- rep("", nrow(df))
  expect_null(ORFik:::runIDsForOutput(df_empty))
  df_part <- df_run; df_part@listData$Run[2] <- ""
  expect_true(is.na(ORFik:::runIDsForOutput(df_part)[2]))
})

test_that("count table colData carries Run only when the experiment has run ids", {
  cd <- ORFik:::colDataFromExp(df)
  expect_false("Run" %in% colnames(cd))
  cd_run <- ORFik:::colDataFromExp(df_run)
  # Existing columns and order unchanged, Run appended last
  expect_equal(colnames(cd_run), c(colnames(cd), "Run"))
  expect_equal(cd_run$Run, df_run$Run)
  expect_equal(rownames(cd_run), rownames(cd))

  SE <- makeSummarizedExperimentFromBam(df_run, region = "mrna")
  expect_equal(SE$Run, df_run$Run)
  SE_old <- makeSummarizedExperimentFromBam(df, region = "mrna")
  expect_false("Run" %in% colnames(colData(SE_old)))
  expect_equal(assay(SE), assay(SE_old))
  # Collapsed columns list all their runs
  collapsed <- scoreSummarizedExperiment(SE, score = "count", collapse = TRUE)
  expected <- tapply(df_run$Run, as.character(SE$SAMPLE), paste, collapse = ";")
  expect_equal(collapsed$Run, unname(expected[colnames(collapsed)]))
  collapsed_all <- scoreSummarizedExperiment(SE, score = "count", collapse = "all")
  expect_equal(collapsed_all$Run, paste(df_run$Run, collapse = ";"))
})

test_that("colData of old (no Run) and new count tables can be combined", {
  a <- ORFik:::colDataFromExp(df_run[1:2, ])
  b <- ORFik:::colDataFromExp(df[3:4, ])
  ab <- ORFik:::rbind_colData_fill(a, b)
  expect_equal(nrow(ab), 4)
  expect_equal(ab$Run, c(df_run$Run[1:2], NA, NA))
  ba <- ORFik:::rbind_colData_fill(b, a)
  expect_equal(ba$Run, c(NA, NA, df_run$Run[1:2]))
})

test_that("readLengthTable adds Run as last column", {
  libs <- 9:10
  dt <- ORFik:::readLengthTable(df[libs, ], BPPARAM = BiocParallel::SerialParam())
  dt_run <- ORFik:::readLengthTable(df_run[libs, ], BPPARAM = BiocParallel::SerialParam())
  expect_false("Run" %in% colnames(dt))
  expect_equal(colnames(dt_run), c(colnames(dt), "Run"))
  expect_equal(dt_run[, colnames(dt), with = FALSE], dt)
  expect_equal(unique(dt_run[, .(sample_id, Run)])$Run, df_run$Run[libs])
})

test_that("orfFrameDistributions adds Run as last column", {
  dt <- orfFrameDistributions(df[9, ], BPPARAM = BiocParallel::SerialParam())
  dt_run <- orfFrameDistributions(df_run[9, ], BPPARAM = BiocParallel::SerialParam())
  expect_false("Run" %in% colnames(dt))
  expect_equal(colnames(dt_run), c(colnames(dt), "Run"))
  expect_equal(dt_run[, colnames(dt), with = FALSE], dt)
  expect_true(all(dt_run$Run == df_run$Run[9]))
})

test_that("alignmentFeatureStatistics adds Run as last column", {
  dt <- ORFik:::alignmentFeatureStatistics(df[9, ], BPPARAM = BiocParallel::SerialParam())
  dt_run <- ORFik:::alignmentFeatureStatistics(df_run[9, ], BPPARAM = BiocParallel::SerialParam())
  expect_false("Run" %in% colnames(dt))
  expect_equal(colnames(dt_run), c(colnames(dt), "Run"))
  expect_equal(dt_run$Run, df_run$Run[9])
})

test_that("QCstats.plot ignores the Run column of STATS.csv", {
  out <- file.path(tempdir(), "run_ids_qc")
  dir.create(out, showWarnings = FALSE)
  QCreport(df_run[9, ], out.dir = out, plot.ext = ".png", create.ofst = FALSE,
           complex.correlation.plots = FALSE, BPPARAM = BiocParallel::SerialParam())
  stats <- fread(file.path(out, "QC_STATS", "STATS.csv"))
  expect_equal(stats$Run, df_run$Run[9])
  expect_equal(tail(colnames(stats), 1), "Run")
  rl <- fread(file.path(out, "QC_STATS", "readLengths.csv"))
  expect_true(all(rl$Run == df_run$Run[9]))
  expect_silent(suppressMessages(QCstats.plot(file.path(out, "QC_STATS"), as_gg_list = TRUE)))
})

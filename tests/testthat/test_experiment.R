context("Experiment")
library(ORFik)
library(data.table)

df <- ORFik.template.experiment()
temp <- ORFik.template.experiment(as.temp = TRUE)
df_z <- ORFik.template.experiment.zf()

dir <- system.file("extdata/Homo_sapiens_sample", "", package = "ORFik")
exper <- "ORFik"
txdb <- system.file("extdata/references/homo_sapiens",
                    "Homo_sapiens_dummy.gtf.db", package = "ORFik")
fa <- system.file("extdata/references/homo_sapiens",
                  "Homo_sapiens_dummy.fasta", package = "ORFik")
org <- "Homo sapiens"
temp_dir <- tempdir()

test_that("Experiment template created as intended", {
  template <- create.experiment(dir = dir, exper, txdb = txdb,
                                saveDir = NULL,
                                fa = fa, organism = org)
  expect_equal(nrow(template), 20)
})

test_that("Experiment class loaded as intended", {
  # test from example table in orfik
  expect_equal(nrow(df), 16)
  expect_equal(ncol(df), 7)
})

test_that("Experiment saved as intended", {
  create.experiment(dir = dir, exper, txdb = txdb,
                    saveDir = temp_dir,
                    fa = fa, organism = org)
  df2 <- read.experiment(exper, temp_dir)
  expect_equal(nrow(df2), nrow(df))
  expect_equal(ncol(df2), ncol(df))
  expect_identical(temp$X6[-seq(4)], filepath(df2, "default"))
})

test_that("Experiment slot access works as intended", {
  # test from example table in orfik
  expect_equal(df$libtype, c(rep("CAGE", 4),rep("PAS", 4),rep("RFP", 4),rep("RNA", 4)))
  expect_equal(name(df), "ORFik")
  expect_equal(df[, "libtype"], c(rep("CAGE", 4),rep("PAS", 4),rep("RFP", 4),rep("RNA", 4)))
  expect_equal(df[df$libtype == "CAGE",]$libtype, rep("CAGE", 4))
})


test_that("output organism correctly", {
  expect_equal(organism(df), "Homo sapiens")
})

test_that("symbols work correctly", {

  cds_names <- names(loadRegion(df, "cds"))
  dt <- data.table(id = cds_names[-1], LFC = seq(5), p.value = 0.05)
  symbols_dt <- data.table(ensembl_transcript_name = cds_names,
   ensembl_gene_id = txNamesToGeneNames(cds_names, df),
   external_gene_name = c("ATF4", "AAT1", "ML4", "AST2", "RPL4", "RPL12"))
  suppressMessages(expect_equal(symbols(df), symbols_dt))
})

test_that("symbols work correctly empty", {
  suppressMessages(expect_equal(symbols(df_z), data.table()))
})

test_that("filepath work as intended", {
  # Summairzed Experiment load
  res <- filepath(df, "default")
  expect_equal(length(res), 16)
})

test_that("libFolder work as intended", {
  # Summairzed Experiment load
  res <- libFolder(df)
  expect_equal(basename(res), "Homo_sapiens_sample")
  res <- libFolder(df, mode = "unique")
  expect_equal(length(res), 1)
  res <- libFolder(df, mode = "all")
  expect_equal(length(res), 16)
})

test_that("validateExperiments() work as intended", {
  # Summairzed Experiment load
  validateExperiments(df)
  # Make it fail
  df2 <- df
  df2$libtype[10] <- "RNA"
  suppressMessages(expect_error(validateExperiments(df2)))
})

test_that("Show experiment correctly", {
  expect_identical(capture_output(df, print = TRUE),
                   "experiment: ORFik with 4 library types and 16 runs \nTjeldnes et al. \n    libtype rep condition\n 1:    CAGE   1    Mutant\n 2:    CAGE   2    Mutant\n 3:    CAGE   1        WT\n 4:    CAGE   2        WT\n 5:     PAS   1    Mutant\n 6:     PAS   2    Mutant\n 7:     PAS   1        WT\n 8:     PAS   2        WT\n 9:     RFP   1    Mutant\n10:     RFP   2    Mutant\n11:     RFP   1        WT\n12:     RFP   2        WT\n13:     RNA   1    Mutant\n14:     RNA   2    Mutant\n15:     RNA   1        WT\n16:     RNA   2        WT")
})

test_that("Experiment class loaded/removed as intended to custom environment", {
  # load file
  envExp(df) <- new.env()
  outputLibs(df[9,])
  expect_equal(exists("RFP", envir = envExp(df)), TRUE)
  expect_equal(exists("RFP", envir = .GlobalEnv), FALSE)
  # Remove
  remove.experiments(df[9,])
  expect_equal(exists("RFP", envir = envExp(df)), FALSE)
  expect_equal(exists("RFP", envir = .GlobalEnv), FALSE)
  envExp(df) <- .GlobalEnv
})

test_that("Experiment class loaded as intended", {
  # load file
  suppressMessages(outputLibs(df, BPPARAM = BiocParallel::SerialParam()))
  expect_equal(exists("CAGE_WT_r1"), TRUE)
})


test_that("Experiment class correct naming", {
  # load file
  names <- bamVarName(df)
  expect_equal(names, c("CAGE_Mutant_r1", "CAGE_Mutant_r2",
                        "CAGE_WT_r1", "CAGE_WT_r2",
                        "PAS_Mutant_r1", "PAS_Mutant_r2",
                        "PAS_WT_r1", "PAS_WT_r2",
                        "RFP_Mutant_r1", "RFP_Mutant_r2",
                        "RFP_WT_r1", "RFP_WT_r2",
                        "RNA_Mutant_r1", "RNA_Mutant_r2",
                        "RNA_WT_r1", "RNA_WT_r2"))
  names <- bamVarName(df, skip.experiment = FALSE)
  expect_equal(names, paste0("ORFik_", c("CAGE_Mutant_r1", "CAGE_Mutant_r2",
                                         "CAGE_WT_r1", "CAGE_WT_r2",
                                         "PAS_Mutant_r1", "PAS_Mutant_r2",
                                         "PAS_WT_r1", "PAS_WT_r2",
                                         "RFP_Mutant_r1", "RFP_Mutant_r2",
                                         "RFP_WT_r1", "RFP_WT_r2",
                                         "RNA_Mutant_r1", "RNA_Mutant_r2",
                                         "RNA_WT_r1", "RNA_WT_r2")))

  names <- bamVarName(df, skip.condition = TRUE)
  expect_equal(names, c("CAGE_r1", "CAGE_r2", "CAGE_r1", "CAGE_r2",
                        "PAS_r1",  "PAS_r2",  "PAS_r1",  "PAS_r2",
                        "RFP_r1",  "RFP_r2",  "RFP_r1",  "RFP_r2",
                        "RNA_r1", "RNA_r2", "RNA_r1", "RNA_r2"))

  names <- bamVarName(df, skip.rep = TRUE)
  expect_equal(names, c("CAGE_Mutant", "CAGE_Mutant", "CAGE_WT", "CAGE_WT",
                        "PAS_Mutant", "PAS_Mutant", "PAS_WT", "PAS_WT",
                        "RFP_Mutant", "RFP_Mutant", "RFP_WT", "RFP_WT",
                        "RNA_Mutant",  "RNA_Mutant", "RNA_WT", "RNA_WT"))
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

test_that("filepath find correct paths", {
  reads <- filepath(df_z[1,], "default")
  reads_as_ofst <- filepath(df_z[1,], "ofst")
  expect_is(reads, "character")
  expect_equal(length(reads), 1)
  expect_false(identical(reads, reads_as_ofst))
})




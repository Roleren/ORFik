context("Optimized txdb region files")
library(ORFik)

make_optimized_ref <- function() {
  dir <- tempfile("optimized_txdb_")
  dir.create(dir)
  src <- system.file("extdata/references/danio_rerio", package = "ORFik")
  gtf <- file.path(dir, "annotations.gtf")
  file.copy(file.path(src, "annotations.gtf"), gtf)
  suppressMessages(suppressWarnings(
    makeTxdbFromGenome(gtf, NULL, "Danio rerio", optimize = TRUE,
                       return = FALSE, gene_symbols = FALSE)))
  list(dir = dir, gtf = gtf, db = paste0(gtf, ".db"),
       opt_dir = file.path(dir, "ORFik_optimized"))
}

test_that("makeTxdbFromGenome(optimize = TRUE) writes files the loader finds", {
  ref <- make_optimized_ref()
  on.exit(unlink(ref$dir, recursive = TRUE))
  files <- list.files(ref$opt_dir)
  expect_true(length(files) > 0)
  # Writer must not name files after the .db path
  expect_false(any(grepl("\\.gtf\\.db_", files)))
  prefix <- ORFik:::optimized_txdb_path(ref$db, stop.error = FALSE)
  expect_true(file.exists(paste0(prefix, "_cds.qs")))
  expect_true(file.exists(paste0(prefix, "_txLengths.fst")))
  # Same prefix whether the txdb is given as path or as loaded object
  expect_equal(ORFik:::optimized_txdb_path(loadTxdb(ref$db), stop.error = FALSE),
               prefix)
  # loadRegion uses the optimized file: make it distinguishable
  cds <- loadRegion(ref$db, "cds")
  ORFik:::save_RDSQS(cds[1], paste0(prefix, "_cds.qs"))
  expect_length(loadRegion(ref$db, "cds"), 1)
  expect_length(loadRegion(ref$db, "cds", skip.optimized = TRUE), length(cds))
})

test_that("loader accepts legacy '<gtf>.db_<time>' and old '<gtf sans ext>_<time>' names", {
  ref <- make_optimized_ref()
  on.exit(unlink(ref$dir, recursive = TRUE))
  prefix <- ORFik:::optimized_txdb_path(ref$db, stop.error = FALSE)
  time <- sub(".*\\.gtf_", "", basename(prefix))
  files <- list.files(ref$opt_dir, full.names = TRUE)
  canonical <- paste0("annotations.gtf_", time)

  for (legacy in c(paste0("annotations.gtf.db_", time),
                   paste0("annotations_", time))) {
    new_names <- file.path(ref$opt_dir,
                           sub(canonical, legacy, basename(files), fixed = TRUE))
    file.rename(files, new_names)
    expect_equal(basename(ORFik:::optimized_txdb_path(ref$db, stop.error = FALSE)),
                 legacy)
    ORFik:::save_RDSQS(GRangesList(), file.path(ref$opt_dir, paste0(legacy, "_cds.qs")))
    expect_length(loadRegion(ref$db, "cds"), 0)
    files <- new_names
    canonical <- legacy
  }
})

# Streaming export (fst::concat_fst): block-aligned per-chromosome pieces joined without
# recompression must give exactly the assembled table.
skip_if_not(.ofst_store_concat_available(), "fst without concat_fst()")

export_fixture <- function(seed, n = 30000, chrs = c("1", "2", "10", "X", "MT", "KI1", "KI2"),
                           weights = c(8, 6, 5, 4, 1, 0.3, 0.2), double_scores = FALSE) {
  set.seed(seed)
  d <- data.table::data.table(seqnames = factor(sample(chrs, n, TRUE, prob = weights[seq_along(chrs)]), levels = chrs),
                              start = sample.int(60000L, n, TRUE),
                              strand = factor(sample(c("+", "-"), n, TRUE), levels = c("+", "-", "*")),
                              size = sample(26:31, n, TRUE),
                              score = if (double_scores) round(runif(n, 0.1, 9), 2) else sample.int(9L, n, TRUE))
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand, size)]
  data.table::setkey(d, NULL)
  d[sample(.N)]
}

export_write <- function(tables, names = paste0("lib", seq_along(tables))) {
  dir <- tempfile("export_in"); dir.create(dir)
  paths <- file.path(dir, paste0(names, ".ofst"))
  for (i in seq_along(tables)) fst::write_fst(tables[[i]], paths[i])
  paths
}

read_out <- function(f) fst::read_fst(f, as.data.table = TRUE)[]
merge_ref <- function(paths) suppressWarnings(suppressMessages(ofst_merge(paths, keep_all_scores = FALSE)))[]

test_that("streaming export equals the assembled table and a full merge", {
  paths <- export_write(lapply(1:4, export_fixture))
  store <- file.path(tempfile("ofst_store"), "store")
  ofst_store_build(store, paths, partition_min_rows = 5000, workers = 2, verbose = FALSE)
  m <- ofst_store_manifest(store)
  expect_true("pooled" %in% m$parts$part)
  expect_gt(max(m$parts$n_keys), 3 * 4096)
  ref <- merge_ref(paths)
  out <- tempfile(fileext = ".ofst")
  expect_identical(ofst_store_export_all(store, out, workers = 2, verbose = FALSE), out)
  expect_identical(read_out(out), ref)
  expect_identical(read_out(out), ofst_store_export_all(store, verbose = FALSE)[])
  expect_null(readRDS(paste0(out, ".removal_summary.rds")))
  expect_false(any(grepl("^\\.ofst_store_export_", list.files(dirname(out), all.files = TRUE))))
  # one worker (no forks) and the fallback (assemble + write_fst) give the same file content
  out1 <- tempfile(fileext = ".ofst")
  ofst_store_export_all(store, out1, workers = 1, verbose = FALSE)
  expect_identical(read_out(out1), ref)
  old <- options(ORFik.ofst_store_concat = FALSE)
  out2 <- tempfile(fileext = ".ofst")
  ofst_store_export_all(store, out2, workers = 2, verbose = FALSE)
  options(old)
  expect_identical(read_out(out2), ref)
  # subset of inputs
  out3 <- tempfile(fileext = ".ofst")
  ofst_store_export_all(store, out3, libs = c("lib2", "lib4"), workers = 2, verbose = FALSE)
  expect_identical(read_out(out3), merge_ref(paths[c(2, 4)]))
  # row ranges read across piece boundaries
  n <- nrow(ref)
  expect_identical(fst::read_fst(out, from = 4000, to = n - 3, as.data.table = TRUE)[], ref[4000:(n - 3)])
})

test_that("streaming export handles dropped chromosomes, double scores and tiny tables", {
  paths <- export_write(c(lapply(1:2, export_fixture, double_scores = TRUE),
                          list(export_fixture(3, n = 2000, chrs = c("Y", "KI9"), weights = c(1, 1)))))
  store <- file.path(tempfile("ofst_store"), "store")
  ofst_store_build(store, paths, partition_min_rows = 5000, workers = 2, verbose = FALSE)
  out <- tempfile(fileext = ".ofst")
  ofst_store_export_all(store, out, workers = 2, verbose = FALSE)
  ref <- merge_ref(paths)
  expect_true(is.double(ref$score))
  expect_identical(read_out(out), ref)
  ofst_store_remove(store, "lib3", workers = 2, verbose = FALSE)   # Y and KI9 disappear
  ofst_store_export_all(store, out, workers = 2, verbose = FALSE)
  got <- read_out(out)
  expect_identical(got, merge_ref(paths[1:2]))
  expect_false(any(c("Y", "KI9") %in% levels(got$seqnames)))
  # scores above int32 stay double
  big <- export_fixture(5, n = 9000)
  big[, score := as.numeric(score)][1, score := 3e9]
  paths2 <- export_write(list(big, export_fixture(6, n = 100)), c("big", "small"))
  store2 <- file.path(tempfile("ofst_store"), "store")
  ofst_store_build(store2, paths2, workers = 2, verbose = FALSE)
  ofst_store_export_all(store2, out, workers = 2, verbose = FALSE)
  expect_identical(read_out(out), merge_ref(paths2))
  expect_true(is.double(read_out(out)$score))
  # a table smaller than one block, and an empty export
  ofst_store_export_all(store2, out, libs = "small", workers = 2, verbose = FALSE)
  expect_identical(read_out(out), merge_ref(paths2[2]))
  ofst_store_remove(store2, c("big", "small"), verbose = FALSE)
  ofst_store_export_all(store2, out, workers = 2, verbose = FALSE)
  expect_identical(nrow(read_out(out)), 0L)
})

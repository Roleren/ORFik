store_fixture <- function(seed, chrs = c("1", "2", "10", "X", "MT", "KI1"), n = 600, weights = c(6, 5, 4, 3, 1, 1)) {
  set.seed(seed)
  w <- weights[seq_along(chrs)]
  d <- data.table::data.table(seqnames = factor(sample(chrs, n, TRUE, prob = w), levels = chrs),
                              start = sample.int(400L, n, TRUE),
                              strand = factor(sample(c("+", "-", "*"), n, TRUE), levels = c("+", "-", "*")),
                              size = sample(26:30, n, TRUE), score = sample.int(6L, n, TRUE))
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand, size)]
  data.table::setkey(d, NULL)
  d[sample(.N)]   # inputs need not be sorted
}

store_write <- function(tables, names = paste0("lib", seq_along(tables))) {
  dir <- tempfile("store_in"); dir.create(dir)
  paths <- file.path(dir, paste0(names, ".ofst"))
  for (i in seq_along(tables)) fst::write_fst(tables[[i]], paths[i])
  paths
}

pooled_merge <- function(paths) suppressWarnings(suppressMessages(ofst_merge(paths, keep_all_scores = FALSE)))

new_store <- function() file.path(tempfile("ofst_store"), "store")

test_that("store export equals a full merge after build and add", {
  paths <- store_write(lapply(1:8, store_fixture))
  store <- new_store()
  # small batches, segments and partitions exercise the segment / pooled-partition paths
  ofst_store_build(store, paths[1:5], batch_size = 2, partition_min_rows = 150, max_contrib_rows = 500,
                   workers = 2, verbose = FALSE)
  m <- ofst_store_manifest(store)
  expect_true("pooled" %in% m$parts$part)
  expect_gt(max(m$segments$seg), 1L)
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], pooled_merge(paths[1:5])[])
  extra <- store_write(list(store_fixture(9, chrs = c("2", "Y", "KI9"))), "lib9")
  ofst_store_add(store, c(paths[6:8], extra), workers = 1, verbose = FALSE)
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], pooled_merge(c(paths[6:8], paths[1:5], extra))[])
  expect_true(ofst_store_verify(store))
  # written file equals mergeLibs output format
  out <- tempfile(fileext = ".ofst")
  ofst_store_export_all(store, out, verbose = FALSE)
  expect_identical(fst::read_fst(out, as.data.table = TRUE)[], pooled_merge(c(paths, extra))[])
  expect_true(file.exists(paste0(out, ".removal_summary.rds")))
})

test_that("add then remove restores the store exactly", {
  paths <- store_write(lapply(1:6, store_fixture))
  store <- new_store()
  ofst_store_build(store, paths[1:4], partition_min_rows = 150, max_contrib_rows = 700, verbose = FALSE)
  before <- ofst_store_export_all(store, verbose = FALSE)
  mat_before <- ofst_store_read_matrix(store, as = "wide")
  ofst_store_add(store, paths[5:6], verbose = FALSE)
  expect_false(identical(ofst_store_export_all(store, verbose = FALSE)[], before[]))
  ofst_store_remove(store, c("lib5", "lib6"), verbose = FALSE)
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], before[])
  expect_identical(ofst_store_read_matrix(store, as = "wide")[], mat_before[])
  expect_true(ofst_store_verify(store))
  m <- ofst_store_manifest(store)
  expect_setequal(m$libs[state == "live"]$lib, paste0("lib", 1:4))
  expect_setequal(m$libs[state == "removed"]$lib, c("lib5", "lib6"))
  # compaction drops the removed rows and changes nothing else
  ofst_store_compact(store, min_dead_fraction = 0, verbose = FALSE)
  expect_false(any(ofst_store_manifest(store)$libs$state == "removed"))
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], before[])
  expect_true(ofst_store_verify(store))
  # removing every input leaves an empty track
  ofst_store_remove(store, paste0("lib", 1:4), verbose = FALSE)
  expect_identical(nrow(ofst_store_export_all(store, verbose = FALSE)), 0L)
  expect_error(ofst_store_remove(store, "lib1", verbose = FALSE), "not live")
})

test_that("replace swaps an input in one transaction; subsets export like a merge of them", {
  paths <- store_write(lapply(1:5, store_fixture))
  store <- new_store()
  ofst_store_build(store, paths, partition_min_rows = 150, verbose = FALSE)
  changed <- store_write(list(store_fixture(42, chrs = c("1", "3"))), "lib2")
  gen <- ofst_store_manifest(store)$generation
  ofst_store_replace(store, changed, verbose = FALSE)
  expect_identical(ofst_store_manifest(store)$generation, gen + 1L)
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], pooled_merge(c(paths[-2], changed))[])
  expect_identical(ofst_store_export_all(store, libs = c("lib1", "lib4"), verbose = FALSE)[],
                   pooled_merge(paths[c(1, 4)])[])
  expect_error(ofst_store_add(store, paths[1], verbose = FALSE), "already in the store")
  expect_true(ofst_store_verify(store))
})

test_that("in-RAM matrix mode reads, exports and writes back the same data", {
  paths <- store_write(lapply(1:4, store_fixture))
  store <- new_store()
  ofst_store_build(store, paths, partition_min_rows = 150, max_contrib_rows = 400, verbose = FALSE)
  x <- ofst_store_read_matrix(store)
  expect_s3_class(x, "ofst_store_matrix")
  expect_identical(ofst_store_export_all(x, verbose = FALSE)[], pooled_merge(paths)[])
  expect_identical(ofst_store_export_all(x, libs = "lib3", verbose = FALSE)[], pooled_merge(paths[3])[])
  copy <- new_store()
  ofst_store_write_matrix(x, copy, max_contrib_rows = 400)
  expect_true(ofst_store_verify(copy))
  expect_identical(ofst_store_export_all(copy, verbose = FALSE)[], pooled_merge(paths)[])
  expect_identical(ofst_store_read_matrix(copy, as = "wide")[], ofst_store_read_matrix(store, as = "wide")[])
  # the written store accepts further updates
  more <- store_write(list(store_fixture(77)), "lib77")
  ofst_store_add(copy, more, verbose = FALSE)
  expect_identical(ofst_store_export_all(copy, verbose = FALSE)[], pooled_merge(c(paths, more))[])
  # wide = per-library scores + total, same rows and values as the keep_all_scores merge
  wide <- ofst_store_read_matrix(store, as = "wide")
  ref <- suppressMessages(ofst_merge(paths, lib_names = paste0("lib", 1:4), keep_all_scores = TRUE))
  for (j in paste0("lib", 1:4)) ref[is.na(get(j)), (j) := 0L]
  data.table::setcolorder(ref, names(wide))
  expect_equal(as.data.frame(wide), as.data.frame(ref), ignore_attr = TRUE)
  # one chromosome only
  one <- ofst_store_read_matrix(store, seqnames = "MT", as = "wide")
  expect_identical(unique(as.character(one$seqnames)), "MT")
  expect_equal(sum(one$score), sum(pooled_merge(paths)[seqnames == "MT"]$score))
  skip_if_not_installed("Matrix")
  dg <- ofst_store_read_matrix(store, as = "dgCMatrix")
  m <- ofst_store_manifest(store)
  sums <- Reduce(`+`, lapply(dg, function(p) Matrix::colSums(p$m)))
  expect_equal(unname(sums[paste0("lib", 1:4)]), m$libs[, .(s = sum(score_sum)), by = lib][match(paste0("lib", 1:4), lib)]$s)
})

test_that("uncommitted rows from a crashed update are ignored", {
  paths <- store_write(lapply(1:3, store_fixture))
  store <- new_store()
  ofst_store_build(store, paths[1:2], partition_min_rows = 150, verbose = FALSE)
  m <- ofst_store_manifest(store)
  p <- m$parts$part[1]
  # orphan key rows, orphan contribution rows and a stray generation file
  keys <- fst::read_fst(file.path(store, p, "keys.fst"), from = 1, to = 3)
  keys$start <- keys$start + 100000L
  fst::write_fst(rbind(fst::read_fst(file.path(store, p, "keys.fst")), keys), file.path(store, p, "keys.fst"))
  seg <- m$segments[part == p][.N]
  sf <- file.path(store, p, sprintf("contrib_%04d.fst", seg$seg))
  fst::write_fst(rbind(fst::read_fst(sf), data.frame(row = 1L, score = 99)), sf)
  file.create(file.path(store, p, "total_g999999.fst"))
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], pooled_merge(paths[1:2])[])
  ofst_store_add(store, paths[3], verbose = FALSE)
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], pooled_merge(paths)[])
  expect_true(ofst_store_verify(store))
  expect_false(file.exists(file.path(store, p, "total_g999999.fst")))
  # a held lock blocks writers
  dir.create(file.path(store, ".lock"))
  expect_error(ofst_store_add(store, store_write(list(store_fixture(5)), "lib5"), verbose = FALSE), "locked")
})

test_that("totals above the int32 range are kept as doubles", {
  a <- store_fixture(1); b <- store_fixture(2)
  a[1, score := 2e9]; b <- rbind(b, a[1])[, .(score = sum(as.numeric(score))), by = .(seqnames, start, strand, size)]
  b[, score := as.integer(pmin(score, 2e9))]
  paths <- store_write(list(a, b))
  store <- new_store()
  ofst_store_build(store, paths, verbose = FALSE)
  got <- ofst_store_export_all(store, verbose = FALSE)
  expect_true(is.double(got$score))
  expect_equal(max(got$score), 4e9)
  expect_identical(got[], pooled_merge(paths)[])
  # an input whose own scores exceed int32 opens a double-score segment
  big <- store_fixture(3)[, score := as.numeric(score)]
  big[1, score := 3e9]
  big_path <- store_write(list(big), "big")
  ofst_store_add(store, big_path, verbose = FALSE)
  segs <- ofst_store_manifest(store)$segments
  expect_gt(nrow(segs), length(unique(segs$part)))
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], pooled_merge(c(paths, big_path))[])
  expect_true(ofst_store_verify(store))
  ofst_store_remove(store, "big", verbose = FALSE)
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], got[])
})

test_that("pre-sorted inputs (merge order) take the no-group path with the same result", {
  paths <- store_write(lapply(1:3, store_fixture))
  sorted <- tempfile(fileext = ".ofst")
  fst::write_fst(pooled_merge(paths[1:2]), sorted)   # sorted, unique keys
  d <- fst::read_fst(sorted, as.data.table = TRUE)
  # a merged file is in merge order: chromosome blocks in label order, strand "*" < "+" < "-"
  expect_true(ofst_keys_sorted_unique(list(as.integer(d$seqnames), d$start,
                                           .ofst_strand_to_sortkey[as.integer(d$strand)], d$size)))
  expect_false(ofst_keys_sorted_unique(list(c(1L, 1L), c(2L, 2L))))
  expect_true(ofst_keys_sorted_unique(list(c(1L, 1L, 2L), c(2L, 3L, 1L))))
  store <- new_store()
  ofst_store_build(store, c(sorted, paths[3]), libs = c("pool12", "lib3"), verbose = FALSE)
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], pooled_merge(paths)[])
})

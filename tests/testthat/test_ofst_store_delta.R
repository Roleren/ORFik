# Delta updates of the OFST store: appended key blocks and pending contributions instead of
# rewriting the partition's order / sorted key / total files.

delta_fixture <- function(seed, chrs = c("1", "2", "X", "KI1"), n = 800) {
  set.seed(seed)
  d <- data.table::data.table(seqnames = factor(sample(chrs, n, TRUE, prob = c(5, 4, 3, 1)[seq_along(chrs)]), levels = chrs),
                              start = sample.int(500L, n, TRUE),
                              strand = factor(sample(c("+", "-"), n, TRUE), levels = c("+", "-", "*")),
                              size = sample(27:30, n, TRUE), score = sample.int(5L, n, TRUE))
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand, size)]
  data.table::setkey(d, NULL)
  d
}
delta_write <- function(tables, names = paste0("lib", seq_along(tables))) {
  dir <- tempfile("delta_in"); dir.create(dir)
  paths <- file.path(dir, paste0(names, ".ofst"))
  for (i in seq_along(tables)) fst::write_fst(tables[[i]], paths[i])
  paths
}
delta_merge <- function(paths) suppressWarnings(suppressMessages(ofst_merge(paths, keep_all_scores = FALSE)))
with_delta <- function(rows, keys, code, on = TRUE) {
  old <- options(ORFik.ofst_store_delta_rows = rows, ORFik.ofst_store_delta_keys = keys, ORFik.ofst_store_delta = on)
  on.exit(options(old))
  force(code)
}

test_that("ofst_store_locate gives the same row ids and new keys as the k-way kernel", {
  set.seed(1)
  mk <- function(n) { k <- unique(data.table::data.table(seq = sample(1:2, n, TRUE), start = sample(1e4L, n, TRUE),
                                                         size = sample(25:30, n, TRUE)))
                      data.table::setorder(k); as.list(k) }
  base <- mk(3000); blk <- mk(500)
  blk_keep <- !paste(blk$seq, blk$start, blk$size) %in% paste(base$seq, base$start, base$size)
  blk <- lapply(blk, `[`, blk_keep)
  nb <- length(base$seq); n1 <- length(blk$seq)
  base_ids <- sample(nb)                    # row ids of the base keys in key order
  blk_ids <- nb + seq_len(n1)
  inputs <- lapply(1:4, function(i) mk(800))
  for (th in c(1L, 3L)) {
    r <- ofst_store_locate(list(base, blk), list(base_ids, blk_ids), inputs, nb + n1, th)
    # reference: k-way on the merged stored keys
    all_keys <- Map(c, base, blk); ids <- c(base_ids, blk_ids)
    o <- do.call(order, unname(all_keys))
    kw <- ofst_store_kway(lapply(all_keys, `[`, o), ids[o], inputs, nb + n1, th, stored_sorted = TRUE)
    expect_identical(r$rows, kw$rows)
    expect_identical(unname(r$new_keys), unname(kw$new_keys))
  }
  expect_error(ofst_store_locate(list(base), list(base_ids), list(lapply(inputs[[1]], rev)), nb, 1L), "not sorted")
})

test_that("delta adds, removes and replaces equal a full merge and fold cleanly", {
  paths <- delta_write(lapply(1:9, delta_fixture))
  store <- file.path(tempfile("delta_store"), "s")
  with_delta(100, 100, {
    ofst_store_build(store, paths[1:4], partition_min_rows = 400, verbose = FALSE)
    m0 <- ofst_store_manifest(store)
    expect_true("pooled" %in% m0$parts$part)
    expect_identical(m0$parts$n_base, m0$parts$n_keys)
    base_files <- sort(list.files(store, "^(order|total|skeys)_g", recursive = TRUE))
    ofst_store_add(store, paths[5:6], verbose = FALSE)
    m <- ofst_store_manifest(store)
    expect_gt(nrow(m$pending), 0)
    expect_gt(nrow(m$blocks), 0)
    expect_true(all(m$parts$n_keys >= m$parts$n_base))
    # no base file was rewritten
    expect_identical(sort(list.files(store, "^(order|total|skeys)_g", recursive = TRUE)), base_files)
    expect_identical(ofst_store_export_all(store, verbose = FALSE)[], delta_merge(paths[1:6])[])
    ofst_store_remove(store, c("lib2", "lib5"), verbose = FALSE)
    expect_identical(sort(list.files(store, "^(order|total|skeys)_g", recursive = TRUE)), base_files)
    ofst_store_add(store, paths[7], verbose = FALSE)
    ofst_store_replace(store, paths[8:9], c("lib3", "lib6"), verbose = FALSE)
    live <- paths[c(1, 4, 7, 8, 9)]
    expect_identical(ofst_store_export_all(store, verbose = FALSE)[], delta_merge(live)[])
    if (ORFik:::.ofst_store_concat_available()) {
      out <- tempfile(fileext = ".ofst")
      ofst_store_export_all(store, out, verbose = FALSE)
      expect_identical(fst::read_fst(out, as.data.table = TRUE)[], delta_merge(live)[])
    }
    sub <- c("lib1", "lib7")
    expect_identical(ofst_store_export_all(store, libs = sub, verbose = FALSE)[], delta_merge(paths[c(1, 7)])[])
    expect_true(ofst_store_verify(store))
    wide <- ofst_store_read_matrix(store, as = "wide")
    expect_identical(wide$score, delta_merge(live)$score)
    # folding (compaction) clears the deltas and keeps everything else
    ofst_store_compact(store, verbose = FALSE)
    m <- ofst_store_manifest(store)
    expect_equal(nrow(m$pending), 0)
    expect_equal(nrow(m$blocks), 0)
    expect_identical(m$parts$n_base, m$parts$n_keys)
    expect_identical(ofst_store_export_all(store, verbose = FALSE)[], delta_merge(live)[])
    expect_identical(ofst_store_read_matrix(store, as = "wide")[], wide[])
    expect_true(ofst_store_verify(store))
  })
})

test_that("delta limits fold partitions and the result equals a store without deltas", {
  paths <- delta_write(lapply(11:18, delta_fixture))
  s1 <- file.path(tempfile("delta_store"), "s"); s2 <- file.path(tempfile("delta_store"), "s")
  steps <- function(store) {
    ofst_store_build(store, paths[1:3], partition_min_rows = 200, verbose = FALSE)
    ofst_store_add(store, paths[4], verbose = FALSE)
    ofst_store_add(store, paths[5:6], verbose = FALSE)
    ofst_store_remove(store, "lib1", verbose = FALSE)
    ofst_store_add(store, paths[7:8], verbose = FALSE)
    ofst_store_manifest(store)
  }
  m1 <- with_delta(0.6, 0.05, steps(s1))
  m2 <- with_delta(1, 1, steps(s2), on = FALSE)
  expect_equal(nrow(m2$pending), 0)
  expect_identical(m2$parts$n_base, m2$parts$n_keys)
  expect_true(any(m1$parts$gen > 1L))       # some partitions were folded by the limits
  expect_identical(ofst_store_export_all(s1, verbose = FALSE)[], ofst_store_export_all(s2, verbose = FALSE)[])
  expect_identical(ofst_store_export_all(s1, verbose = FALSE)[], delta_merge(paths[2:8])[])
  expect_identical(ofst_store_read_matrix(s1, as = "wide")[], ofst_store_read_matrix(s2, as = "wide")[])
  expect_true(ofst_store_verify(s1))
  # a full update (deltas off) on a store with pending deltas folds them
  with_delta(1, 1, ofst_store_remove(s1, "lib2", verbose = FALSE), on = FALSE)
  m <- ofst_store_manifest(s1)
  expect_equal(nrow(m$pending), 0)
  expect_identical(ofst_store_export_all(s1, verbose = FALSE)[], delta_merge(paths[3:8])[])
  expect_true(ofst_store_verify(s1))
})

test_that("stores written before delta updates are read and updated", {
  paths <- delta_write(lapply(21:24, delta_fixture))
  store <- file.path(tempfile("delta_store"), "s")
  ofst_store_build(store, paths[1:3], partition_min_rows = 200, verbose = FALSE)
  m <- readRDS(file.path(store, "manifest.rds"))
  m$pending <- NULL; m$blocks <- NULL; m$parts$n_base <- NULL   # manifest of an older version
  saveRDS(m, file.path(store, "manifest.rds"))
  expect_identical(ofst_store_manifest(store)$parts$n_base, ofst_store_manifest(store)$parts$n_keys)
  with_delta(100, 100, ofst_store_add(store, paths[4], verbose = FALSE))
  expect_gt(nrow(ofst_store_manifest(store)$pending), 0)
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], delta_merge(paths)[])
  expect_true(ofst_store_verify(store))
})

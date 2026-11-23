# Sidecar chromosome-block index <file>.ofst.idx
idx_fixture <- function(seed = 1, chrs = c("1", "10", "2", "MT", "X"), n = 500) {
  set.seed(seed)
  d <- data.table::data.table(seqnames = factor(sample(chrs, n, TRUE), levels = chrs),
                              start = sample.int(300L, n, TRUE),
                              strand = factor(sample(c("+", "-", "*"), n, TRUE), levels = c("+", "-", "*")),
                              size = sample(26:30, n, TRUE), score = sample.int(5L, n, TRUE))
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand, size)]
  data.table::setkey(d, NULL)
  d
}

scan_blocks <- function(path) {
  old <- options(ORFik.ofst_index = FALSE); on.exit(options(old))
  .ofst_chr_blocks(path)
}

test_that("export.ofst writes a valid index that matches the seqnames scan", {
  d <- idx_fixture()
  path <- tempfile(fileext = ".ofst")
  export.ofst(as.data.frame(d), path)
  expect_true(file.exists(paste0(path, ".idx")))
  idx <- .ofst_index_read(path)
  expect_false(is.null(idx))
  expect_identical(idx$sorted_by, c("seqnames", "start"))
  expect_equal(idx$blocks$score, d[, sum(score), by = seqnames]$V1)
  expect_identical(.ofst_chr_blocks(path), scan_blocks(path))
  # GRanges method, and import is unaffected
  gr <- GRanges(d$seqnames, IRanges(d$start, width = 1L), d$strand, size = d$size, score = d$score)
  path_gr <- tempfile(fileext = ".ofst")
  export.ofst(gr, path_gr)
  expect_identical(.ofst_chr_blocks(path_gr), scan_blocks(path_gr))
  expect_identical(length(fimport(path_gr)), length(gr))
})

test_that("a stale or broken index is ignored and readers fall back to the scan", {
  d <- idx_fixture(2)
  path <- tempfile(fileext = ".ofst")
  export.ofst(as.data.frame(d), path)
  idx_copy <- tempfile(); file.copy(paste0(path, ".idx"), idx_copy)
  # file rewritten without the index (other writer): row count / size differ
  fst::write_fst(d[1:100], path)
  file.copy(idx_copy, paste0(path, ".idx"), overwrite = TRUE)
  expect_null(.ofst_index_read(path))
  expect_identical(.ofst_chr_blocks(path), scan_blocks(path))
  # garbage sidecar
  writeLines("not an index", paste0(path, ".idx"))
  expect_null(.ofst_index_read(path))
  expect_identical(.ofst_chr_blocks(path), scan_blocks(path))
  # tampered blocks
  export.ofst(as.data.frame(d), path)
  idx <- readRDS(paste0(path, ".idx")); idx$blocks$to[1] <- idx$blocks$to[1] + 1
  saveRDS(idx, paste0(path, ".idx"))
  expect_null(.ofst_index_read(path))
})

test_that("index of a split chromosome gives NULL blocks, like the scan", {
  d <- idx_fixture(3)
  d <- rbind(d[seqnames == "1"][1:5], d[seqnames != "1"], d[seqnames == "1"][-(1:5)])
  path <- tempfile(fileext = ".ofst")
  export.ofst(as.data.frame(d), path)
  idx <- .ofst_index_read(path)
  expect_false(is.null(idx))
  expect_identical(idx$sorted_by, character())
  expect_null(.ofst_chr_blocks(path))
  expect_null(scan_blocks(path))
})

test_that("pooled merge output gets an index; index option and rename helper", {
  tabs <- list(idx_fixture(4), idx_fixture(5))
  paths <- vapply(1:2, function(i) tempfile(fileext = ".ofst"), "")
  for (i in 1:2) fst::write_fst(tabs[[i]], paths[i])
  dt <- ofst_merge(paths, keep_all_scores = FALSE)
  out <- file.path(tempfile("merged")); dir.create(out); out <- file.path(out, "all.ofst")
  .ofst_save_merge(dt, out)
  idx <- .ofst_index_read(out)
  expect_false(is.null(idx))
  expect_identical(.ofst_chr_blocks(out), scan_blocks(out))
  expect_equal(sum(idx$blocks$score), sum(fst::read_fst(out, "score")$score))
  # option off: no index written, an old one removed
  old <- options(ORFik.ofst_index = FALSE)
  export.ofst(as.data.frame(tabs[[1]]), out)
  options(old)
  expect_false(file.exists(paste0(out, ".idx")))
  # rename helper moves the index along with its file
  tmp <- tempfile(fileext = ".ofst"); final <- tempfile(fileext = ".ofst")
  export.ofst(as.data.frame(tabs[[1]]), tmp)
  file.rename(tmp, final); .ofst_index_rename(tmp, final)
  expect_false(file.exists(paste0(tmp, ".idx")))
  expect_false(is.null(.ofst_index_read(final)))
})

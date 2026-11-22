# ofst_chr_rle(): Rle runs per chromosome block, against the position-sum path it replaces.
chr_rle_fixture <- function(seed = 1, n = 4000, chrs = c("chr1", "chr2", "chr3", "chrM"),
                            score_type = "integer") {
  set.seed(seed)
  d <- data.table::data.table(seqnames = factor(sample(chrs, n, TRUE), levels = chrs),
                              start = sample.int(380L, n, TRUE),
                              strand = factor(sample(c("+", "-", "*"), n, TRUE, prob = c(.45, .45, .1)),
                                              levels = c("+", "-", "*")),
                              size = sample(c(20L, 26L, 28L, 31L), n, TRUE),
                              score = sample.int(9L, n, TRUE))
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand, size)]
  data.table::setkey(d, NULL)
  if (score_type == "double") d[, score := score + 0.25]
  d
}

chr_rle_inputs <- function(d, si) {
  sizes <- sort(unique(d$size))
  list(chr = match(levels(d$seqnames), seqlevels(si))[as.integer(d$seqnames)],
       start = d$start, strand = match(levels(d$strand), c("+", "-", "*"))[as.integer(d$strand)],
       size_index = match(d$size, sizes), n_size = length(sizes), score = d$score,
       integer_weights = is.integer(d$score))
}

old_groups <- function(x, si) {
  sums <- ofst_point_sums(x$chr, x$start, x$strand, x$size_index, x$n_size, x$score)
  lapply(sums, function(s) .rlelist_from_points(data.table::setDT(s), si, x$integer_weights))
}

new_groups <- function(x, si) .point_rlelists_by_chr(x$chr, x$start, x$strand, x$size_index, x$n_size, x$score,
                                                      si, x$integer_weights)

test_that("ofst_chr_rle groups are identical to the position-sum path", {
  si <- Seqinfo(c("chr1", "chr2", "chr3", "chrM", "chrUn"), c(400L, 390L, 385L, 381L, 50L))  # chrUn empty
  for (score_type in c("integer", "double")) {
    d <- chr_rle_fixture(score_type = score_type)
    x <- chr_rle_inputs(d, si)
    expect_identical(new_groups(x, si), old_groups(x, si))
  }
  # a single position per chromosome at both ends, one read length, chromosome order not seqinfo order
  d <- data.table::data.table(seqnames = factor(c("chrM", "chrM", "chr1"), levels = c("chrM", "chr1")),
                              start = c(1L, 381L, 400L), strand = factor(c("+", "-", "*"), levels = c("+", "-", "*")),
                              size = 30L, score = c(2L, 3L, 4L))
  x <- chr_rle_inputs(d, si)
  expect_identical(new_groups(x, si), old_groups(x, si))
})

test_that("ofst_chr_rle returns NULL for unsorted or split blocks, and handles NA and overflow", {
  si <- Seqinfo(c("chr1", "chr2"), c(100L, 100L))
  expect_null(ofst_chr_rle(c(5L, 3L), c(1L, 1L), c(1L, 1L), 1L, c(1L, 1L), 100L, 0, 2))
  expect_null(ofst_chr_rle(c(5L, 101L), c(1L, 1L), c(1L, 1L), 1L, c(1L, 1L), 100L, 0, 2))
  expect_null(.point_rlelists_by_chr(c(1L, 2L, 1L), c(1L, 2L, 3L), c(1L, 1L, 1L), c(1L, 1L, 1L), 1L,
                                     c(1L, 1L, 1L), si, TRUE))
  r <- ofst_chr_rle(c(2L, 2L, 7L), c(1L, 1L, 2L), c(1L, 1L, 1L), 1L, c(3L, NA, 4L), 10L, 0, 3)
  expect_identical(r[[1]], list(values = c(0L, 3L, 0L), lengths = c(1L, 1L, 8L)))
  expect_identical(r[[2]], list(values = c(0L, 4L, 0L), lengths = c(6L, 1L, 3L)))
  big <- .Machine$integer.max
  r <- ofst_chr_rle(c(2L, 2L), c(1L, 1L), c(1L, 1L), 1L, c(big, big), 10L, 0, 2)
  expect_true(is.double(r[[1]]$values))
  expect_equal(r[[1]]$values[2], 2 * big)
  expect_true(is.integer(r[[2]]$values))   # reverse strand: only zeros
  # a sub-range of rows (one chromosome block of a larger table)
  r <- ofst_chr_rle(c(9L, 1L, 4L, 4L), c(1L, 1L, 1L, 1L), c(1L, 1L, 1L, 1L), 1L, c(5L, 1L, 2L, 3L), 6L, 1, 4)
  expect_identical(r[[1]], list(values = c(1L, 0L, 5L, 0L), lengths = c(1L, 2L, 1L, 2L)))
})

test_that("point OFST covRle via ofst_chr_rle equals the fallback path on files", {
  si <- Seqinfo(c("chr1", "chr2", "chr3", "chrM"), c(400L, 390L, 385L, 381L))
  d <- chr_rle_fixture(3)
  path <- tempfile(fileext = ".ofst"); fst::write_fst(d, path)
  out1 <- tempfile(); out1m <- tempfile()
  expect_true(.convert_point_ofst_cov(path, out1, out1m, si, "qs"))
  # force the fallback: chromosome split over two blocks gives the old path
  d2 <- rbind(d[seqnames == "chr1"][1:10], d[seqnames != "chr1"], d[seqnames == "chr1"][-(1:10)])
  path2 <- tempfile(fileext = ".ofst"); fst::write_fst(d2, path2)
  out2 <- tempfile(); out2m <- tempfile()
  expect_true(.convert_point_ofst_cov(path2, out2, out2m, si, "qs"))
  expect_identical(read_RDSQS(paste0(out1, ".covqs")), read_RDSQS(paste0(out2, ".covqs")))
  expect_identical(read_RDSQS(paste0(out1m, ".covqs")), read_RDSQS(paste0(out2m, ".covqs")))
})

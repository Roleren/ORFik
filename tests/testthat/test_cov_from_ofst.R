point_ofst_fixture <- function(seed = 1, n = 3000, score_type = "integer", sorted = TRUE, chrs = c("chr1", "chr2")) {
  set.seed(seed)
  d <- data.table::data.table(seqnames = factor(sample(chrs, n, TRUE), levels = chrs),
                              start = sample.int(900L, n, TRUE),
                              strand = factor(sample(c("+", "-", "*"), n, TRUE, prob = c(.45, .45, .1)),
                                              levels = c("+", "-", "*")),
                              size = sample(c(26L, 28L, 31L), n, TRUE),
                              score = sample.int(9L, n, TRUE))
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand, size)]
  data.table::setkey(d, NULL)
  if (!sorted) d <- d[sample(.N)]
  if (score_type == "double") d[, score := score + 0.5]
  d
}

point_seqinfo <- function() Seqinfo(c("chr1", "chr2", "chr3"), c(1000L, 950L, 400L))   # chr3 empty

point_cov_reference <- function(d, si) {
  gr <- GRanges(d$seqnames, IRanges(d$start, width = 1L), strand = d$strand, seqinfo = si)
  mcols(gr) <- data.frame(size = d$size, score = d$score)
  list(list = covRleListFromGR(gr, weight = "score", verbose = FALSE),
       merged = {x <- covRleFromGR(gr, weight = "score"); seqinfo(x) <- si; x})
}

point_cov_new <- function(d, si) {
  path <- tempfile(fileext = ".ofst"); fst::write_fst(d, path)
  out <- tempfile(); out_m <- tempfile()
  ok <- .convert_point_ofst_cov(path, out, out_m, si, "qs")
  list(ok = ok, list = if (ok) read_RDSQS(paste0(out, ".covqs")), merged = if (ok) read_RDSQS(paste0(out_m, ".covqs")))
}

test_that("point OFST coverage equals covRleFromGR / covRleListFromGR", {
  si <- point_seqinfo()
  for (score_type in c("integer", "double")) for (sorted in c(TRUE, FALSE)) {
    d <- point_ofst_fixture(score_type = score_type, sorted = sorted)
    ref <- point_cov_reference(d, si)
    got <- point_cov_new(d, si)
    expect_true(got$ok)
    expect_identical(got$list, ref$list)
    expect_identical(got$merged, ref$merged)
  }
})

test_that("integer overflow at a position gives numeric coverage for that chromosome", {
  si <- point_seqinfo()
  d <- point_ofst_fixture(2)
  d[1:2, `:=`(seqnames = factor("chr1", levels = levels(seqnames)), start = 10L, strand = factor("+", levels = levels(strand)),
              size = 26L, score = .Machine$integer.max)]
  d <- d[, .(score = sum(as.numeric(score))), keyby = .(seqnames, start, strand, size)]
  d[, score := as.integer(pmin(score, .Machine$integer.max))]
  d <- rbind(d, d[seqnames == "chr1" & start == 10L & strand == "+" & size == 26L][, size := 28L])
  d <- d[order(seqnames, start, strand, size)]
  got <- point_cov_new(d, si)
  expect_true(is.numeric(S4Vectors::runValue(f(got$merged)[["chr1"]])))
  expect_equal(as.numeric(f(got$merged)[["chr1"]][10]), 2 * .Machine$integer.max)
  expect_true(is.integer(S4Vectors::runValue(f(got$merged)[["chr2"]])))
})

test_that("non-point files are left to the general converter", {
  si <- point_seqinfo()
  d <- point_ofst_fixture(3)
  path <- tempfile(fileext = ".ofst")
  fst::write_fst(cbind(d, cigar = "28M"), path)
  expect_false(.convert_point_ofst_cov(path, tempfile(), tempfile(), si, "qs"))
  fst::write_fst(cbind(d, width = 2L), path)
  expect_false(.convert_point_ofst_cov(path, tempfile(), tempfile(), si, "qs"))
  d2 <- data.table::copy(d)[1, start := 2000L]      # outside chr1
  fst::write_fst(d2, path)
  expect_false(.convert_point_ofst_cov(path, tempfile(), tempfile(), si, "qs"))
})

# ofst_derive(): single-pass covRleList, covRle, bigWig and count tables of a sorted point OFST
# file, against the separate steps (old covRle path, export.bigWig, countOverlapsW).

derive_fixture <- function(dir, seed = 1, n = 60000, score_type = "integer", with_index = TRUE,
                           chrs = c("chr2", "chr1", "chrM"), lens = c(chr1 = 30000L, chr2 = 25000L, chrM = 9000L)) {
  set.seed(seed)
  d <- data.table::data.table(seqnames = factor(sample(chrs, n, TRUE, prob = c(.5, .4, .1)), levels = chrs),
                              strand = factor(sample(c("+", "-", "*"), n, TRUE, prob = c(.45, .45, .1)),
                                              levels = c("+", "-", "*")),
                              size = sample(c(20L, 26L, 28L, 31L), n, TRUE),
                              score = sample.int(9L, n, TRUE))
  d[, start := sample.int(lens[as.character(seqnames[1])], .N, TRUE), by = seqnames]
  # merge/collapse order: chromosome label blocks, start, then strand and size inside a position
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand, size)]
  data.table::setkey(d, NULL)
  data.table::setcolorder(d, c("seqnames", "start", "strand", "size", "score"))
  if (score_type == "double") d[, score := score + 0.25]
  file <- file.path(dir, paste0("lib_", seed, "_", score_type, if (!with_index) "_noidx", ".ofst"))
  fst::write_fst(d, file)
  if (with_index) ORFik:::.ofst_index_write(file, d)
  file
}

derive_seqinfo <- function() Seqinfo(c("chr1", "chr2", "chrM", "chrUn"), c(30000L, 25000L, 9000L, 500L))

old_cov <- function(file, si, dir) {
  op <- options(ORFik.ofst_derive = FALSE); on.exit(options(op))
  o1 <- file.path(dir, "old_list"); o2 <- file.path(dir, "old_merged")
  stopifnot(ORFik:::.convert_point_ofst_cov(file, o1, o2, si, "qs"))
  list(list = read_RDSQS(paste0(o1, ".covqs")), merged = read_RDSQS(paste0(o2, ".covqs")))
}

test_that("ofst_derive covRleList and covRle are identical to the old path", {
  dir <- tempfile("derive_"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  si <- derive_seqinfo()
  for (score_type in c("integer", "double")) for (idx in c(TRUE, FALSE)) {
    file <- derive_fixture(dir, score_type = score_type, with_index = idx)
    ref <- old_cov(file, si, dir)
    # one wave and one task per block; several waves; small tasks inside blocks, several threads
    for (cfg in list(list(wave = 1e9, task = 1e9, th = 1L), list(wave = 1, task = 1e4, th = 3L),
                     list(wave = 3e4, task = 1e4, th = 2L))) {
      covl <- file.path(dir, "new_list.covqs"); covm <- file.path(dir, "new_merged.covqs")
      res <- ofst_derive(file, si, covlist_file = covl, cov_file = covm, threads = cfg$th,
                         wave_rows = cfg$wave, task_rows = cfg$task)
      expect_false(isFALSE(res))
      expect_identical(res$sizes, c(20L, 26L, 28L, 31L))
      expect_identical(read_RDSQS(covl), ref$list)
      expect_identical(read_RDSQS(covm), ref$merged)
      unlink(c(covl, covm))
    }
  }
})

test_that(".convert_point_ofst_cov uses ofst_derive and matches the old path; no size column works", {
  dir <- tempfile("derive_"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  si <- derive_seqinfo()
  file <- derive_fixture(dir, seed = 2)
  ref <- old_cov(file, si, dir)
  expect_true(ORFik:::.convert_point_ofst_cov(file, file.path(dir, "a"), file.path(dir, "b"), si, "qs"))
  expect_identical(read_RDSQS(file.path(dir, "a.covqs")), ref$list)
  expect_identical(read_RDSQS(file.path(dir, "b.covqs")), ref$merged)
  d <- fst::read_fst(file, as.data.table = TRUE)[, size := NULL]
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand)]
  data.table::setkey(d, NULL)
  f2 <- file.path(dir, "nosize.ofst"); fst::write_fst(d, f2)
  ref2 <- old_cov(f2, si, dir)
  covl <- file.path(dir, "n_list.covqs"); covm <- file.path(dir, "n_merged.covqs")
  expect_false(isFALSE(ofst_derive(f2, si, covl, covm, task_rows = 1e4, threads = 2L)))
  expect_identical(read_RDSQS(covl), ref2$list)
  expect_identical(read_RDSQS(covm), ref2$merged)
})

test_that("ofst_derive bigWig files import as export.bigWig files do", {
  dir <- tempfile("derive_"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  si <- derive_seqinfo()
  file <- derive_fixture(dir, seed = 3)
  covm <- file.path(dir, "m.covqs")
  expect_false(isFALSE(ofst_derive(file, si, cov_file = covm, bigwig_file = file.path(dir, "new"),
                                   threads = 2L, bigwig_threads = 2L, task_rows = 1e4, wave_rows = 2e4)))
  cov <- read_RDSQS(covm)
  export.bigWig(cov, file.path(dir, "rt.bigWig"), writer = "rtracklayer")
  export.bigWig(cov, file.path(dir, "orfik.bigWig"), writer = "ORFik", threads = 2L)
  for (s in c("forward", "reverse")) {
    new <- file.path(dir, paste0("new_", s, ".bigWig"))
    rt <- file.path(dir, paste0("rt_", s, ".bigWig"))
    expect_identical(rtracklayer::import.bw(new, as = "RleList"), rtracklayer::import.bw(rt, as = "RleList"))
    expect_identical(rtracklayer::import.bw(new), rtracklayer::import.bw(rt))
    expect_identical(seqinfo(rtracklayer::BigWigFile(new)), seqinfo(rtracklayer::BigWigFile(rt)))
    expect_identical(rtracklayer::import.bw(new, as = "RleList"),
                     rtracklayer::import.bw(file.path(dir, paste0("orfik_", s, ".bigWig")), as = "RleList"))
  }
  expect_length(list.files(dir, pattern = "^\\.new"), 0)   # no temporary files left
})

test_that("count tables from ofst_derive equal countOverlapsW on the covRle", {
  dir <- tempfile("derive_"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  si <- derive_seqinfo()
  file <- derive_fixture(dir, seed = 4)
  set.seed(5)
  make_tx <- function(n, strands = c("+", "-", "*")) {   # 1-3 exons, exons may overlap inside a transcript
    k <- sample.int(3L, n, TRUE)
    chr <- sample(c("chr1", "chr2", "chrM"), n, TRUE)
    lens <- c(chr1 = 30000L, chr2 = 25000L, chrM = 9000L)
    st <- rep(sample(strands, n, TRUE, prob = c(.45, .45, .1)[seq_along(strands)]), k)
    seqn <- rep(chr, k)
    s <- unlist(lapply(seq_len(n), function(i) sort(sample.int(lens[chr[i]] - 200L, k[i]))))
    w <- sample(0:150, length(s), TRUE)
    gr <- GRanges(seqn, IRanges(s, width = w), strand = st, seqinfo = si)
    grl <- split(gr, rep(seq_len(n), k))
    names(grl) <- paste0("tx", seq_len(n))
    grl
  }
  # "*" exons count 0 in the large-region path; the small-region path does not accept them at all
  regions <- list(big = make_tx(6000), small = make_tx(300, c("+", "-")))
  covm <- file.path(dir, "m.covqs")
  res <- ofst_derive(file, si, cov_file = covm, regions = regions, threads = 2L)
  cov <- read_RDSQS(covm)
  for (nm in names(regions)) expect_identical(res$counts[[nm]], countOverlapsW(regions[[nm]], cov, weight = NULL))
  expect_type(res$counts$big, "integer")
})

test_that("ofst_derive returns FALSE and writes nothing for data it cannot take", {
  dir <- tempfile("derive_"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  si <- derive_seqinfo()
  file <- derive_fixture(dir, seed = 6)
  d <- fst::read_fst(file, as.data.table = TRUE)
  check <- function(dd, name, seq_info = si) {
    f <- file.path(dir, paste0(name, ".ofst")); fst::write_fst(dd, f)
    out <- file.path(dir, paste0(name, c("_l.covqs", "_m.covqs")))
    expect_false(ofst_derive(f, seq_info, out[1], out[2], bigwig_file = file.path(dir, name), task_rows = 1e4))
    expect_false(any(file.exists(out)))
    expect_length(list.files(dir, pattern = name, all.files = TRUE), 1L)   # only the input
  }
  u <- data.table::copy(d); u[c(5, 9), start := start[c(9, 5)]]; check(u, "unsorted")
  u <- data.table::copy(d); u[100, score := NA]; check(u, "na_score")
  u <- data.table::copy(d); u[nrow(u), start := 9001L]; check(u, "out_of_range")
  check(d, "no_chr", Seqinfo(c("chr1", "chr2"), c(30000L, 25000L)))
  u <- data.table::copy(d); u[, width := 1L]; u[7, width := 2L]; check(u, "width")
  u <- data.table::copy(d)[, seqnames := factor(as.character(seqnames), levels = c("chr1", "chr2", "chrM"))]
  u <- rbind(u[seqnames == "chr2"][1:10], u[seqnames != "chr2"], u[seqnames == "chr2"][-(1:10)])
  check(u, "split_chr")
})

test_that("bigWig stream: any add order, same content; id order gives the bytes of bw_write_runs", {
  dir <- tempfile("bwstream_"); dir.create(dir); on.exit(unlink(dir, recursive = TRUE))
  set.seed(7)
  chrs <- c(sprintf("s%03d", sample.int(300L)))   # 300 sequences, not sorted (above 256: one leaf)
  lens <- sample(100:3000, length(chrs), TRUE)
  rl <- lapply(lens, function(l) { p <- sort(sample.int(l, min(l, 40))); x <- integer(l); x[p] <- sample.int(5, length(p), TRUE); Rle(x) })
  names(rl) <- chrs
  rl <- IRanges::RleList(rl, compress = FALSE)
  a <- file.path(dir, "a.bigWig"); b <- file.path(dir, "b.bigWig"); s <- file.path(dir, "s.bigWig")
  bw_write_runs(a, chrs, lens, lapply(as.list(rl), runValue), lapply(as.list(rl), runLength), threads = 2L)
  for (async in c(FALSE, TRUE)) {
    w <- bw_stream_open(b, chrs, lens, threads = 2L, async = async)
    for (i in seq_along(chrs)) bw_stream_add(w, i, runValue(rl[[i]]), runLength(rl[[i]]))
    bw_stream_close(w)
    expect_identical(unname(tools::md5sum(b)), unname(tools::md5sum(a)))
  }
  w <- bw_stream_open(s, chrs, lens, threads = 3L)
  for (i in rev(seq_along(chrs)[-1])) bw_stream_add(w, i, runValue(rl[[i]]), runLength(rl[[i]]))
  expect_error(bw_stream_add(w, 2L, runValue(rl[[2]]), runLength(rl[[2]])), "twice")
  bw_stream_close(w)   # sequence 1 never added: written without data
  ref <- rtracklayer::import.bw(a, as = "RleList")
  got <- rtracklayer::import.bw(s, as = "RleList")
  expect_identical(seqinfo(rtracklayer::BigWigFile(s)), seqinfo(rtracklayer::BigWigFile(a)))
  keep <- names(ref) != chrs[1]
  expect_identical(got[keep], ref[keep])
  expect_identical(as.vector(got[[chrs[1]]]), numeric(lens[1]))
  expect_identical(as.vector(rtracklayer::import.bw(s, which = GRanges(chrs[5], IRanges(1, lens[5])), as = "NumericList")[[1]]),
                   as.numeric(rl[[5]]))
  expect_error(bw_stream_close(w), "closed")
  w <- bw_stream_open(file.path(dir, "x.bigWig"), chrs, lens); bw_stream_abort(w)
  expect_error(bw_stream_add(w, 1L, 0L, lens[1]), "closed")
})

test_that("rle_range_sums equals sums of views", {
  set.seed(8)
  x <- Rle(sample(c(0L, 0L, 1L, 7L), 5000, TRUE))
  s <- sample.int(5000, 300, TRUE); w <- sample(0:200, 300, TRUE); e <- pmin(s + w - 1L, 5000L)
  expect_equal(rle_range_sums(runValue(x), runLength(x), s, e), as.numeric(sum(Views(x, IRanges(s, e)))))
  expect_null(rle_range_sums(runValue(x), runLength(x), 1L, 5001L))
})

test_that("count fast path follows countOverlapsW types: double strands, integer overflow, fallbacks", {
  set.seed(9)
  si <- Seqinfo(c("chr1", "chr2"), c(5000L, 4000L))
  mk <- function(l, big = FALSE) { x <- sample(0:4, l, TRUE); if (big) { x <- as.numeric(x); x[100] <- 3e9 }; Rle(x) }
  fw <- IRanges::RleList(chr1 = mk(5000, big = TRUE), chr2 = mk(4000), compress = FALSE)   # chr1 forward double
  rv <- IRanges::RleList(chr1 = mk(5000), chr2 = mk(4000), compress = FALSE)
  seqinfo(fw) <- seqinfo(rv) <- si
  cov <- covRle(fw, rv)
  n <- 6000
  tx_of <- function(chrs, strands) {
    gr <- GRanges(sample(chrs, n, TRUE), IRanges(sample.int(3800, n, TRUE), width = sample(1:150, n, TRUE)),
                  strand = sample(strands, n, TRUE), seqinfo = si)
    grl <- split(gr, seq_len(n)); names(grl) <- paste0("t", seq_len(n)); grl
  }
  regions <- list(dbl = tx_of(c("chr1", "chr2"), c("+", "-")),      # touches the double strand
                  int = tx_of("chr2", c("+", "-", "*")),           # integer only
                  rev = tx_of("chr1", "-"))                         # chr1 reverse only: integer
  got <- ORFik:::.cov_regions_counts(regions, cov)
  ref <- lapply(regions, countOverlapsW, subject = cov, weight = NULL)
  expect_identical(got, ref)
  expect_type(got$dbl, "double"); expect_type(got$int, "integer"); expect_type(got$rev, "integer")
  expect_false(any(vapply(ORFik:::.cov_regions_counts_fast(regions, cov), is.null, logical(1))))
})

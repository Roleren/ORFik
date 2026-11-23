# ORFik's bigWig writer (src/bigwig_writer.cpp) against rtracklayer::export.bw.
bw_fixture <- function(seed = 1, integer = TRUE) {
  set.seed(seed)
  chrs <- c("chr2", "chr10", "1", "chrM", "chrX", "chr1")     # not in bytewise order
  lens <- c(5000L, 3000L, 2500L, 300L, 1200L, 4000L)
  one <- function(n) {
    v <- integer(n); hit <- sample.int(n, max(1, n %/% 7)); v[hit] <- sample.int(50L, length(hit), TRUE)
    if (!integer) v <- v * 0.37
    S4Vectors::Rle(v)
  }
  f <- IRanges::RleList(setNames(lapply(lens, one), chrs), compress = FALSE)
  r <- IRanges::RleList(setNames(lapply(lens, one), chrs), compress = FALSE)
  r[["chrM"]] <- S4Vectors::Rle(if (integer) 0L else 0, 300L)      # an all-zero sequence
  si <- Seqinfo(chrs, lens)
  seqinfo(f) <- seqinfo(r) <- si
  covRle(f, r)
}

bw_imports <- function(file) {
  list(gr = rtracklayer::import.bw(file), rle = rtracklayer::import.bw(file, as = "RleList"),
       nl = rtracklayer::import.bw(file, as = "NumericList"))
}

test_that("ORFik bigWig writer imports identical to rtracklayer's files", {
  for (integer in c(TRUE, FALSE)) {
    x <- bw_fixture(integer = integer)
    a <- file.path(tempdir(), paste0("bw_rt_", integer)); b <- file.path(tempdir(), paste0("bw_orfik_", integer))
    export.bigWig(x, a, writer = "rtracklayer")
    export.bigWig(x, b, writer = "ORFik", threads = 2)
    for (s in c("_forward.bigWig", "_reverse.bigWig")) {
      expect_identical(bw_imports(paste0(b, s)), bw_imports(paste0(a, s)))
      expect_identical(seqinfo(rtracklayer::BigWigFile(paste0(b, s))),
                       seqinfo(rtracklayer::BigWigFile(paste0(a, s))))
    }
    # values as written (float precision)
    expect_equal(as.list(rtracklayer::import.bw(paste0(b, "_forward.bigWig"), as = "RleList")[seqlevels(x)]),
                 lapply(as.list(f(x)), function(r) S4Vectors::Rle(as.numeric(r))), tolerance = 1e-6)
  }
})

test_that("ORFik bigWig writer: single-file mode, threads, level, summaries", {
  x <- bw_fixture(2)
  a <- file.path(tempdir(), "bw_single_rt"); b <- file.path(tempdir(), "bw_single_orfik")
  export.bigWig(x, a, split.by.strand = FALSE, writer = "rtracklayer")
  export.bigWig(x, b, split.by.strand = FALSE, writer = "ORFik", threads = 1)
  expect_identical(bw_imports(paste0(b, ".bigWig")), bw_imports(paste0(a, ".bigWig")))
  # threads do not change the bytes; level 1 changes only the compression
  b4 <- file.path(tempdir(), "bw_single_orfik4")
  export.bigWig(x, b4, split.by.strand = FALSE, threads = 4)
  expect_identical(unname(tools::md5sum(paste0(b4, ".bigWig"))), unname(tools::md5sum(paste0(b, ".bigWig"))))
  b1 <- file.path(tempdir(), "bw_single_orfik_l1")
  export.bigWig(x, b1, split.by.strand = FALSE, compression_level = 1)
  expect_identical(bw_imports(paste0(b1, ".bigWig")), bw_imports(paste0(b, ".bigWig")))
  # whole-sequence and zoomed summaries (mean over all bases, max) are exact
  cov <- f(x) + r(x)
  bf <- rtracklayer::BigWigFile(paste0(b, ".bigWig"))
  w <- GRanges(names(cov), IRanges(1, lengths(cov)))
  expect_equal(unlist(lapply(summary(bf, w, type = "mean"), function(g) g$score)),
               unname(vapply(as.list(cov), function(r) mean(as.numeric(r)), 0)), tolerance = 1e-5)
  expect_equal(unlist(lapply(summary(bf, w, type = "max"), function(g) g$score)),
               unname(vapply(as.list(cov), function(r) max(as.numeric(r)), 0)), tolerance = 1e-6)
  w2 <- GRanges("chr2", IRanges(1, 4096)); n2 <- as.numeric(cov[["chr2"]][1:4096])
  expect_equal(summary(bf, w2, size = 4L, type = "max")[[1]]$score, apply(matrix(n2, ncol = 4), 2, max))
  expect_equal(summary(bf, w2, size = 4L, type = "mean")[[1]]$score, colMeans(matrix(n2, ncol = 4)), tolerance = 1e-5)
})

test_that("ORFik bigWig writer: multi-level chromosome tree and input checks", {
  n <- 40L
  chrs <- sprintf("scaffold_%03d", seq_len(n)); lens <- rep(120L, n)
  set.seed(3)
  rl <- IRanges::RleList(setNames(lapply(lens, function(l) S4Vectors::Rle(sample(0:3, l, TRUE))), chrs), compress = FALSE)
  file <- tempfile(fileext = ".bigWig")
  ord <- order(chrs, method = "radix")
  # block size 3: 40 sequences need a 4-level B+ tree (the default is one leaf of up to 65,535)
  bw_write_runs(file, chrs[ord], lens[ord], lapply(as.list(rl)[ord], S4Vectors::runValue),
                lapply(as.list(rl)[ord], S4Vectors::runLength), threads = 2L, chrom_block_size = 3L)
  got <- rtracklayer::import.bw(file, as = "RleList")
  expect_identical(names(got), chrs)
  expect_equal(as.list(got), lapply(as.list(rl), function(r) S4Vectors::Rle(as.numeric(r))))
  one <- rtracklayer::import.bw(file, which = GRanges("scaffold_017", IRanges(1, 120)), as = "NumericList")[[1]]
  expect_equal(one, as.numeric(rl[["scaffold_017"]]))
  # inputs
  expect_error(bw_write_runs(file, c("b", "a", "c"), c(1L, 1L, 1L), list(0L, 0L, 0L), list(1L, 1L, 1L),
                             chrom_block_size = 2L), "sorted")
  expect_error(bw_write_runs(file, c("b", "b"), c(1L, 1L), list(0L, 0L), list(1L, 1L)), "unique")
  expect_error(bw_write_runs(file, "a", 5L, list(1L), list(4L)), "sequence length")
  expect_error(bw_write_runs(file, "a", 2L, list(c(1L, NA)), list(c(1L, 1L))), "NA")
})

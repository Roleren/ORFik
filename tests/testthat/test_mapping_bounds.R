library(ORFik)

test_that("native to-transcript mapping stays inside its paired transcript", {
  for (direction in c("+", "-")) {
    x <- ORFik:::pmapToTranscriptsCPP(c(205L, 999L, 100L), c(210L, 1000L, 120L),
      c(100L, 200L), c(120L, 220L), c(1L, 2L), direction, c(0L, 1L, 0L))
    expect_identical(x$ranges, list(c(0L, 0L, 1L), c(-1L, -1L, 21L)))
  }
  expect_identical(ORFik:::pmapToTranscriptsCPP(integer(), integer(), integer(),
    integer(), integer(), "+", integer())$ranges, list(integer(), integer()))
  expect_identical(ORFik:::pmapToTranscriptsCPP(1L, 2L, integer(), integer(),
    integer(), "+", 0L)$ranges, list(0L, -1L))
})

test_that("native mapping validates coordinates and vector lengths", {
  expect_error(ORFik:::pmapToTranscriptsCPP(1L, 2L, 1L, 10L, integer(), "+", 0L), "lengths")
  expect_error(ORFik:::pmapToTranscriptsCPP(1L, 2L, 1L, 10L, 1L, "x", 0L), "direction")
  expect_error(ORFik:::pmapToTranscriptsCPP(1L, 2L, NA_integer_, 10L, 1L, "+", 0L), "exon")
  expect_error(ORFik:::pmapToTranscriptsCPP(1L, 2L, c(1L, 20L), c(10L, 30L),
    c(1L, 1L), "+", 1L), "first exon")
  expect_error(ORFik:::pmapFromTranscriptsCPP(1L, 22L, 100L, 120L, 1L, "+", TRUE), "out of bounds")
  expect_error(ORFik:::pmapFromTranscriptsCPP(0L, 2L, 100L, 120L, 1L, "+", TRUE), "out of bounds")
  expect_error(ORFik:::pmapFromTranscriptsCPP(NA_integer_, 2L, 100L, 120L, 1L, "+", TRUE), "out of bounds")
  expect_error(ORFik:::pmapFromTranscriptsCPP(c(1L, 1L), c(2L, 2L), 100L, 120L,
    1L, "+", TRUE), "group counts")
  expect_error(ORFik:::pmapToTranscriptsCPP(1L, 2L, c(1L, 1L),
    c(.Machine$integer.max, 1L), c(1L, 1L), "+", 0L), "integer range")
  invalid <- ORFik:::pmapToTranscriptsCPP(c(0L, NA_integer_, 105L),
    c(1L, NA_integer_, 104L), 100L, 120L, 1L, "+", c(0L, 0L, 0L))
  expect_identical(invalid$ranges, list(integer(3L), rep(-1L, 3L)))
})

test_that("large genomic positions do not overflow transcript coordinates", {
  begin <- .Machine$integer.max - 30L
  for (direction in c("+", "-")) {
    res <- ORFik:::pmapToTranscriptsCPP(begin + 2L, begin + 5L, begin,
      begin + 20L, 1L, direction, 0L)
    expect_identical(res$ranges, if (direction == "+") list(3L, 6L) else list(16L, 19L))
  }
})

test_that("irregular exon ordering uses a bounded linear fallback", {
  s <- c(300L, 100L, 200L, seq(400L, 900L, 100L))
  res <- ORFik:::pmapToTranscriptsCPP(205L, 210L, s, s + 20L,
    rep(1L, length(s)), "+", 0L)
  expect_identical(res$ranges, list(48L, 53L))
})

test_that("single-exon GRanges references recycle and map each paired row", {
  tx <- GRanges("1", IRanges(c(100, 200, 300), width = 21), "+")
  q <- GRanges("1", IRanges(c(102, 202, 302), width = 4), "+")
  res <- pmapToTranscriptF(q, tx)
  expect_identical(start(res), rep(3L, 3L))
  expect_identical(end(res), rep(6L, 3L))
})

test_that("genomic sequence identity is checked within each pair", {
  tx <- GRangesList(a = GRanges("1", IRanges(100, 120), "+"),
                    b = GRanges("1", IRanges(100, 120), "+"))
  seqlevels(tx) <- c("1", "2")
  q <- GRangesList(a = GRanges("1", IRanges(105, 110), "+"),
                   b = GRanges("2", IRanges(105, 110), "+"))
  res <- pmapToTranscriptF(q, tx)
  expect_identical(start(res[[1]]), 6L)
  expect_identical(start(res[[2]]), 0L)
  expect_identical(end(res[[2]]), -1L)
  expect_identical(as.character(strand(res[[2]])), "*")
})

test_that("empty inputs and named native indices retain their schemas", {
  tx <- GRangesList(tx = GRanges("1", IRanges(100, 120), "+"))
  expect_identical(pmapToTranscriptF(GRangesList(), tx), GRangesList())
  expect_length(pmapToTranscriptF(GRangesList(empty = GRanges()), tx)[[1]], 0L)
  expect_identical(pmapFromTranscriptF(IRanges(), tx), tx[integer()])
  idx <- c(first = 1L)
  expect_identical(ORFik:::pmapToTranscriptsCPP(105L, 110L, 100L, 120L, idx, "+", 0L)$index, 1L)
  expect_identical(ORFik:::pmapFromTranscriptsCPP(2L, 5L, 100L, 120L, idx, "+", FALSE)$index, 1L)
})

test_that("both mapping directions match independent exon-by-exon coordinates", {
  set.seed(510)
  for (direction in c("+", "-")) for (iteration in seq_len(30L)) {
    n <- sample.int(25L, 1L)
    widths <- sample.int(50L, n, replace = TRUE)
    s <- cumsum(c(100L, head(widths, -1L) + sample.int(40L, n - 1L)))
    en <- s + widths - 1L
    if (direction == "-") { s <- rev(s); en <- rev(en); widths <- rev(widths) }
    pick <- sample.int(n, 20L, replace = TRUE)
    qs <- s[pick] + sample.int(10L, 20L, replace = TRUE) %% widths[pick]
    qe <- pmin(qs + sample.int(10L, 20L, replace = TRUE) - 1L, en[pick])
    prefix <- c(0L, head(cumsum(widths), -1L))[pick]
    expected <- if (direction == "+") list(prefix + (qs - s[pick] + 1L),
      prefix + (qe - s[pick] + 1L)) else list(prefix + (en[pick] - qe + 1L),
      prefix + (en[pick] - qs + 1L))
    res <- ORFik:::pmapToTranscriptsCPP(qs, qe, s, en, rep(1L, n), direction, rep(0L, 20L))
    expect_identical(res$ranges, expected)
    first <- sample.int(sum(widths), 1L)
    last <- first + sample.int(sum(widths) - first + 1L, 1L) - 1L
    bounds <- c(0L, cumsum(widths)); os <- oe <- integer(n)
    for (j in seq_len(n)) {
      use <- intersect(first:last, seq.int(bounds[j] + 1L, bounds[j + 1L]))
      if (length(use)) {
        pos <- seq.int(s[j], en[j]); if (direction == "-") pos <- rev(pos)
        chosen <- pos[use - bounds[j]]; os[j] <- min(chosen); oe[j] <- max(chosen)
      }
    }
    for (drop in c(FALSE, TRUE)) {
      keep <- if (drop) os > 0 else rep(TRUE, n)
      res <- ORFik:::pmapFromTranscriptsCPP(first, last, s, en, rep(1L, n), direction, drop)
      expect_identical(res$ranges, list(os[keep], oe[keep]))
    }
  }
})

test_that("adaptive group accessors preserve compressed metadata on both paths", {
  for (exons in c(17L, 33L, 200L)) {
    one <- GRangesList(tx = GRanges("1", IRanges(seq.int(100L, by = 30L,
      length.out = exons), width = 20L), "+"))
    tx <- rep(one, if (exons == 17L) 1000L else if (exons == 33L) 2048L else 5000L)
    positions <- tx@partitioning@end
    expect_identical(strandPerGroup(tx, FALSE), as.character(tx@unlistData@strand)[positions])
    expect_identical(seqnamesPerGroup(tx, FALSE), as.character(tx@unlistData@seqnames)[positions])
    expect_identical(strandPerGroup(tx, TRUE), tx@unlistData@strand[positions])
    expect_identical(seqnamesPerGroup(tx, TRUE), tx@unlistData@seqnames[positions])
  }
})

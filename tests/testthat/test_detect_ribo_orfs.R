context("detect_ribo_orfs")
library(ORFik)

test_that("covRle statistics path equals the data.table path", {
  tx <- GRangesList(a = GRanges("1", IRanges(c(10, 50), c(20, 60)), "+"),
                    b = GRanges("1", IRanges(c(200, 100), c(230, 120)), "-"))
  seqlevels(tx) <- "1"; seqlengths(tx) <- 300
  set.seed(1)
  fw <- RleList("1" = Rle(rpois(300, 2))); rv <- RleList("1" = Rle(rpois(300, 3)))
  seqinfo(fw) <- seqinfo(tx); seqinfo(rv) <- seqinfo(tx)
  cov <- covRle(fw, rv)
  # ORFs/windows incl. a split ORF, a width-2 window (no frame score) and an empty window
  grl <- GRangesList(GRanges("1", IRanges(c(12, 50), c(20, 55)), "+"),
                     GRanges("1", IRanges(c(205, 100), c(230, 110)), "-"),
                     GRanges("1", IRanges(15, 16), "+"),
                     GRanges("1", IRanges(10, 9), "+"))
  fast <- ORFik:::coveragePerORFStatistics(grl, cov)
  # reference: data.table path on the same coverage (an RleList pair passed as covRle without strand mode)
  ref_cov <- coveragePerTiling(grl, cov, is.sorted = TRUE, as.data.table = TRUE, withFrames = TRUE)
  ref <- ref_cov[, .(sum = sum(count), median = median(count), mean = mean(count)), by = genes]
  expect_equal(fast$sum[1:3], ref$sum)
  expect_equal(fast$median[1:3], ref$median)
  expect_equal(fast$mean[1:3], ref$mean)
  sc <- orfScore(grl[1:3], cov, TRUE, coverage = ref_cov, stop3 = FALSE)
  expect_equal(fast$F1[1:3], sc$frame_zero_RP); expect_equal(fast$ORFScores[1:3], sc$ORFScores)
  expect_equal(unlist(fast[4, -1]), c(sum = 0, median = 0, mean = 0, F1 = 0, F2 = 0, F3 = 0, ORFScores = 0))
  expect_equal(fast$F1[3], 0)   # width <= 2: no frame scoring, as orfScore()
})

test_that("detect_ribo_orfs upstream_filter", {
  df <- ORFik.template.experiment()
  df <- df[df$libtype == "RFP", ][1, ]
  d <- tempfile()
  detect_ribo_orfs(df, file.path(d, "all"), c("uORF", "annotated"))
  detect_ribo_orfs(df, file.path(d, "na"), c("uORF", "annotated"), upstream_filter = "non_annotated")
  a <- riboORFs(df, type = "table", file.path(d, "all")); b <- riboORFs(df, type = "table", file.path(d, "na"))
  expect_equal(nrow(a), nrow(b))
  expect_true(all(b$predicted >= a$predicted))                 # only adds predictions
  expect_true(all(a$predicted[b$type != "annotated"] == b$predicted[b$type != "annotated"]))
  expect_error(detect_ribo_orfs(df, d, "uORF", upstream_filter = "x"))
})

test_that("detect_ribo_orfs with several libraries equals one call per library", {
  df <- ORFik.template.experiment()
  df <- df[df$libtype == "RFP", ][1:2, ]
  d <- tempfile()
  detect_ribo_orfs(df, file.path(d, "both"), c("uORF", "annotated"))
  detect_ribo_orfs(df[1, ], file.path(d, "one"), c("uORF", "annotated"))
  detect_ribo_orfs(df[2, ], file.path(d, "two"), c("uORF", "annotated"))
  tab <- function(dir) lapply(sort(list.files(dir, "_prediction_table.rds$", full.names = TRUE)), readRDS)
  both <- tab(file.path(d, "both"))
  expect_length(both, 2)
  expect_equal(both[[1]], tab(file.path(d, "one"))[[1]])
  expect_equal(both[[2]], tab(file.path(d, "two"))[[1]])
})

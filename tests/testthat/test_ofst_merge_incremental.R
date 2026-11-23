incremental_fixture <- function(seed, chrs = c("1", "2", "10", "X", "MT"), n = 400) {
  set.seed(seed)
  d <- data.table::data.table(seqnames = factor(sample(chrs, n, TRUE), levels = chrs),
                              start = sample.int(300L, n, TRUE),
                              strand = factor(sample(c("+", "-", "*"), n, TRUE), levels = c("+", "-", "*")),
                              size = sample(26:30, n, TRUE), score = sample.int(5L, n, TRUE))
  d <- d[, .(score = sum(score)), keyby = .(seqnames, start, strand, size)]
  data.table::setkey(d, NULL)
  d
}

incremental_write <- function(tables) {
  paths <- vapply(seq_along(tables), function(i) tempfile(fileext = ".ofst"), "")
  for (i in seq_along(tables)) fst::write_fst(tables[[i]], paths[i])
  paths
}

test_that("adding inputs to an existing merged base equals a full merge", {
  tabs <- lapply(1:8, incremental_fixture)
  tabs[[8]] <- incremental_fixture(8, chrs = c("2", "Y", "KI270"))   # chromosomes missing from the base
  paths <- incremental_write(tabs)
  full <- ofst_merge(paths, keep_all_scores = FALSE)
  base <- tempfile(fileext = ".ofst")
  fst::write_fst(ofst_merge(paths[1:5], keep_all_scores = FALSE), base)
  got <- ofst_merge(paths[6:8], keep_all_scores = FALSE, base = base)
  expect_identical(got[], full[])
  # small units and several workers: same result
  keys <- c("seqnames", "start", "strand", "size")
  for (unit_rows in c(1, 300)) for (workers in 1:2) {
    d <- .ofst_finalize(.ofst_merge_sorted(paths[6:8], keys, unit_rows = unit_rows, workers = workers,
                                           verbose = FALSE, base = base, subtract = rep(FALSE, 3)))
    expect_identical(d[], full[])
  }
})

test_that("subtracting inputs contained in the base restores the smaller merge", {
  paths <- incremental_write(lapply(1:6, incremental_fixture))
  full <- tempfile(fileext = ".ofst")
  fst::write_fst(ofst_merge(paths, keep_all_scores = FALSE), full)
  got <- ofst_merge(paths[1], keep_all_scores = FALSE, base = full, subtract = paths[c(1, 5, 6)])
  # base (1..6) + 1 - 1 - 5 - 6 = 1..4: rows only in inputs 5 and 6 drop out
  expect_identical(got[], ofst_merge(paths[1:4], keep_all_scores = FALSE)[])
  # replace: remove one input, add a changed version of it
  changed <- incremental_fixture(99)
  new_path <- incremental_write(list(changed))
  got <- ofst_merge(new_path, keep_all_scores = FALSE, base = full, subtract = paths[6])
  expect_identical(got[], ofst_merge(c(paths[1:5], new_path), keep_all_scores = FALSE)[])
  # subtracting something that is not in the base is an error, not a negative score
  expect_error(ofst_merge(paths[1], keep_all_scores = FALSE, base = paths[2], subtract = paths[3]),
               "negative")
  # a base that is not in merge order (strand "+" < "-" < "*" here) is grouped, not trusted
  raw <- incremental_fixture(7)
  expect_identical(ofst_merge(paths[1:2], keep_all_scores = FALSE, base = incremental_write(list(raw)))[],
                   ofst_merge(c(paths[1:2], incremental_write(list(raw))), keep_all_scores = FALSE)[])
})

test_that("mergeLibs adds libraries to an existing all.ofst in place", {
  paths <- incremental_write(lapply(1:5, incremental_fixture))
  df <- data.frame(libtype = "RFP", rep = 1:5)
  out_dir <- tempfile("merge_inc"); dir.create(out_dir)
  # a plain data.frame is enough: lib_names_full and paths are given explicitly
  mergeLibs(df[1:3, ], out_dir, mode = "all", keep_all_scores = FALSE, paths = paths[1:3],
            lib_names_full = paste0("l", 1:3))
  out <- file.path(out_dir, "all.ofst")
  mergeLibs(df[4:5, ], out_dir, mode = "all", keep_all_scores = FALSE, paths = paths[4:5],
            lib_names_full = paste0("l", 4:5), base = out)
  expect_identical(fst::read_fst(out, as.data.table = TRUE)[],
                   fst::read_fst({ fst::write_fst(ofst_merge(paths, keep_all_scores = FALSE), f <- tempfile()); f },
                                 as.data.table = TRUE)[])
  expect_error(mergeLibs(df[4:5, ], out_dir, mode = "all", keep_all_scores = TRUE, paths = paths[4:5],
                         lib_names_full = paste0("l", 4:5), base = out), "keep_all_scores")
})

test_that("incremental merge refuses unsorted inputs and non-pooled modes", {
  paths <- incremental_write(lapply(1:3, incremental_fixture))
  expect_error(ofst_merge(paths[2:3], keep_all_scores = TRUE, base = paths[1]), "keep_all_scores = FALSE")
  # chromosome split over two blocks in the base
  b <- incremental_fixture(5)
  b <- rbind(b[seqnames == "1"][1:5], b[seqnames == "2"], b[seqnames == "1"][-(1:5)])
  bad <- incremental_write(list(b))
  expect_error(ofst_merge(paths[2:3], keep_all_scores = FALSE, base = bad), "contiguous")
})

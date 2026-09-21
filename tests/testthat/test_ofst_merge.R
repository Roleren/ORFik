context("OFST merging and memory-safe normalization")
library(ORFik)

ofst_test_fixture <- function() {
  data.table::data.table(
    seqnames = factor(c("chr2", "chr1", "chr1", "chr1", "chr2", "chr10"),
                      levels = c("unused", "chr2", "chr10", "chr1")),
    start = c(10L, 20L, 20L, 20L, 10L, 30L),
    strand = factor(c("+", "+", "+", "-", "+", "*"), levels = c("-", "*", "+")),
    cigar = c("10M", "10M", "5M100N5M", "10M", "10M", "2S8M"),
    score = c(1L, 2L, 3L, 4L, 5L, 6L))
}

ofst_test_plain <- function(x) {
  x <- data.table::copy(x)
  for (col in names(x)) if (is.factor(x[[col]]))
    data.table::set(x, j = col, value = as.character(x[[col]]))
  x
}

ofst_test_canonical <- function(x) {
  x <- ofst_test_plain(x)
  data.table::setcolorder(x, sort(names(x)))
  data.table::setorderv(x, names(x))
  as.data.frame(x)
}

ofst_test_reference <- function(tables, keepCigar = TRUE) {
  x <- data.table::rbindlist(lapply(tables, ofst_test_plain), use.names = TRUE)
  keys <- setdiff(names(x), c("score", if (!keepCigar) "cigar"))
  x[, .(score = sum(score, na.rm = TRUE)), by = keys]
}

ofst_test_merge <- function(tables, ...) {
  paths <- vapply(seq_along(tables), function(i) tempfile("orfik-merge-", fileext = ".ofst"), "")
  on.exit(unlink(paths), add = TRUE)
  for (i in seq_along(tables)) fst::write_fst(tables[[i]], paths[i])
  suppressMessages(ORFik:::ofst_merge(paths, lib_names = paste0("lib", seq_along(paths)), ...))
}

ofst_test_expect_total <- function(result, tables, keepCigar = TRUE) {
  result <- data.table::copy(result)
  scores <- grep("^lib[0-9]+$", names(result), value = TRUE)
  if (length(scores)) result[, (scores) := NULL]
  expect_equal(ofst_test_canonical(result),
               ofst_test_canonical(ofst_test_reference(tables, keepCigar)))
}

test_that("normalization retains compact factors and accepts character columns", {
  x <- ofst_test_fixture()
  y <- ORFik:::.normalize_dt(data.table::copy(x))
  expect_identical(y$seqnames, x$seqnames)
  expect_identical(y$strand, x$strand)
  expect_identical(y$cigar, x$cigar)
  z <- ORFik:::.normalize_dt(ofst_test_plain(x))
  expect_true(is.factor(z$seqnames))
  expect_true(is.factor(z$strand))
  expect_equal(ofst_test_canonical(z), ofst_test_canonical(x))
  x[, cigar := factor(cigar, levels = rev(unique(cigar)))]
  expect_identical(ORFik:::.normalize_dt(x)$cigar, z$cigar)
})

test_that("ordered, empty, missing and malformed factors are handled by labels", {
  variants <- list(
    ordered(c("chr2", "chr1", NA), levels = c("chr2", "chr1", "unused")),
    factor(character()), factor(c(NA_character_, NA_character_)),
    structure(c(1L, 0L, -1L, 5L, NA_integer_, 2L), levels = c("chr1", "chr2"), class = "factor"),
    structure(c(1L, 2L, 3L), levels = c("chr1", NA_character_, "chr2"), class = "factor"),
    structure(c(1L, 2L, 3L), levels = c("chr1", "chr1", "chr2"), class = "factor"))
  expected <- list(c("chr2", "chr1", NA), character(), c(NA, NA),
                   c("chr1", NA, NA, NA, NA, "chr2"), c("chr1", NA, "chr2"), c("chr1", "chr1", "chr2"))
  for (i in seq_along(variants)) {
    z <- ORFik:::.ofst_factor(variants[[i]], "seqnames")
    expect_true(is.factor(z))
    expect_false(is.ordered(z))
    expect_equal(as.character(z), as.character(expected[[i]]))
  }
  x <- ofst_test_fixture()[1:3]
  x[, cigar := structure(c(1L, 0L, 8L), levels = "10M", class = "factor")]
  expect_equal(ORFik:::.normalize_dt(x)$cigar, c("10M", NA, NA))
})

test_that("mixed and disjoint factor levels preserve all CIGAR alignment keys", {
  a <- ofst_test_fixture()
  b <- ofst_test_plain(a)
  b[, seqnames := c("chr3", "chr1", "chr1", "chr1", "chr3", "chr10")]
  b[, score := 2L * score]
  c <- data.table::copy(b)
  c[, seqnames := ordered(seqnames, levels = c("chr10", "chr3", "chr1", "unused"))]
  c[, strand := factor(strand, levels = c("*", "+", "-"))]
  c[, cigar := factor(cigar, levels = rev(unique(cigar)))]
  for (inputs in list(list(a, b), list(a, c), list(b, c), list(a, b, c))) {
    for (keep in c(FALSE, TRUE)) {
      actual <- ofst_test_merge(inputs, keep_all_scores = keep)
      ofst_test_expect_total(actual, inputs)
      expect_identical(levels(actual$strand), c("+", "-", "*"))
      expect_false("unused" %in% levels(actual$seqnames))
      expect_type(actual$cigar, "character")
      if (keep) {
        for (i in seq_along(inputs))
          expect_equal(sum(actual[[paste0("lib", i)]], na.rm = TRUE), sum(inputs[[i]]$score))
      }
    }
  }
})

test_that("sort TRUE preserves character-label ordering and sort FALSE preserves counts", {
  a <- ofst_test_fixture()
  for (keep in c(FALSE, TRUE)) {
    sorted <- ofst_test_merge(list(a, a), keep_all_scores = keep)
    expected <- ofst_test_plain(sorted)
    data.table::setorderv(expected, setdiff(names(expected), "score"))
    expect_equal(as.data.frame(ofst_test_plain(sorted)), as.data.frame(expected))
    unsorted <- ofst_test_merge(list(a, a), keep_all_scores = keep, sort = FALSE)
    expect_equal(ofst_test_canonical(unsorted), ofst_test_canonical(sorted))
  }
})

test_that("singletons collapse duplicate reads in both score modes", {
  a <- ofst_test_fixture()
  for (keep in c(FALSE, TRUE)) {
    result <- ofst_test_merge(list(a), keep_all_scores = keep)
    expect_s3_class(result, "data.table")
    ofst_test_expect_total(result, list(a))
    expect_equal(nrow(result), 5L)
  }
})

test_that("singleton and mixed-size split chunks preserve keys and scores", {
  a <- ofst_test_fixture()
  scenarios <- list(list(inputs = list(a, a), limit = 7),
                    list(inputs = list(a, a[1:2], a[3:4]), limit = 7),
                    list(inputs = list(a, a, a, a), limit = 13))
  for (case in scenarios) for (keep in c(FALSE, TRUE)) {
    result <- ofst_test_merge(case$inputs, keep_all_scores = keep, dt_max_index_size = case$limit)
    ofst_test_expect_total(result, case$inputs)
    unsplit <- ofst_test_merge(case$inputs, keep_all_scores = keep)
    expect_equal(ofst_test_canonical(result), ofst_test_canonical(unsplit))
  }
})

test_that("empty files work alone, together and beside nonempty inputs", {
  a <- ofst_test_fixture()
  empty <- a[0]
  for (inputs in list(list(empty), list(empty, empty), list(empty, a), list(a, empty))) {
    for (keep in c(FALSE, TRUE)) {
      result <- ofst_test_merge(inputs, keep_all_scores = keep)
      ofst_test_expect_total(result, inputs)
      expect_true(is.factor(result$seqnames))
      expect_identical(levels(result$strand), c("+", "-", "*"))
    }
  }
})

test_that("missing keys match consistently and missing scores contribute zero", {
  a <- ofst_test_fixture()
  a[1, seqnames := NA]
  a[2, strand := NA]
  a[3, cigar := NA_character_]
  a[4, start := NA_integer_]
  a[5:6, score := NA_integer_]
  for (keep in c(FALSE, TRUE)) {
    result <- ofst_test_merge(list(a, a), keep_all_scores = keep)
    ofst_test_expect_total(result, list(a, a))
    expect_equal(sum(result$score), 20)
  }
})

test_that("points, widths, ends and optional metadata remain merge keys", {
  a <- ofst_test_fixture()[1:2]
  a[, `:=`(seqnames = factor("chr1"), start = 100L, strand = factor("+"), cigar = NULL)]
  for (geometry in c("width", "end", "size")) {
    b <- data.table::copy(a)
    b[, (geometry) := c(10L, 20L)]
    b[, annotation := c("first", "second")]
    for (keep in c(FALSE, TRUE)) {
      result <- ofst_test_merge(list(b, b), keep_all_scores = keep)
      ofst_test_expect_total(result, list(b, b))
      expect_equal(nrow(result), 2L)
    }
  }
  a[, width := 1L]
  ofst_test_expect_total(ofst_test_merge(list(a), keep_all_scores = FALSE), list(a))
})

test_that("explicit CIGAR removal uses the new key set in both score modes", {
  a <- ofst_test_fixture()
  for (inputs in list(list(a), list(a, a))) for (keep in c(FALSE, TRUE)) {
    result <- ofst_test_merge(inputs, keep_all_scores = keep, keepCigar = FALSE)
    expect_false("cigar" %in% names(result))
    ofst_test_expect_total(result, inputs, keepCigar = FALSE)
  }
  for (keep in c(FALSE, TRUE)) {
    result <- ofst_test_merge(list(a, a), keep_all_scores = keep, keepCigar = FALSE,
                              dt_max_index_size = 7)
    ofst_test_expect_total(result, list(a, a), keepCigar = FALSE)
  }
})

test_that("paired CIGAR geometry and factor labels survive the generic reducer", {
  a <- ofst_test_fixture()
  data.table::setnames(a, c("start", "cigar"), c("start1", "cigar1"))
  a[, `:=`(start2 = start1 + 100L, cigar2 = factor(cigar1))]
  for (keep in c(FALSE, TRUE)) {
    result <- ofst_test_merge(list(a, a), keep_all_scores = keep)
    ofst_test_expect_total(result, list(a, a))
    expect_true(all(c("cigar1", "cigar2", "start1", "start2") %in% names(result)))
  }
})

test_that("scores preserve fractions and sums above the integer limit", {
  a <- ofst_test_fixture()[1]
  for (score in c(.Machine$integer.max, 0.25, 0, 2^33)) {
    data.table::set(a, j = "score", value = score)
    for (keep in c(FALSE, TRUE)) {
      result <- ofst_test_merge(list(a, a), keep_all_scores = keep)
      expect_equal(result$score, score * 2)
      if (score * 2 > .Machine$integer.max || score == .25) expect_type(result$score, "double")
    }
  }
  # Keep the source column genuinely integer here: a mixed numeric vector in
  # the loop above would otherwise only exercise double-input summation.
  a <- ofst_test_fixture()
  a[, score := rep.int(.Machine$integer.max, .N)]
  expect_type(a$score, "integer")
  for (keep in c(FALSE, TRUE)) {
    result <- suppressWarnings(ofst_test_merge(list(a, a), keep_all_scores = keep))
    suppressWarnings(ofst_test_expect_total(result, list(a, a)))
    expect_equal(sum(result$score), 12 * as.numeric(.Machine$integer.max))
    expect_type(result$score, "double")
  }
})

test_that("whole-valued double coordinates are accepted without precision loss", {
  a <- ofst_test_fixture()
  a[, start := c(1, 2, .Machine$integer.max, NA_real_, 0, 100)]
  y <- ORFik:::.normalize_dt(a)
  expect_type(y$start, "integer")
  expect_equal(y$start, c(1L, 2L, .Machine$integer.max, NA_integer_, 0L, 100L))
  for (value in c(0.5, Inf, -Inf, 2^31, -2^31)) {
    bad <- ofst_test_fixture()
    bad[, start := as.numeric(start)]
    bad[1, start := value]
    expect_error(ORFik:::.normalize_dt(bad), "32-bit integer-valued")
  }
})

test_that("unsupported column types and incomplete schemas fail clearly", {
  expect_error(ORFik:::.normalize_dt(as.data.frame(ofst_test_fixture())), "data.table")
  for (col in c("seqnames", "strand", "cigar")) {
    for (bad_value in list(1:6, rep(TRUE, 6), as.list(1:6))) {
      bad <- ofst_test_fixture()
      data.table::set(bad, j = col, value = bad_value)
      expect_error(ORFik:::.normalize_dt(bad), "must be character or factor")
    }
  }
  for (col in c("start", "score")) {
    bad <- ofst_test_fixture()
    data.table::set(bad, j = col, value = rep("1", 6))
    expect_error(ORFik:::.normalize_dt(bad), "must be integer or numeric")
  }
  for (col in c("seqnames", "start", "strand", "score")) {
    bad <- ofst_test_fixture()
    bad[, (col) := NULL]
    expect_error(ORFik:::.normalize_dt(bad), "Missing required OFST columns")
  }
  bad <- ofst_test_fixture()
  bad[, strand := "?"]
  expect_error(ORFik:::.normalize_dt(bad), "strand values")
  bad <- ofst_test_fixture()
  data.table::setnames(bad, "cigar", "start")
  expect_error(ORFik:::.normalize_dt(bad), "column names must be unique")
})

test_that("public validation rejects bad paths, flags, names and split limits", {
  for (paths in list(character(), NA_character_, "", 1, list("x")))
    expect_error(ORFik:::ofst_merge(paths), "file_paths")
  p <- tempfile("orfik-validate-", fileext = ".ofst")
  on.exit(unlink(p))
  fst::write_fst(ofst_test_fixture(), p)
  for (arg in c("keep_all_scores", "keepCigar", "sort")) {
    for (value in list(NA, logical(), c(TRUE, FALSE), "TRUE", 1)) {
      args <- list(file_paths = p)
      args[[arg]] <- value
      expect_error(do.call(ORFik:::ofst_merge, args), "must be TRUE or FALSE")
    }
  }
  for (names in list(NA_character_, "", 1))
    expect_error(ORFik:::ofst_merge(p, lib_names = names), "lib_names")
  for (name in c("score", "seqnames", "cigar"))
    expect_error(ORFik:::ofst_merge(p, lib_names = name), "must not collide")
  expect_error(ORFik:::ofst_merge(c(p, p), c("same", "same")), "must be unique")
  # Duplicate labels are harmless when no per-library columns are requested.
  result <- suppressMessages(ORFik:::ofst_merge(c(p, p), c("same", "same"), keep_all_scores = FALSE))
  expect_equal(sum(result$score), 42)
  for (limit in list(0, -1, NA_real_, Inf, 2^32, "5", c(1, 2)))
    expect_error(ORFik:::ofst_merge(p, dt_max_index_size = limit), "dt_max_index_size")
  for (splits in list(0, -1, .5, NA_real_, Inf, "2", c(1, 2)))
    expect_error(ORFik:::ofst_merge(p, max_splits = splits), "max_splits")
  expect_error(ORFik:::ofst_merge(p, dt_max_index_size = 2), "not enough for safe chunking")
  expect_error(ORFik:::ofst_merge(p, dt_max_index_size = 2, max_splits = .Machine$integer.max),
                "not enough for safe chunking")
  expect_error(ORFik:::ofst_merge_internal(list(), character()), "at least one.*table")
})

test_that("file schemas may reorder columns but may not omit required fields", {
  a <- ofst_test_fixture()
  b <- data.table::copy(a)
  data.table::setcolorder(b, rev(names(b)))
  for (keep in c(FALSE, TRUE)) ofst_test_expect_total(
    ofst_test_merge(list(a, b), keep_all_scores = keep), list(a, b))
  b[, cigar := NULL]
  expect_error(ofst_test_merge(list(a, b)), "identical columns")
  a[, score := NULL]
  expect_error(ofst_test_merge(list(a)), "Invalid OFST file.*Missing required")
})

test_that("the slow join works with factors, duplicate labels and chunk totals", {
  # Exercise the actual fallback body at a small threshold, without changing
  # the namespace or allocating billions of rows.
  sandbox <- new.env(parent = asNamespace("ORFik"))
  sandbox$.int32_max <- 7
  reducer <- ORFik:::.ofst_merge_lossless
  environment(reducer) <- sandbox
  a <- ORFik:::.normalize_dt(ofst_test_fixture())
  b <- ORFik:::.normalize_dt(ofst_test_plain(ofst_test_fixture()))
  # Pre-collapsed inputs, as produced by the first round.
  a <- a[, .(score = sum(score)), by = .(seqnames, start, strand, cigar)]
  b <- b[, .(score = sum(score)), by = .(seqnames, start, strand, cigar)]
  for (chunkified in c(FALSE, TRUE)) {
    labels <- if (chunkified) c("same", "same", "extra") else c("same", "same")
    actual <- suppressMessages(reducer(list(data.table::copy(a), data.table::copy(b)),
                                       labels, keep_all_scores = FALSE, chunkified = chunkified))
    ofst_test_expect_total(actual, list(a, b))
    expect_false(any(grepl("ofst_score", names(actual))))
  }
})

test_that("lossless reducer errors never invoke local lossy filtering", {
  sandbox <- new.env(parent = asNamespace("ORFik"))
  sandbox$.int32_max <- 7
  sandbox$.merge_by_keys_reduce <- function(...) stop("Total rows in the list is 8 which exceeds the limit")
  reducer <- ORFik:::.ofst_merge_lossless
  environment(reducer) <- sandbox
  a <- ofst_test_fixture()
  expect_error(suppressMessages(reducer(list(data.table::copy(a), data.table::copy(a)),
                                       c("a", "b"), keep_all_scores = FALSE)), "Total rows in the list")
  expect_identical(formals(ORFik:::ofst_merge_internal)$allow_filtering, TRUE)
  expect_identical(formals(ORFik:::ofst_merge)$allow_filtering, TRUE)
})

test_that("raw duplicate keys cannot multiply scores in the slow join", {
  sandbox <- new.env(parent = asNamespace("ORFik"))
  sandbox$.int32_max <- 7
  reducer <- ORFik:::.ofst_merge_lossless
  environment(reducer) <- sandbox
  a <- ORFik:::.normalize_dt(ofst_test_fixture())
  result <- suppressMessages(reducer(list(data.table::copy(a), data.table::copy(a)),
                                     c("same", "same"), keep_all_scores = FALSE))
  ofst_test_expect_total(result, list(a, a))
})

test_that("unused library labels never remove alignment keys from total-only chunks", {
  sandbox <- new.env(parent = asNamespace("ORFik"))
  sandbox$.int32_max <- 7
  reducer <- ORFik:::.ofst_merge_lossless
  environment(reducer) <- sandbox
  a <- ORFik:::.normalize_dt(ofst_test_fixture()[1:3])
  result <- suppressMessages(reducer(list(data.table::copy(a), data.table::copy(a)),
                                     c("start", "strand"), keep_all_scores = FALSE,
                                     chunkified = TRUE))
  ofst_test_expect_total(result, list(a, a))
})

test_that("consumed inputs and the old accumulator are released before grouping", {
  sandbox <- new.env(parent = asNamespace("ORFik"))
  sandbox$.int32_max <- 7
  remaining <- integer()
  sandbox$gc <- function(...) {
    caller <- parent.frame()
    if (exists("combined", caller, inherits = FALSE) && !is.null(caller$combined)) {
      remaining <<- c(remaining, length(caller$dt_list))
      expect_null(caller$dt)
    }
    base::gc(...)
  }
  reducer <- ORFik:::.ofst_merge_lossless
  environment(reducer) <- sandbox
  a <- ORFik:::.normalize_dt(ofst_test_fixture()[1:3])
  result <- suppressMessages(reducer(lapply(1:4, function(i) data.table::copy(a)),
                                      paste0("lib", 1:4), keep_all_scores = FALSE))
  expect_identical(remaining, c(2L, 1L, 0L))
  ofst_test_expect_total(result, rep(list(a), 4))
})

test_that("mergeLibs writes valid OFST output for both RNA and point fixtures", {
  df <- ORFik::ORFik.template.experiment()
  df <- df[df$libtype == "RNA", ][1:2, ]
  inputs <- list(ofst_test_fixture(), ofst_test_fixture())
  paths <- vapply(1:2, function(i) tempfile("orfik-merge-libs-", fileext = ".ofst"), "")
  out <- tempfile("orfik-merge-output-")
  dir.create(out)
  on.exit(unlink(c(paths, file.path(out, c("all.ofst", "all.ofst.removal_summary.rds")))), add = TRUE)
  on.exit(unlink(out, recursive = FALSE), add = TRUE)
  for (point in c(FALSE, TRUE)) {
    if (point) inputs <- lapply(inputs, function(x) { x[, `:=`(cigar = NULL, width = 1L)]; x })
    for (i in 1:2) fst::write_fst(inputs[[i]], paths[i])
    for (keep in c(FALSE, TRUE)) {
      suppressMessages(ORFik::mergeLibs(df, out_dir = out, mode = "all", paths = paths,
                                       keep_all_scores = keep, lib_names_full = c("lib1", "lib2")))
      result <- fst::read_fst(file.path(out, "all.ofst"), as.data.table = TRUE)
      ofst_test_expect_total(result, inputs)
      expect_s4_class(ORFik::import.ofst(file.path(out, "all.ofst")),
                      if (point) "GRanges" else "GAlignments")
    }
  }
})

test_that("OFST round trips preserve RNA CIGARs and scores", {
  a <- ofst_test_fixture()
  result <- ofst_test_merge(list(a, a), keep_all_scores = FALSE)
  p <- tempfile("orfik-roundtrip-", fileext = ".ofst")
  on.exit(unlink(p))
  ORFik::export.ofst(result, p)
  reads <- ORFik::import.ofst(p)
  expect_s4_class(reads, "GAlignments")
  expect_equal(GenomicAlignments::cigar(reads), result$cigar)
  expect_equal(as.character(GenomicRanges::seqnames(reads)), as.character(result$seqnames))
  expect_equal(as.character(GenomicRanges::strand(reads)), as.character(result$strand))
  expect_equal(as.numeric(S4Vectors::mcols(reads)$score), as.numeric(result$score))
})

test_that("randomized factor encodings and chunking conserve every alignment key", {
  set.seed(947)
  for (iteration in seq_len(4)) {
    tables <- lapply(seq_len(4), function(i) {
      x <- data.table::data.table(seqnames = sample(c("chr1", "chr2", "chr10", NA), 20, TRUE),
                                 start = sample(c(1:4, NA_integer_), 20, TRUE),
                                 strand = sample(c("+", "-", "*", NA), 20, TRUE),
                                 cigar = sample(c("10M", "5M100N5M", "2S8M", NA), 20, TRUE),
                                 score = sample(c(0:5, NA_integer_), 20, TRUE))
      for (col in c("seqnames", "strand", "cigar")) if (runif(1) < .7) {
        lv <- sample(unique(stats::na.omit(x[[col]])))
        data.table::set(x, j = col, value = factor(x[[col]], levels = lv))
      }
      x
    })
    for (keep in c(FALSE, TRUE)) for (keep_cigar in c(FALSE, TRUE)) {
      result <- ofst_test_merge(tables, keep_all_scores = keep, keepCigar = keep_cigar,
                                dt_max_index_size = 41)
      ofst_test_expect_total(result, tables, keep_cigar)
    }
  }
})

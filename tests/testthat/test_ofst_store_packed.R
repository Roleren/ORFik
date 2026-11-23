# Packed-key kernels of the OFST store (src/ofst_store.cpp): k-way union, match, merge order
# and sum merge must equal a plain column-wise reference for every packing width, including
# the extreme int32 values.

int_max <- .Machine$integer.max
int_min <- -.Machine$integer.max   # NA_integer_ is -2^31, so this is the smallest valid key

kt <- function(...) data.table::data.table(...)

# random sorted unique key table (seq, start, strand, size) drawn from given value pools
rand_keys <- function(n, start_pool, seq_pool = 1:3, strand_pool = 1:3, size_pool = 20:35, seed = 1) {
  set.seed(seed)
  d <- unique(kt(seq = sample(seq_pool, n, TRUE), start = sample(start_pool, n, TRUE),
                 strand = sample(strand_pool, n, TRUE), size = sample(size_pool, n, TRUE)))
  data.table::setorderv(d, names(d))
  d
}

# reference for ofst_store_kway: ids of stored keys from `ord`, new keys n_old + rank in key order
kway_reference <- function(stored_sorted, ord, inputs) {
  kc <- names(stored_sorted)
  all <- unique(data.table::rbindlist(c(list(stored_sorted), inputs)))
  data.table::setorderv(all, kc)
  st <- data.table::copy(stored_sorted)[, id := ord]
  all <- st[all, on = kc]
  new <- is.na(all$id)
  all$id[new] <- nrow(stored_sorted) + seq_len(sum(new))
  list(order = all$id,
       new_keys = all[new, kc, with = FALSE],
       rows = lapply(inputs, function(x) all[x, on = kc]$id))
}

check_kway <- function(stored_sorted, ord, inputs, threads) {
  kc <- names(stored_sorted)
  # the kernel takes the stored keys in row-id order: row ord[i] holds the i-th smallest key
  by_id <- stored_sorted[order(ord)]
  got <- ofst_store_kway(as.list(by_id), as.integer(ord), lapply(inputs, as.list), nrow(stored_sorted), threads)
  ref <- kway_reference(stored_sorted, ord, inputs)
  expect_identical(got$order, ref$order)
  expect_identical(got$rows, ref$rows)
  expect_identical(as.data.frame(stats::setNames(got$new_keys, kc)), as.data.frame(ref$new_keys))
  got
}

test_that("keys pack into the expected number of 64-bit words", {
  typical <- list(list(1L, 5020345L, 2L, 21L), list(1L, 46699833L, 3L, 33L))
  expect_identical(ofst_packed_key_words(typical), 1L)
  # full int32 range in start: 32 bits; seq + strand + size still fit in word 1
  wide_start <- list(list(1L, int_min, 1L, 20L), list(3L, int_max, 3L, 35L))
  expect_identical(ofst_packed_key_words(wide_start), 1L)
  # two full-range columns fill one word, the rest spill into word 2
  two_full <- list(list(1L, int_min, 1L, int_min), list(3L, int_max, 3L, int_max))
  expect_identical(ofst_packed_key_words(two_full), 2L)
  # six full-range columns: 3 words (a column never straddles words)
  six <- list(as.list(rep(int_min, 6)), as.list(rep(int_max, 6)))
  expect_identical(ofst_packed_key_words(six), 3L)
  # constant columns take no bits
  expect_identical(ofst_packed_key_words(list(list(7L, 7L), list(7L, 7L))), 1L)
})

test_that("k-way union equals the reference for every packing width and thread count", {
  cases <- list(
    one_word = list(start = 1:5000, size = 20:35),
    extreme_start = list(start = c(int_min, int_min + 1L, -1L, 0L, 1:3000, int_max - 1L, int_max), size = 20:35),
    two_words = list(start = c(int_min, 0L, 1:3000, int_max), size = c(int_min, 1L, 20:35, int_max)))
  for (nm in names(cases)) {
    cs <- cases[[nm]]
    stored <- rand_keys(4e4, cs$start, size_pool = cs$size, seed = 1)
    ord <- sample(nrow(stored))   # stored row ids need not follow key order
    inputs <- lapply(2:6, function(s) rand_keys(5e4, cs$start, size_pool = cs$size, seed = s))
    inputs[[6]] <- stored[0]                      # empty input
    inputs[[7]] <- stored[c(1, .N)]               # only stored keys (first and last)
    for (th in c(1L, 4L)) {
      got <- check_kway(stored, ord, inputs, th)
      if (nm == "two_words") expect_identical(got$words, 2L) else expect_identical(got$words, 1L)
    }
    # empty store: fresh build, ids in key order
    got <- check_kway(stored[0], integer(), inputs, 4L)
    expect_identical(got$order, seq_along(got$order))
  }
})

test_that("k-way kernel rejects unsorted, duplicated and missing keys", {
  s <- kt(seq = 1L, start = c(1L, 2L, 3L), strand = 2L, size = 28L)
  bad_order <- s[c(2, 1, 3)]
  dup <- s[c(1, 1, 2)]
  expect_error(ofst_store_kway(as.list(s[0]), integer(), list(as.list(bad_order)), 0, 1L), "not sorted and unique")
  expect_error(ofst_store_kway(as.list(s[0]), integer(), list(as.list(dup)), 0, 1L), "not sorted and unique")
  expect_error(ofst_store_kway(as.list(bad_order), 1:3, list(as.list(s)), 3, 1L), "stored keys are not sorted")
  expect_error(ofst_store_kway(as.list(s), c(1L, 1L, 2L), list(as.list(s)), 3, 1L), "stored keys are not sorted")
  expect_error(ofst_store_kway(as.list(s), c(1L, 2L, 4L), list(as.list(s)), 3, 1L), "not a permutation")
  na <- data.table::copy(s); na$start[2] <- NA_integer_
  expect_error(ofst_store_kway(as.list(s[0]), integer(), list(as.list(na)), 0, 1L), "missing value")
  expect_error(ofst_store_kway(as.list(s), 1:3, list(as.list(s)), 2, 1L), "n_old")
  expect_error(ofst_store_kway(as.list(s), 1:2, list(as.list(s)), 3, 1L), "order length")
  expect_error(ofst_store_kway(as.list(s), 1:3, list(as.list(s[, 1:3])), 3, 1L), "key column counts")
})

test_that("packed match, merge order and sum merge equal the column-wise reference at the int32 limits", {
  for (sz in list(20:35, c(int_min, 0L, int_max))) {
    t <- rand_keys(3000, c(int_min, -5L, 0L, 7L, int_max), size_pool = sz, seed = 11)
    q <- rand_keys(2000, c(int_min, -5L, 0L, 8L, int_max), size_pool = sz, seed = 12)
    kc <- c("seq", "start", "strand", "size")
    ref <- t[, i := .I][q, on = kc]$i
    ref[is.na(ref)] <- 0L
    expect_identical(ofst_sorted_match(as.list(t[, ..kc]), as.list(q)), ref)
    t[, i := NULL]
    newk <- q[!t, on = kc]
    ord_a <- sample(nrow(t))
    got <- ofst_sorted_merge_order(as.list(t), ord_a, as.list(newk), nrow(t) + seq_len(nrow(newk)))
    both <- data.table::rbindlist(list(t[, id := ord_a], newk[, id := nrow(t) + seq_len(.N)]))
    data.table::setorderv(both, kc)
    expect_identical(got, both$id)
    t[, id := NULL]; newk[, id := NULL]
    expect_error(ofst_sorted_merge_order(as.list(t), ord_a, as.list(t[1]), 1L), "already in the table")
    a_s <- as.numeric(seq_len(nrow(t))); b_s <- as.numeric(seq_len(nrow(q)))
    sm <- ofst_sorted_sum_merge(as.list(t), a_s, as.list(q), b_s, FALSE, FALSE)
    rs <- data.table::rbindlist(list(t[, score := a_s], q[, score := b_s]))[, .(score = sum(score)), keyby = kc]
    expect_identical(as.data.frame(data.table::as.data.table(sm$keys)[, score := sm$score]), as.data.frame(rs))
    t[, score := NULL]; q[, score := NULL]
  }
  expect_error(ofst_sorted_match(list(1:3), list(c(2L, 1L))), "not sorted and unique")
})

test_that("store build and add are exact with keys at the int32 limits", {
  lib <- function(seed) {
    set.seed(seed)
    n <- 400
    d <- data.table::data.table(seqnames = factor(sample(c("1", "2"), n, TRUE)),
                                start = sample(c(1L, 2L, 1e6L, int_max - 1L, int_max), n, TRUE),
                                strand = factor(sample(c("+", "-", "*"), n, TRUE), levels = c("+", "-", "*")),
                                size = sample(c(1L, 28L, 1000L), n, TRUE), score = sample.int(9L, n, TRUE))
    d[, .(score = sum(score)), by = .(seqnames, start, strand, size)][sample(.N)]
  }
  dir <- tempfile("packed_in"); dir.create(dir)
  paths <- file.path(dir, paste0("lib", 1:5, ".ofst"))
  for (i in 1:5) fst::write_fst(lib(i), paths[i])
  store <- file.path(tempfile("packed_store"), "s")
  ofst_store_build(store, paths[1:3], batch_size = 2, partition_min_rows = 1, verbose = FALSE)
  ofst_store_add(store, paths[4:5], verbose = FALSE)
  ref <- suppressWarnings(suppressMessages(ofst_merge(paths, keep_all_scores = FALSE)))
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], ref[])
  expect_true(ofst_store_verify(store))
  ofst_store_remove(store, "lib2", verbose = FALSE)
  ref <- suppressWarnings(suppressMessages(ofst_merge(paths[-2], keep_all_scores = FALSE)))
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], ref[])
})

test_that("k-way kernel reads sorted stored keys sequentially and returns the sorted union", {
  for (two in c(FALSE, TRUE)) {
    sp <- if (two) c(int_min, 1L, 20:35, int_max) else 20:35
    stored <- rand_keys(4e4, 1:5000, size_pool = sp, seed = 11)   # in key order
    ord <- sample(nrow(stored))                                     # their row ids
    by_id <- stored[order(ord)]                                     # the keys in row-id order
    inputs <- lapply(12:15, function(s) rand_keys(3e4, 1:6000, size_pool = sp, seed = s))
    for (th in c(1L, 4L)) {
      a <- ofst_store_kway(as.list(by_id), ord, lapply(inputs, as.list), nrow(stored), th, FALSE, TRUE)
      b <- ofst_store_kway(as.list(stored), ord, lapply(inputs, as.list), nrow(stored), th, TRUE, TRUE)
      expect_identical(a[c("rows", "order", "new_keys", "keys")], b[c("rows", "order", "new_keys", "keys")])
      ref <- kway_reference(stored, ord, inputs)
      expect_identical(b$order, ref$order)
      union <- unique(data.table::rbindlist(c(list(stored), inputs)))
      data.table::setorderv(union, names(union))
      expect_identical(lapply(b$keys, identity), unname(as.list(union)))
    }
  }
})

test_that("the sorted key file is written, linked on removals and recreated for older stores", {
  lib <- function(seed) {
    set.seed(seed)
    d <- data.table::data.table(seqnames = factor(sample(c("1", "2", "MT"), 3000, TRUE)), start = sample.int(5000L, 3000, TRUE),
                                strand = factor(sample(c("+", "-"), 3000, TRUE), levels = c("+", "-", "*")),
                                size = sample(26:30, 3000, TRUE), score = sample.int(5L, 3000, TRUE))
    d[, .(score = sum(score)), by = .(seqnames, start, strand, size)]
  }
  dir <- tempfile("skeys_in"); dir.create(dir)
  paths <- file.path(dir, paste0("lib", 1:5, ".ofst"))
  for (i in 1:5) fst::write_fst(lib(i), paths[i])
  store <- file.path(tempfile("ofst_store"), "store")
  ofst_store_build(store, paths[1:3], workers = 2, verbose = FALSE)
  skeys_of <- function() {
    m <- ofst_store_manifest(store)
    .ofst_store_gen_file(store, m$parts$part, "skeys", m$parts$gen)
  }
  expect_true(all(file.exists(skeys_of())))
  expect_true(ofst_store_verify(store))
  ofst_store_remove(store, "lib2", verbose = FALSE)        # no new keys: linked
  expect_true(all(file.exists(skeys_of())))
  # a store written before the sorted key file existed: the next add creates it
  unlink(skeys_of())
  ofst_store_add(store, paths[4], workers = 1, verbose = FALSE)
  expect_true(all(file.exists(skeys_of())))
  ofst_store_add(store, paths[5], workers = 2, verbose = FALSE)
  expect_true(ofst_store_verify(store))
  ref <- suppressMessages(ofst_merge(paths[-2], keep_all_scores = FALSE))
  expect_identical(ofst_store_export_all(store, verbose = FALSE)[], ref[])
  out <- tempfile(fileext = ".ofst")
  ofst_store_export_all(store, out, workers = 2, verbose = FALSE)
  expect_identical(fst::read_fst(out, as.data.table = TRUE)[], ref[])
  # a corrupted sorted key file is caught by verify
  f <- skeys_of()[1]
  k <- fst::read_fst(f, as.data.table = TRUE); k[1, start := start + 1L]; fst::write_fst(k, f)
  expect_error(ofst_store_verify(store), "sorted key file")
})

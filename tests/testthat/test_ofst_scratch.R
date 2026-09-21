context("OFST scratch fallback, ownership and abort cleanup")

cache_fixture <- function() data.table::data.table(
  seqnames = rep(c("chr1", "chr2"), 4), start = 1:8, strand = "+", cigar = "10M", score = 1L)

cache_setup <- function() {
  root <- tempfile("cache-test-")
  dir.create(root)
  primary <- file.path(root, "primary")
  output <- file.path(root, "output")
  dir.create(primary)
  dir.create(output)
  paths <- file.path(root, c("input1.ofst", "input2.ofst"))
  for (path in paths) fst::write_fst(cache_fixture(), path)
  list(root = root, primary = primary, output = output, paths = paths)
}

cache_merge <- function(x, ...) suppressMessages(ORFik:::ofst_merge(x$paths,
  c("lib1", "lib2"), keep_all_scores = FALSE, filter_target_rows = 4,
  filter_chunk_rows = 4, filter_tmpdir = x$primary, filter_fallback_dir = x$output, ...))

test_that("missing primary storage falls back and ordinary success cleans both caches", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  x$primary <- file.path(x$root, "missing")
  result <- cache_merge(x)
  expect_equal(nrow(result), 4L)
  expect_true(attr(result, "removal_summary")$scratch_fallback_used)
  expect_equal(attr(result, "removal_summary")$scratch_parent, normalizePath(x$output))
  expect_length(list.files(x$output, all.files = TRUE, no.. = TRUE), 0L)
})

test_that("late scratch-write failure restarts globally and preserves selection and inputs", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  expected <- cache_merge(x)
  checksum <- tools::md5sum(x$paths)
  original <- ORFik:::.ofst_spill
  primary_writes <- 0L
  failed_cache <- NULL
  testthat::local_mocked_bindings(.ofst_spill = function(dt, scratch, prefix = "part-") {
    if (identical(dirname(scratch), normalizePath(x$primary))) {
      primary_writes <<- primary_writes + 1L
      if (primary_writes == 5L) {
        failed_cache <<- scratch
        ORFik:::.ofst_storage_abort("simulated primary disk full after several writes")
      }
    }
    original(dt, scratch, prefix)
  }, .package = "ORFik")
  actual <- cache_merge(x)
  expect_equal(primary_writes, 5L)
  expect_equal(actual$start, expected$start)
  expect_equal(actual$score, expected$score)
  expect_true(attr(actual, "removal_summary")$scratch_fallback_used)
  expect_false(dir.exists(failed_cache))
  expect_length(list.files(x$primary), 0L)
  expect_length(list.files(x$output), 0L)
  expect_identical(tools::md5sum(x$paths), checksum)
})

test_that("both storage failures report both causes and remove private caches", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  fst::write_fst(cache_fixture(), file.path(x$output, "all.ofst"))
  checksum <- tools::md5sum(file.path(x$output, "all.ofst"))
  attempts <- 0L
  testthat::local_mocked_bindings(.ofst_spill = function(...) {
    attempts <<- attempts + 1L
    ORFik:::.ofst_storage_abort(paste0("storage failure ", attempts))
  }, .package = "ORFik")
  expect_error(cache_merge(x), "both scratch locations failed.*storage failure 1.*storage failure 2")
  expect_equal(attempts, 2L)
  expect_length(list.files(x$primary), 0L)
  expect_identical(list.files(x$output), "all.ofst")
  expect_identical(tools::md5sum(file.path(x$output, "all.ofst")), checksum)
})

test_that("validation errors and interrupts do not trigger another attempt", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  for (kind in c("error", "interrupt")) {
    attempts <- 0L
    fun <- function(scratch) {
      attempts <<- attempts + 1L
      ORFik:::.ofst_spill(cache_fixture(), scratch)
      if (kind == "error") stop("invalid data deliberately")
      stop(structure(list(message = "user interrupt", call = NULL), class = c("interrupt", "condition")))
    }
    if (kind == "error") expect_error(ORFik:::.ofst_with_scratch(x$primary, x$output, fun), "invalid data") else {
      interrupted <- tryCatch(ORFik:::.ofst_with_scratch(x$primary, x$output, fun), interrupt = function(e) TRUE)
      expect_true(interrupted)
    }
    expect_equal(attempts, 1L)
    expect_length(list.files(x$primary), 0L)
    expect_length(list.files(x$output), 0L)
  }
})

test_that("fallback can be disabled and identical paths are not retried", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  for (fallback in list(NULL, x$primary)) {
    attempts <- 0L
    expect_error(ORFik:::.ofst_with_scratch(x$primary, fallback, function(scratch) {
      attempts <<- attempts + 1L
      ORFik:::.ofst_storage_abort("full disk")
    }), "No distinct filter_fallback_dir")
    expect_equal(attempts, 1L)
  }
  expect_error(ORFik:::.ofst_with_scratch(x$primary, file.path(x$root, "missing"),
    function(scratch) ORFik:::.ofst_storage_abort("primary full")), "Fallback directory.*unavailable")
  for (bad in list(NA_character_, "", 4, c(x$output, x$primary), matrix(x$output)))
    expect_error(suppressMessages(ORFik:::ofst_merge(x$paths, filter_fallback_dir = bad)), "filter_fallback_dir")
})

test_that("memory allocation failures are not mistaken for disk failures", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  attempts <- 0L
  testthat::local_mocked_bindings(write_fst = function(...) {
    attempts <<- attempts + 1L
    stop("cannot allocate vector of size 1 GB")
  }, .package = "fst")
  expect_error(cache_merge(x), "exhausted memory.*lower filter_chunk_rows")
  expect_equal(attempts, 1L)
  expect_length(list.files(x$primary), 0L)
  expect_length(list.files(x$output), 0L)
})

test_that("interrupting the fallback attempt removes its cache without another retry", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  attempts <- 0L
  result <- tryCatch(suppressMessages(ORFik:::.ofst_with_scratch(x$primary, x$output, function(scratch) {
    attempts <<- attempts + 1L
    ORFik:::.ofst_spill(cache_fixture(), scratch)
    if (attempts == 1L) ORFik:::.ofst_storage_abort("primary full")
    stop(structure(list(message = "stop fallback", call = NULL), class = c("interrupt", "condition")))
  })), interrupt = function(e) TRUE)
  expect_true(result)
  expect_equal(attempts, 2L)
  expect_length(list.files(x$primary), 0L)
  expect_length(list.files(x$output), 0L)
})

test_that("an interrupt while preparing output leaves previous output and summary intact", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  path <- file.path(x$output, "all.ofst")
  fst::write_fst(cache_fixture(), path)
  saveRDS("previous summary", paste0(path, ".removal_summary.rds"))
  checksum <- tools::md5sum(c(path, paste0(path, ".removal_summary.rds")))
  testthat::local_mocked_bindings(write_fst = function(...) {
    stop(structure(list(message = "interrupt writer", call = NULL), class = c("interrupt", "condition")))
  }, .package = "fst")
  result <- tryCatch(ORFik:::.ofst_save_merge(cache_fixture(), path), interrupt = function(e) TRUE)
  expect_true(result)
  expect_identical(tools::md5sum(names(checksum)), checksum)
  expect_setequal(list.files(x$output), c("all.ofst", "all.ofst.removal_summary.rds"))
})

test_that("mergeLibs defaults to output-folder fallback and leaves no scratch payload", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  df <- ORFik::ORFik.template.experiment()[1:2, ]
  suppressMessages(ORFik::mergeLibs(df, out_dir = x$output, paths = x$paths,
    lib_names_full = c("a", "b"), filter_target_rows = 4,
    filter_tmpdir = file.path(x$root, "unavailable")))
  expect_setequal(list.files(x$output), c("all.ofst", "all.ofst.removal_summary.rds"))
  summary <- readRDS(file.path(x$output, "all.ofst.removal_summary.rds"))
  expect_true(summary$scratch_fallback_used)
  expect_equal(summary$rows_after, 4)
})

test_that("stale cleanup never claims ordinary folders, symlinks or active caches", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  ordinary <- file.path(x$output, "orfik-ofst-filter-user-folder")
  dir.create(ordinary)
  saveRDS("user data", file.path(ordinary, "keep.rds"))
  live <- ORFik:::.ofst_new_cache(x$output)
  foreign <- ORFik:::.ofst_new_cache(x$output)
  owner <- ORFik:::.ofst_cache_owner(foreign)
  owner$process$host <- "another-host"
  ORFik:::.ofst_write_cache_owner(foreign, owner)
  link <- file.path(x$output, "orfik-ofst-filter-symlink")
  linked <- file.symlink(ordinary, link)
  suppressMessages(ORFik:::.ofst_cleanup_stale_caches(x$output))
  expect_true(dir.exists(live))
  expect_true(dir.exists(foreign))
  expect_true(file.exists(file.path(ordinary, "keep.rds")))
  if (linked) expect_true(nzchar(Sys.readlink(link)))
  expect_true(ORFik:::.ofst_cleanup_cache(live, current_run = TRUE))
})

test_that("PID reuse and mutually linked dead caches are handled without trusting arbitrary paths", {
  testthat::skip_if(is.null(ORFik:::.ofst_process_owner()$start), "Linux process identity unavailable")
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  anchor <- ORFik:::.ofst_new_cache(x$output, "anchor")
  primary <- ORFik:::.ofst_new_cache(x$primary)
  ORFik:::.ofst_link_cache(anchor, primary)
  for (path in c(anchor, primary)) {
    owner <- ORFik:::.ofst_cache_owner(path)
    owner$process$start <- "0"
    ORFik:::.ofst_write_cache_owner(path, owner)
  }
  expect_false(ORFik:::.ofst_owner_alive(ORFik:::.ofst_cache_owner(anchor)$process))
  suppressMessages(ORFik:::.ofst_cleanup_stale_caches(x$output))
  expect_false(dir.exists(anchor))
  expect_false(dir.exists(primary))
  fake <- ORFik:::.ofst_new_cache(x$output, "anchor")
  owner <- ORFik:::.ofst_cache_owner(fake)
  owner$process$start <- "0"
  owner$children <- list(list(path = x$root, id = "not-an-owned-cache"))
  ORFik:::.ofst_write_cache_owner(fake, owner)
  expect_false(ORFik:::.ofst_cleanup_cache(fake))
  expect_true(all(file.exists(x$paths)))
})

test_that("in-memory merging uses the same fallback and cleanup", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  result <- suppressMessages(ORFik:::ofst_merge_internal(list(cache_fixture(), cache_fixture()),
    c("a", "b"), filter_target_rows = 4, filter_tmpdir = file.path(x$root, "missing"),
    filter_fallback_dir = x$output))
  expect_equal(nrow(result), 4L)
  expect_true(attr(result, "removal_summary")$scratch_fallback_used)
  expect_length(list.files(x$output), 0L)
})

test_that("interrupted output publication retains recoverable files, not stale summary claims", {
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  path <- file.path(x$output, "all.ofst")
  # A directory at the sidecar target reliably makes the final rename fail.
  dir.create(paste0(path, ".removal_summary.rds"))
  expect_error(suppressWarnings(suppressMessages(ORFik:::.ofst_save_merge(cache_fixture(), path))),
               "summary could not replace.*correct summary is preserved")
  cache <- list.dirs(x$output, recursive = FALSE, full.names = TRUE)
  cache <- cache[grepl("^orfik-ofst-output-", basename(cache))]
  expect_length(cache, 1L)
  expect_true(file.exists(file.path(cache, "removal_summary.rds")))
  owner <- ORFik:::.ofst_cache_owner(cache)
  owner$process$start <- "0"
  ORFik:::.ofst_write_cache_owner(cache, owner)
  suppressMessages(ORFik:::.ofst_cleanup_stale_caches(x$output))
  expect_true(dir.exists(cache))
  expect_true(file.exists(path))
})

test_that("malformed process identities never constitute evidence that an owner died", {
  process <- ORFik:::.ofst_process_owner()
  for (change in list(list(start = NA_character_), list(start = "garbled"),
                      list(pid = Inf), list(pid = NA_real_), list(pid = c(1L, 2L)),
                      list(boot = "not-a-boot-id"), list(pid_namespace = NA_character_))) {
    candidate <- process
    candidate[names(change)] <- change
    expect_true(is.na(ORFik:::.ofst_owner_alive(candidate)))
  }
})

cache_wait_ready <- function(process, path) {
  deadline <- Sys.time() + 45
  while (!file.exists(path) && process$is_alive() && Sys.time() < deadline) Sys.sleep(.05)
  if (!file.exists(path)) stop("child did not reach the cache checkpoint: ", process$read_all_error())
}

test_that("real SIGINT cleans caches and SIGKILL caches are recovered on the next run", {
  testthat::skip_if_not_installed("callr")
  testthat::skip_if(Sys.info()[["sysname"]] != "Linux")
  repo <- system.file(package = "ORFik")
  x <- cache_setup()
  on.exit(unlink(x$root, recursive = TRUE))
  for (signal in c(2L, 9L)) {
    checkpoint <- file.path(x$root, paste0("ready-", signal, ".rds"))
    child <- callr::r_bg(function(repo, primary, output, checkpoint) {
      devtools::load_all(repo, recompile = FALSE)
      ORFik:::.ofst_with_scratch(primary, output, function(scratch) {
        saveRDS(1:10, file.path(scratch, "payload.rds"))
        saveRDS(list(scratch = scratch), checkpoint)
        repeat Sys.sleep(.1)
      })
    }, args = list(repo, x$primary, x$output, checkpoint),
    stdout = file.path(x$root, paste0("child-", signal, ".log")),
    stderr = file.path(x$root, paste0("child-", signal, ".err")))
    on.exit(if (child$is_alive()) child$kill(), add = TRUE)
    cache_wait_ready(child, checkpoint)
    primary_cache <- readRDS(checkpoint)$scratch
    output_cache <- list.dirs(x$output, recursive = FALSE, full.names = TRUE)
    expect_length(output_cache, 1L)
    suppressMessages(ORFik:::.ofst_cleanup_stale_caches(x$output))
    expect_true(dir.exists(primary_cache))
    expect_true(dir.exists(output_cache))
    child$signal(signal)
    child$wait(timeout = 10000)
    expect_false(child$is_alive())
    if (signal == 9L) {
      expect_true(dir.exists(primary_cache))
      expect_true(dir.exists(output_cache))
      suppressMessages(ORFik:::.ofst_cleanup_stale_caches(x$output))
    }
    expect_false(dir.exists(primary_cache))
    expect_false(dir.exists(output_cache))
    expect_true(all(file.exists(x$paths)))
  }
})

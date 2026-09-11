context("Explicit OFST row budgets and optional memory sizing")

budget_inputs <- function(fun) {
  root <- tempfile("budget-test-")
  dir.create(root)
  on.exit(unlink(root, recursive = TRUE))
  dt <- data.table::data.table(seqnames = "chr1", start = 1:12,
                              strand = "+", cigar = "10M", score = 1L)
  path <- file.path(root, "input.ofst")
  fst::write_fst(dt, path)
  fun(path, root, dt)
}

test_that("manual 5e8 budget bypasses the estimator in all public entry points", {
  testthat::local_mocked_bindings(.ofst_batch_budget = function(...) stop("estimator must not run"),
                                 .package = "ORFik")
  budget_inputs(function(path, root, dt) {
    direct <- suppressMessages(ORFik:::ofst_merge(path, keep_all_scores = FALSE,
      filter_target_rows = 6, filter_chunk_rows = 5e8, filter_tmpdir = root))
    internal <- suppressMessages(ORFik:::ofst_merge_internal(list(dt), "a", keep_all_scores = FALSE,
      filter_target_rows = 6, filter_chunk_rows = 5e8, filter_tmpdir = root))
    df <- ORFik::ORFik.template.experiment()[1, ]
    out <- file.path(root, "output")
    suppressMessages(ORFik::mergeLibs(df, out_dir = out, paths = path, lib_names_full = "a",
      keep_all_scores = FALSE, filter_target_rows = 6, filter_chunk_rows = 5e8))
    summaries <- list(attr(direct, "removal_summary"), attr(internal, "removal_summary"),
                      readRDS(file.path(out, "all.ofst.removal_summary.rds")))
    for (s in summaries) {
      expect_identical(s$filter_auto_memory, FALSE)
      expect_equal(s$requested_chunk_rows, 5e8)
      expect_equal(s$filter_chunk_rows, 5e8)
      expect_equal(s$batches, 1L)
    }
  })
})

test_that("automatic sizing is forwarded only after opt-in and records both budgets", {
  calls <- 0L
  testthat::local_mocked_bindings(.ofst_batch_budget = function(paths, rows, requested) {
    calls <<- calls + 1L
    expect_equal(requested, 5e8)
    4
  }, .package = "ORFik")
  budget_inputs(function(path, root, dt) {
    direct <- suppressMessages(ORFik:::ofst_merge(path, keep_all_scores = FALSE,
      filter_target_rows = 6, filter_chunk_rows = 5e8, filter_tmpdir = root, filter_auto_memory = TRUE))
    internal <- suppressMessages(ORFik:::ofst_merge_internal(list(dt), "a", keep_all_scores = FALSE,
      filter_target_rows = 6, filter_chunk_rows = 5e8, filter_tmpdir = root, filter_auto_memory = TRUE))
    out <- file.path(root, "output")
    suppressMessages(ORFik::mergeLibs(ORFik::ORFik.template.experiment()[1, ], out_dir = out,
      paths = path, lib_names_full = "a", keep_all_scores = FALSE, filter_target_rows = 6,
      filter_chunk_rows = 5e8, filter_auto_memory = TRUE))
    summaries <- list(attr(direct, "removal_summary"), attr(internal, "removal_summary"),
                      readRDS(file.path(out, "all.ofst.removal_summary.rds")))
    for (s in summaries) {
      expect_true(s$filter_auto_memory)
      expect_equal(s$requested_chunk_rows, 5e8)
      expect_equal(s$filter_chunk_rows, 4)
      expect_equal(s$batches, 3L)
    }
    expect_equal(direct$start, internal$start)
    expect_equal(calls, 3L)
  })
})

test_that("cache-heavy cgroup reproduces the old collapse but cannot override manual budgets", {
  budget_inputs(function(path, root, dt) {
    bytes <- as.double(object.size(dt)) / nrow(dt)
    headroom <- ceiling(8 * bytes * 33820.5)
    sandbox <- new.env(parent = asNamespace("ORFik"))
    sandbox$file.exists <- function(paths) paths %in% c("/proc/meminfo",
      "/sys/fs/cgroup/memory.max", "/sys/fs/cgroup/memory.current")
    usage_calls <- 0L
    sandbox$get_system_usage <- function(drive) {
      expect_true(is.na(drive))
      usage_calls <<- usage_calls + 1L
      list(Memory_Available_Bytes = 250 * 1024^3)
    }
    sandbox$readLines <- function(path, ...) switch(path,
      "/sys/fs/cgroup/memory.max" = format(500 * 1024^3, scientific = FALSE),
      "/sys/fs/cgroup/memory.current" = format(500 * 1024^3 - headroom, scientific = FALSE))
    sandbox$.ofst_read_part <- function(...) dt
    estimator <- ORFik:::.ofst_batch_budget
    environment(estimator) <- sandbox
    messages <- character()
    value <- withCallingHandlers(estimator(path, nrow(dt), 5e8), message = function(m) {
      messages <<- c(messages, conditionMessage(m))
      invokeRestart("muffleMessage")
    })
    expect_equal(value, 33820)
    expect_true(any(grepl("raw cgroup headroom", messages)))
    expect_true(any(grepl("filesystem cache", messages)))
    expect_true(any(grepl("sampled bytes/row", messages)))
    testthat::local_mocked_bindings(.ofst_batch_budget = estimator, .package = "ORFik")
    result <- suppressMessages(ORFik:::ofst_merge(path, keep_all_scores = FALSE,
      filter_target_rows = 6, filter_chunk_rows = 5e8, filter_tmpdir = root))
    expect_equal(attr(result, "removal_summary")$filter_chunk_rows, 5e8)
    expect_equal(usage_calls, 1L) # Only the explicitly invoked estimator.
  })
})

test_that("invalid opt-in flags fail even when no filtering is necessary", {
  budget_inputs(function(path, root, dt) {
    for (bad in list(NA, 1, "TRUE", logical(), c(TRUE, FALSE), matrix(TRUE))) {
      expect_error(ORFik:::ofst_merge(path, filter_auto_memory = bad), "filter_auto_memory must be TRUE or FALSE")
      expect_error(ORFik:::ofst_merge_internal(list(dt), "a", filter_auto_memory = bad), "filter_auto_memory must be TRUE or FALSE")
      expect_error(ORFik::mergeLibs(ORFik::ORFik.template.experiment()[1, ], paths = path,
        out_dir = root, filter_auto_memory = bad), "filter_auto_memory must be TRUE or FALSE")
    }
  })
})

test_that("available-memory parsing retains byte precision and handles unavailable values", {
  sandbox <- new.env(parent = asNamespace("ORFik"))
  parser <- ORFik:::.system_available_memory_bytes
  environment(parser) <- sandbox
  for (text in list("MemAvailable:   262144001 kB", "MemAvailable: 0 kB")) {
    sandbox$readLines <- function(...) text
    expected <- if (grepl("262144001", text)) 262144001 * 1024 else 0
    expect_equal(parser(), expected)
  }
  for (text in list(character(), "MemTotal: 500 kB", "MemAvailable: broken kB",
                   rep("MemAvailable: 1 kB", 2))) {
    sandbox$readLines <- function(...) text
    expect_true(is.na(parser()))
  }
  sandbox$readLines <- function(...) stop("unavailable")
  expect_true(is.na(parser()))
})

test_that("get_system_usage exposes available memory without changing existing fields", {
  testthat::skip_if(Sys.info()[["sysname"]] != "Linux")
  testthat::local_mocked_bindings(.system_available_memory_bytes = function() 250 * 1024^3,
                                 .package = "ORFik")
  usage <- ORFik::get_system_usage(drive = NA_character_)
  expect_equal(usage$Memory_Available_Bytes, 250 * 1024^3)
  expect_equal(usage$Memory_Available_GB, 250)
  expect_identical(names(usage)[1:9], c("CPU_Usage_Percent", "Memory_Usage_GB", "Memory_Total_GB",
    "Memory_Usage_Percent", "Drive", "Drive_Total", "Drive_Used", "Drive_Free", "Drive_Usage_Percent"))
  expect_true(all(c("CPU_Usage_Percent", "Memory_Usage_GB", "Memory_Total_GB",
    "Memory_Usage_Percent", "Drive", "Drive_Free") %in% names(usage)))
})

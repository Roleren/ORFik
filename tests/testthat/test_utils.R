context("Utils")
library(ORFik)


test_that("NGS pair detection works as intended", {
  # As data.table
  paths <- c("/aligned_GRCz10/wigs/Bud_fwd.wig", "/aligned_GRCz10/wigs/Bud_MZDicer_fwd.wig",
             "/aligned_GRCz10/wigs/Bud_MZDicer_rev.wig", "/aligned_GRCz10/wigs/Bud_rev.wig")
  res <- findNGSPairs(paths)
  expect_equal(res$match[1], TRUE)
  # Compression
  paths <- paste0(paths, ".gz")
  res <- findNGSPairs(paths)
  expect_equal(res$match[1], TRUE)
  # No match
  paths[4] <- "/aligned_GRCz10/wigs/Bud1_rev.wig"
  res <- findNGSPairs(paths)
  expect_equal(res[4], paths[4])

  # Bed files
  paths <- c("1KCell_fwd.bed.gz", "2KCell_fwd.bed.gz",
             "2KCell_rev.bed.gz","1KCell_rev.bed.gz")
  res <- findNGSPairs(paths, format = "bed")
  expect_equal(res$match[1], TRUE)
})

test_that("detect_drive() picks the longest (most specific) matching mount, not just excludes root", {
  if (.Platform$OS.type != "unix") skip("detect_drive() is unix-only")
  # Portable (no mocking, no assumption about this machine's specific
  # layout): find whatever non-root mount is CURRENTLY deepest on this
  # machine and confirm a path nested under it resolves to that same
  # mount. Catches the exact bug this guards against: a path sitting
  # under mount A, itself nested inside a DIFFERENT, non-root mount B
  # (e.g. a second drive mounted as a subdirectory of an already
  # separate data mount) used to leave 2 candidates after only
  # excluding root, giving up with "Could not determine the correct
  # unique Filesystem drive" even though a single real answer exists.
  df_lines <- system("df -h", intern = TRUE)[-1]
  mounts <- vapply(df_lines, function(l) strsplit(trimws(l), " +")[[1]][6], character(1))
  non_root <- unique(mounts[mounts != "/"])
  skip_if(length(non_root) == 0, "no non-root mount available on this machine to test against")
  deepest <- non_root[which.max(nchar(non_root))]

  probe_path <- file.path(deepest, "a_subdir_that_neednt_exist")
  expect_false(is.na(detect_drive(probe_path)))
  expect_identical(detect_drive(probe_path), detect_drive(deepest))
})

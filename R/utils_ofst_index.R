# Sidecar chromosome-block index of an OFST file: "<file>.idx" (e.g. "all.ofst.idx").
#
# Written by export.ofst() and the merge (mergeLibs / .ofst_save_merge) from the table they
# just wrote, so finding the chromosome blocks of a file needs no scan of its seqnames column.
# It holds the format version, the row count, file size and column names of the OFST file it
# describes, and one row per chromosome run (chr, from, to, score sum, start sorted inside the
# run). Readers accept it only if row count, file size and columns match the file; otherwise
# (stale, missing, unreadable) they ignore it and scan. The fst file itself is unchanged, so
# older ORFik versions simply ignore the sidecar.

.ofst_index_version <- 1L

.ofst_index_path <- function(path) paste0(path, ".idx")

#' Write the chromosome-block index of an OFST file
#'
#' Internal. \code{dt} must be the table that was just written to \code{path} (all rows; it may
#' hold a subset of the columns, at least seqnames; start and score add the sort flag and sums).
#' Never fails the caller: errors become a warning and no (or no stale) index.
#' @param path OFST file path (already written)
#' @param dt the data.frame / data.table written to path
#' @return invisible TRUE if written
#' @noRd
.ofst_index_write <- function(path, dt) {
  idx_path <- .ofst_index_path(path)
  if (file.exists(idx_path)) unlink(idx_path)
  if (!isTRUE(getOption("ORFik.ofst_index", TRUE))) return(invisible(FALSE))
  ok <- tryCatch({
    s <- dt[["seqnames"]]
    if (!is.factor(s) || !file.exists(path)) return(invisible(FALSE))
    meta <- fst::metadata_fst(path)
    # dt must be the table in the file (it may hold a subset of the columns, e.g. seqnames only)
    if (as.numeric(meta$nrOfRows) != as.numeric(nrow(dt)) || !all(names(dt) %in% meta$columnNames))
      stop("table does not match the file")
    start <- dt[["start"]]
    score <- dt[["score"]]
    r <- ofst_block_index(s, if (is.integer(start)) start else NULL,
                          if (is.numeric(score)) score else NULL)
    lv <- levels(s)
    v <- r$value
    chr <- rep(NA_character_, length(v))
    ok_code <- !is.na(v) & v >= 1L & v <= length(lv)
    chr[ok_code] <- lv[v[ok_code]]
    blocks <- data.table::data.table(chr = chr, from = c(1, head(r$end, -1) + 1), to = r$end,
                                     score = if (is.numeric(score)) r$score else NA_real_,
                                     start_sorted = r$start_sorted)
    contiguous <- !anyNA(chr) && !anyDuplicated(chr)
    sorted_by <- if (contiguous) {
      if (nrow(blocks) && all(blocks$start_sorted)) c("seqnames", "start") else "seqnames"
    } else character()
    idx <- list(format = "ORFik.ofst.idx", version = .ofst_index_version, nrow = as.numeric(meta$nrOfRows),
                file_size = as.numeric(file.size(path)), columns = as.character(meta$columnNames), sorted_by = sorted_by,
                blocks = blocks)
    tmp <- paste0(idx_path, ".tmp", Sys.getpid())
    saveRDS(idx, tmp, compress = FALSE)
    if (!file.rename(tmp, idx_path)) { unlink(tmp); stop("cannot rename ", tmp) }
    TRUE
  }, error = function(e) {
    warning("OFST index not written for '", path, "': ", conditionMessage(e), call. = FALSE)
    FALSE
  })
  invisible(isTRUE(ok))
}

#' Read the chromosome-block index of an OFST file, if present and valid
#'
#' Internal. Returns NULL when the sidecar is missing, unreadable or stale (row count, file
#' size or columns differ from the file), so callers fall back to scanning the file.
#' @param path OFST file path
#' @return list(version, nrow, file_size, columns, sorted_by, blocks) or NULL
#' @noRd
.ofst_index_read <- function(path) {
  idx_path <- .ofst_index_path(path)
  if (!isTRUE(getOption("ORFik.ofst_index", TRUE)) || !file.exists(idx_path)) return(NULL)
  idx <- tryCatch(readRDS(idx_path), error = function(e) NULL)
  if (!is.list(idx) || !identical(idx$format, "ORFik.ofst.idx") ||
      !identical(idx$version, .ofst_index_version)) return(NULL)
  meta <- tryCatch(fst::metadata_fst(path), error = function(e) NULL)
  if (is.null(meta)) return(NULL)
  b <- idx$blocks
  valid <- isTRUE(as.numeric(meta$nrOfRows) == idx$nrow) &&
    isTRUE(as.numeric(file.size(path)) == idx$file_size) &&
    identical(as.character(meta$columnNames), as.character(idx$columns)) &&
    data.table::is.data.table(b) && all(c("chr", "from", "to") %in% names(b)) &&
    (if (nrow(b)) b$from[1] == 1 && b$to[nrow(b)] == idx$nrow &&
       all(b$from[-1] == b$to[-nrow(b)] + 1) && all(b$to >= b$from) else idx$nrow == 0)
  if (!valid) return(NULL)
  idx
}

#' Move an OFST index along with its file (after export to a temporary path + rename)
#' @noRd
.ofst_index_rename <- function(from, to) {
  src <- .ofst_index_path(from); dst <- .ofst_index_path(to)
  if (file.exists(dst)) unlink(dst)
  if (file.exists(src) && !file.rename(src, dst)) unlink(src)
  invisible(NULL)
}

# .ofst_chr_blocks() result from a valid index: NULL if a chromosome is split.
.ofst_index_blocks <- function(idx) {
  b <- idx$blocks
  if (anyNA(b$chr) || anyDuplicated(b$chr)) return(NULL)
  data.table::data.table(chr = as.character(b$chr), from = as.numeric(b$from), to = as.numeric(b$to))
}

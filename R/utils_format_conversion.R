#' Convert libraries to coverage RLEs
#'
#' Saved in folder "cov_RLE" relative to default libraries of experiment
#' @inheritParams outputLibs
#' @param in_files paths to input files, default pshifted files:
#' \code{filepath(df, "pshifted")} in ofst format
#' @param seqinfo SeqInfo object, default \code{seqinfo(findFa(df))}
#' @param weight integer, numeric or single length character. Default "score".
#' Use score column in loaded in_files.
#' @param split.by.strand logical, default TRUE, split into forward and reverse
#' strand RleList inside covRle object.
#' @return invisible(NULL), files saved to disc
convert_to_covRle <- function(df, in_files =  filepath(df, "pshifted"),
                              out_dir = file.path(libFolder(df), "cov_RLE"),
                              split.by.strand = TRUE,
                              split.by.readlength = FALSE,
                              seq_info = seqinfo(df), weight = "score",
                              verbose = TRUE) {
  if (length(in_files) != nrow(df))
    stop("'df' and 'in_files must have equal size!")
  lib_names <- remove.file_ext(filepath(df, "default", basename = TRUE))
  out_filepaths <- file.path(out_dir, lib_names)
  dir.create(out_dir, showWarnings = FALSE)
  if (verbose) message("-- Converting to Coverage RleList")
  if (verbose) message("Output to dir: ", out_dir)
  for (i in seq_along(out_filepaths)) {
    if (verbose) message("- Library: ", i)
    out_file <- out_filepaths[i]
    in_file <- in_files[i]
    if (split.by.readlength) {
      x <- fimport(in_file)
      all_readl_lengths <- readWidths(x)
      read_lengths <- sort(unique(all_readl_lengths))
      if (verbose) message("Readlength:", appendLF = FALSE)
      for (i in read_lengths) {
        if (verbose) message(", ", i, appendLF = FALSE)
        out_file_rl <- paste0(out_file, "_", i)
        export.cov(x = x[all_readl_lengths == i], file = out_file_rl,
                   split.by.strand = split.by.strand,
                   seqinfo = seq_info, weight = weight)
      }
      if (verbose) message(", All readlengths merged")
      export.cov(x = x, file = out_file,
                 split.by.strand = split.by.strand,
                 seqinfo = seq_info, weight = weight)
    } else {
      export.cov(x = fimport(in_file), file = out_file,
                 split.by.strand = split.by.strand,
                 seqinfo = seq_info, weight = weight)
    }
  }
  if (verbose) message("Done")
  return(invisible(NULL))
}

#' Convert to BigWig
#'
#' @inheritParams convert_to_covRle
#' @return invisible(NULL), files saved to disc
convert_to_bigWig <- function(df, in_files =  filepath(df, "pshifted"),
                              out_dir = file.path(libFolder(df), "bigwig"),
                              split.by.strand = TRUE,
                              split.by.readlength = FALSE,
                              seq_info = seqinfo(df), weight = "score",
                              is_pre_collapsed = FALSE, verbose = TRUE) {
  if (length(in_files) != nrow(df))
    stop("'df' and 'in_files must have equal size!")
  lib_names <- remove.file_ext(filepath(df, "default", basename = TRUE))
  out_filepaths <- file.path(out_dir, lib_names)
  dir.create(out_dir, showWarnings = FALSE)
  if (verbose) message("-- Converting to BigWig")
  if (verbose) message("Output to dir: ", out_dir)
  for (i in seq_along(out_filepaths)) {
    if (verbose) message("- Library: ", i)
    out_file <- out_filepaths[i]
    in_file <- in_files[i]
    if (split.by.readlength) {
      x <- fimport(in_file, chrStyle = seq_info)
      all_readl_lengths <- readWidths(x)
      read_lengths <- sort(unique(all_readl_lengths))
      if (verbose) message("Readlength:", appendLF = FALSE)
      for (i in read_lengths) {
        if (verbose) message(", ", i, appendLF = FALSE)
        out_file_rl <- paste0(out_file, "_", i)
        export.bigWig(x = x[all_readl_lengths == i], file = out_file_rl,
                      split.by.strand = split.by.strand,
                      seq_info = seq_info, is_pre_collapsed = is_pre_collapsed)
      }
      if (verbose) message(", All readlengths merged")
      export.bigWig(x = x, file = out_file,
                   split.by.strand = split.by.strand,
                   seq_info = seq_info, is_pre_collapsed = is_pre_collapsed)
    } else {
      export.bigWig(x = fimport(in_file), file = out_file,
                    split.by.strand = split.by.strand,
                    seq_info = seq_info, is_pre_collapsed = is_pre_collapsed)
    }
  }
  if (verbose) message("Done")
  return(invisible(NULL))
}

#' Convert to fstwig
#'
#' Will split files by chromosome for faster loading for now.
#' This feature might change in the future!
#' @inheritParams convert_to_covRle
#' @return invisible(NULL), files saved to disc
convert_to_fstWig <- function(df, in_files =  filepath(df, "pshifted"),
                              out_dir = file.path(libFolder(df), "bigwig"),
                              split.by.strand = TRUE,
                              split.by.readlength = FALSE,
                              seq_info = seqinfo(df), weight = "score",
                              is_pre_collapsed = FALSE, verbose = TRUE) {
  if (length(in_files) != nrow(df))
    stop("'df' and 'in_files must have equal size!")
  lib_names <- remove.file_ext(filepath(df, "default", basename = TRUE))
  out_filepaths <- file.path(out_dir, lib_names)
  dir.create(out_dir, showWarnings = FALSE)
  if (verbose) message("-- Converting to fstwig")
  if (verbose) message("Output to dir: ", out_dir)
  for (i in seq_along(out_filepaths)) {
    if (verbose) message("- Library: ", i)
    out_file <- out_filepaths[i]
    in_file <- in_files[i]
    export.fstwig(x = fimport(in_file, chrStyle = seq_info), file = out_file,
                  by.readlength = split.by.readlength,
                  seqinfo = seq_info)
  }
  if (verbose) message("Done")
  return(invisible(NULL))
}

RiboCrypt::multiOmicsPlot_ORFikExp

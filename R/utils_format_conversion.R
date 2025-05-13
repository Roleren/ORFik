#' Convert libraries to ofst
#'
#' Saved by default in folder "ofst" relative to default
#' libraries of experiment.
#' Speeds up loading of full files compared to bam by large margins.
#'
#' If you want to keep bam files loaded or faster conversion if
#' you already have them loaded, use ORFik::convertLibs instead
#' @inheritParams outputLibs
#' @param df an ORFik \code{\link{experiment}}, or NULL is allowed if
#' both in_files and out_dir is specified manually.
#' @param in_files paths to input files, default:
#'  \code{filepath(df, "default")} with bam format files.
#' @param out_dir paths to output files, default
#'  \code{file.path(libFolder(df), "cov_RLE")}.
#' @return invisible(NULL), files saved to disc
#' @export
#' @family lib_converters
#' @examples
#' df <- ORFik.template.experiment.zf()
#' ## Usually do default folder, here we use tmpdir
#' folder_to_save <- file.path(tempdir(), "ofst")
#' convert_bam_to_ofst(df, out_dir = folder_to_save)
#' fimport(file.path(folder_to_save, "ribo-seq.ofst"))
convert_bam_to_ofst <- function(df, in_files =  filepath(df, "default"),
                            out_dir = file.path(libFolder(df), "ofst"),
                            verbose = TRUE, strandMode = rep(0, length(in_files))) {
  stopifnot(all(is(strandMode, "numeric")))
  stopifnot(length(in_files) > 0)
  if (!is.null(df)) {
    stopifnot(is(df, "experiment"))
    if (length(in_files) != nrow(df))
      stop("'df' and 'in_files must have equal size!")
  }
  lib_names <- remove.file_ext(basename(in_files))
  out_filepaths <- file.path(out_dir, paste0(lib_names, ".ofst"))
  dir.create(out_dir, showWarnings = FALSE)
  if (verbose) message("-- Converting bam to ofst")
  if (verbose) message("Output to dir: ", out_dir)
  total_files <- length(out_filepaths)
  for (i in seq_along(out_filepaths)) {
    index <- paste0("(", i, "/", total_files, ")")
    if (verbose) message("- Library", index, ": ", lib_names[i])
    out_file <- out_filepaths[i]
    in_file <- in_files[i]
    export.ofst(readBam(in_file, strandMode = strandMode[i]), out_file)
  }
  if (verbose) message("Done")
  return(invisible(NULL))
}

#' Convert libraries to covRle
#'
#' Saved by default in folder "cov_RLE" relative to default
#' libraries of experiment
#' @inheritParams convert_bam_to_ofst
#' @param in_files paths to input files, default pshifted files:
#'  \code{filepath(df, "pshifted")} in ofst format
#' @param out_dir paths to output files, default
#'  \code{file.path(libFolder(df), "cov_RLE")}.
#' @param split.by.strand logical, default TRUE, split into forward and reverse
#' strand RleList inside covRle object.
#' @param split.by.readlength logical, default FALSE, split into files
#' for each readlength, defined by readWidths(x) for each file.
#' @param seq_info SeqInfo object, default \code{seqinfo(findFa(df))}
#' @param format chatacter, default "qs", alternative "rds". File format to
#' save R object.
#' @param weight integer, numeric or single length character. Default "score".
#' Use score column in loaded in_files.
#' @return invisible(NULL), files saved to disc
#' @export
#' @family lib_converters
#' @examples
#' df <- ORFik.template.experiment()[10,]
#' ## Usually do default folder, here we use tmpdir
#' folder_to_save <- file.path(tempdir(), "cov_RLE")
#' convert_to_covRle(df, out_dir = folder_to_save)
#' fimport(file.path(folder_to_save, "RFP_Mutant_rep2.covrds"))
convert_to_covRle <- function(df, in_files =  filepath(df, "pshifted"),
                              out_dir = file.path(libFolder(df), "cov_RLE"),
                              split.by.strand = TRUE,
                              split.by.readlength = FALSE,
                              seq_info = seqinfo(df), weight = "score",
                              format = "qs",
                              verbose = TRUE) {
  if (!is.null(df)) {
    stopifnot(is(df, "experiment"))
    if (length(in_files) != nrow(df))
      stop("'df' and 'in_files must have equal size!")
  }
  lib_names <- remove.file_ext(basename(in_files))
  out_filepaths <- file.path(out_dir, lib_names)
  dir.create(out_dir, showWarnings = FALSE)
  if (verbose) message("-- Converting to covRle objects")
  if (verbose) message("Output to dir: ", out_dir)
  for (i in seq_along(out_filepaths)) {
    if (verbose) message("- Library: ", lib_names[i])
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
                   seqinfo = seq_info, weight = weight, format = format)
      }
      if (verbose) message(", All readlengths merged")
      export.cov(x = x, file = out_file,
                 split.by.strand = split.by.strand,
                 seqinfo = seq_info, weight = weight, format = format)
    } else {
      export.cov(x = fimport(in_file), file = out_file,
                 split.by.strand = split.by.strand,
                 seqinfo = seq_info, weight = weight, format = format)
    }
  }
  if (verbose) message("Done")
  return(invisible(NULL))
}

#' Convert libraries to covRleList objects
#'
#' Useful to store reads separated by readlength, for much faster
#' coverage calculation.
#' Saved by default in folder "cov_RLE_List" relative to default
#' libraries of experiment
#' @inheritParams convert_to_covRle
#' @param out_dir paths to output files, default
#' \code{file.path(libFolder(df), "cov_RLE_List")}.
#' @param out_dir_merged character vector of paths, default:
#'  \code{file.path(libFolder(df), "cov_RLE")}.
#'  Paths to merged output files, Set to NULL to skip making merged covRle.
#' @return invisible(NULL), files saved to disc
#' @export
#' @family lib_converters
#' @examples
#' df <- ORFik.template.experiment()[10,]
#' ## Usually do default folder, here we use tmpdir
#' folder_to_save <- file.path(tempdir(), "cov_RLE_List")
#' folder_to_save_merged <- file.path(tempdir(), "cov_RLE")
#' ORFik:::convert_to_covRleList(df, out_dir = folder_to_save,
#' out_dir_merged = folder_to_save_merged)
#' fimport(file.path(folder_to_save, "RFP_Mutant_rep2.covrds"))
convert_to_covRleList <- function(df, in_files =  filepath(df, "pshifted"),
                              out_dir = file.path(libFolder(df), "cov_RLE_List"),
                              out_dir_merged = file.path(libFolder(df), "cov_RLE"),
                              split.by.strand = TRUE,
                              seq_info = seqinfo(df), weight = "score",
                              format = format,
                              verbose = TRUE) {
  if (!is.null(df)) {
    stopifnot(is(df, "experiment"))
    if (length(in_files) != nrow(df))
      stop("'df' and 'in_files must have equal size!")
  }
  lib_names <- remove.file_ext(basename(in_files))
  out_filepaths <- file.path(out_dir, lib_names)
  out_filepaths_merged <- file.path(out_dir_merged, lib_names)
  dir.create(out_dir, showWarnings = FALSE)
  if (!is.null(out_dir_merged)) dir.create(out_dir_merged, showWarnings = FALSE)
  if (verbose) message("-- Converting to covRleList objects")
  if (verbose) message("Output to dir: ", out_dir)
  for (i in seq_along(out_filepaths)) {
    if (verbose) message("- Library: ", lib_names[i])
    out_file <- out_filepaths[i]
    in_file <- in_files[i]
    x <- fimport(in_file)
    export.covlist(x, file = out_file,
                   split.by.strand = split.by.strand,
                   seqinfo = seq_info, weight = weight,
                   verbose = verbose, format = format)

    if (!is.null(out_dir_merged)) {
      if (verbose) message(", All readlengths merged")
      export.cov(x = x, file = out_filepaths_merged[i],
                 split.by.strand = split.by.strand,
                 seqinfo = seq_info, weight = weight, format = format)
    }

  }
  if (verbose) message("Done")
  return(invisible(NULL))
}

#' Convert to BigWig
#'
#' @inheritParams convert_to_covRle
#' @inheritParams export.bigWig
#' @param out_dir paths to output files, default
#'  \code{file.path(libFolder(df), "bigwig")}.
#' @return invisible(NULL), files saved to disc
#' @export
#' @family lib_converters
#' @examples
#' df <- ORFik.template.experiment()[10,]
#' ## Usually do default folder, here we use tmpdir
#' folder_to_save <- file.path(tempdir(), "bigwig")
#' convert_to_bigWig(df, out_dir = folder_to_save)
#' fimport(file.path(folder_to_save, c("RFP_Mutant_rep2_forward.bigWig",
#'  "RFP_Mutant_rep2_reverse.bigWig")))
convert_to_bigWig <- function(df, in_files =  filepath(df, "pshifted"),
                              out_dir = file.path(libFolder(df), "bigwig"),
                              split.by.strand = TRUE,
                              split.by.readlength = FALSE,
                              seq_info = seqinfo(df), weight = "score",
                              is_pre_collapsed = FALSE, verbose = TRUE) {
  if (!is.null(df)) {
    stopifnot(is(df, "experiment"))
    if (length(in_files) != nrow(df))
      stop("'df' and 'in_files must have equal size!")
  }
  lib_names <- remove.file_ext(basename(in_files))
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
#' @inheritParams convert_to_bigWig
#' @param out_dir paths to output files, default
#'  \code{file.path(libFolder(df), "bigwig")}.
#' @return invisible(NULL), files saved to disc
convert_to_fstWig <- function(df, in_files =  filepath(df, "pshifted"),
                              out_dir = file.path(libFolder(df), "fstwig"),
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
                  by.readlength = split.by.readlength)
  }
  if (verbose) message("Done")
  return(invisible(NULL))
}

#' Convert libraries to ofst
#'
#' Saved by default in folder "ofst" relative to default
#' libraries of experiment.
#' Speeds up loading of full files compared to bam by large margins.
#'
#' If you want to keep bam files loaded or faster conversion if
#' you already have them loaded, use ORFik::convertLibs instead
#' @inheritParams outputLibs
#' @param read_strand Single-end RNA orientation passed to readBam. Auto uses
#' the validated per-BAM decision; existing paired strandMode semantics remain.
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
#' folder_to_save_unique <- file.path(tempdir(), "unique_mappers", "ofst")
#' convert_bam_to_ofst(df, out_dir = folder_to_save_unique, only_unique_mappers = TRUE)
#' fimport(file.path(folder_to_save_unique, "ribo-seq.ofst"))
convert_bam_to_ofst <- function(df,
                            in_files =  filepath(df, "default"),
                            out_dir = file.path(libFolder(df, unique_mappers = only_unique_mappers), "ofst"),
                            verbose = TRUE, strandMode = rep(0, length(in_files)),
                            only_unique_mappers = uniqueMappers(df), read_strand = "auto") {
  stopifnot(all(is(strandMode, "numeric")))
  if (is.character(in_files) && length(in_files) == 2L && identical(in_files[2], "paired-end")) in_files <- list(in_files)
  sources <- if (is.list(in_files)) in_files else as.list(in_files)
  stopifnot(length(sources) > 0)
  if (!is.null(df)) {
    stopifnot(is(df, "experiment"))
    if (length(sources) != nrow(df))
      stop("'df' and 'in_files must have equal size!")
  }
  read_strand <- rep_len(read_strand, length(sources))
  lib_names <- vapply(sources, function(path) remove.file_ext(basename(path[1])), character(1))
  out_filepaths <- file.path(out_dir, paste0(lib_names, ".ofst"))
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  if (verbose) message("-- Converting bam to ofst")
  if (verbose) message("Output to dir: ", out_dir)
  total_files <- length(out_filepaths)
  for (i in seq_along(out_filepaths)) {
    index <- paste0("(", i, "/", total_files, ")")
    if (verbose) message("- Library", index, ": ", lib_names[i])
    out_file <- out_filepaths[i]
    in_file <- sources[[i]]
    source_file <- in_file[1]
    .publish_bam_ofst(readBam(in_file, strandMode = strandMode[i], only_unique_mappers = only_unique_mappers, read_strand = read_strand[i]),
                      out_file, source_file, read_strand[i])
  }
  if (verbose) message("Done")
  return(invisible(NULL))
}

# Publish each mapper mode independently, preserving strand sidecars and atomic OFST writes.
.publish_bam_ofst <- function(reads, out_file, source_file, read_strand) {
  out_dir <- dirname(out_file)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  temporary <- tempfile(pattern = ".rna_ofst_", tmpdir = out_dir, fileext = ".ofst")
  on.exit(unlink(temporary))
  export.ofst(reads, temporary)
  if (!file.rename(temporary, out_file)) stop("Cannot publish complete OFST: ", out_file)
  decision <- rna_strand_record(source_file)
  orientation <- rna_strand_resolve(source_file, read_strand)
  if (!is.null(decision) || orientation != "alignment") {
    if (is.null(decision)) decision <- list(version = 1L, status = "explicit", read_strand = orientation)
    decision$read_strand <- orientation
    decision$normalized <- TRUE
    st <- file.info(out_file)
    decision$product_signature <- c(size = st$size, mtime = as.numeric(st$mtime))
    sidecar <- rna_strand_sidecar(out_file)
    dir.create(dirname(sidecar), recursive = TRUE, showWarnings = FALSE)
    saveRDS(decision, sidecar)
  } else if (file.exists(rna_strand_sidecar(out_file))) unlink(rna_strand_sidecar(out_file))
  invisible(out_file)
}

.ofst_read_length_distribution <- function(reads) {
  if (!length(reads)) return(data.table::data.table(read_length = integer(), count = numeric(), percent = numeric()))
  score <- S4Vectors::mcols(reads)$score
  if (is.null(score)) score <- rep(1, length(reads))
  dt <- data.table::data.table(read_length = readWidths(reads), score = score)
  result <- dt[, list(count = sum(score)), by = "read_length"]
  data.table::setorder(result, read_length)
  result[["percent"]] <- round(100 * result[["count"]] / sum(result[["count"]]), 2)
  result
}

#' Convert one BAM sample to requested OFST mapper modes with one alignment import
#'
#' For an ordinary single-BAM alignment representation, NH tags are imported
#' with the alignments and used before exporting unique mappers. RNA strand
#' normalization, selected-mate handling, weights and atomic output sidecars
#' are identical to the ordinary converter. Complex paired inputs fall back
#' to the ordinary converter. Successful modes invoke after_mode immediately,
#' so callers can record independent resume flags even if a later mode fails.
#' @param df One-row ORFik experiment.
#' @param modes Character subset of c("all", "unique").
#' @param out_dirs Named output directories for the requested modes.
#' @param read_strand Orientation selection, default "auto".
#' @param after_mode Optional function(mode, path, distribution), invoked after
#' publishing each output. distribution is a score-weighted length table,
#' or NULL for complex paired input (caller may import the published file).
#' @return Invisible named output paths.
#' @export
convert_bam_to_ofst_modes <- function(df, modes = c("all", "unique"),
    out_dirs = NULL, read_strand = "auto", after_mode = NULL) {
  stopifnot(is(df, "experiment"), nrow(df) == 1L,
            length(modes) > 0L, !anyDuplicated(modes), all(modes %in% c("all", "unique")),
            length(read_strand) == 1L, is.null(after_mode) || is.function(after_mode))
  input <- filepath(df, "default")
  if (is.null(out_dirs)) out_dirs <- setNames(vapply(modes, function(mode)
    file.path(libFolder(df, unique_mappers = mode == "unique"), "ofst"), character(1)), modes)
  stopifnot(all(modes %in% names(out_dirs)))
  simple <- is.character(input) && length(input) == 1L && grepl("\\.bam$", input, ignore.case = TRUE)
  source <- if (is(input, "data.table")) input$forward[1] else if (is.list(input)) input[[1]][1] else input[1]
  paths <- setNames(file.path(out_dirs[modes], paste0(remove.file_ext(basename(source)), ".ofst")), modes)
  if (!simple) {
    for (mode in modes) {
      convert_bam_to_ofst(df, out_dir = out_dirs[[mode]], only_unique_mappers = mode == "unique", read_strand = read_strand)
      if (!is.null(after_mode)) after_mode(mode, paths[[mode]], NULL)
    }
    return(invisible(paths))
  }
  param <- if ("unique" %in% modes) Rsamtools::ScanBamParam(tag = "NH") else NULL
  reads <- readBam(input, param = param, read_strand = read_strand)
  unique_index <- NULL; unique_nh <- NULL
  decision <- rna_strand_record(input)
  keep_unique_nh <- !is.null(decision) && identical(decision$selected_mate, 1L)
  if ("unique" %in% modes) {
    nh <- if (!length(reads)) integer() else S4Vectors::mcols(reads)$NH
    if (is.null(nh) || length(nh) != length(reads) || anyNA(nh))
      stop("Unique mapper conversion requires complete NH tags")
    unique_index <- which(nh == 1L)
    if (keep_unique_nh) unique_nh <- nh[unique_index]
    S4Vectors::mcols(reads)$NH <- NULL
  }
  # Always publish all first if requested, releasing its full alignment view
  # before constructing the unique subset. Do not hold two large copies.
  for (mode in intersect(c("all", "unique"), modes)) {
    if (mode == "unique") {
      reads <- reads[unique_index]
      # The established selected-mate reader retains NH in unique output.
      if (keep_unique_nh) S4Vectors::mcols(reads)$NH <- unique_nh
    }
    .publish_bam_ofst(reads, paths[[mode]], input, read_strand)
    if (!is.null(after_mode)) after_mode(mode, paths[[mode]], tryCatch(.ofst_read_length_distribution(reads), error = function(e) NULL))
  }
  invisible(paths)
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
#' fimport(file.path(folder_to_save, "RFP_Mutant_rep2.covqs"))
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
                   split.by.strand = split.by.strand && rna_strand_resolve(in_file) != "unstranded",
                   seqinfo = seq_info, weight = weight, format = format)
      }
      if (verbose) message(", All readlengths merged")
      export.cov(x = x, file = out_file,
                 split.by.strand = split.by.strand && rna_strand_resolve(in_file) != "unstranded",
                 seqinfo = seq_info, weight = weight, format = format)
    } else {
      export.cov(x = fimport(in_file), file = out_file,
                 split.by.strand = split.by.strand && rna_strand_resolve(in_file) != "unstranded",
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
#' fimport(file.path(folder_to_save, "RFP_Mutant_rep2.covqs"))
convert_to_covRleList <- function(df, in_files =  filepath(df, "pshifted"),
                              out_dir = file.path(libFolder(df), "cov_RLE_List"),
                              out_dir_merged = file.path(libFolder(df), "cov_RLE"),
                              split.by.strand = TRUE,
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
                   split.by.strand = split.by.strand && rna_strand_resolve(in_file) != "unstranded",
                   seqinfo = seq_info, weight = weight,
                   verbose = verbose, format = format)

    if (!is.null(out_dir_merged)) {
      if (verbose) message(", All readlengths merged")
      export.cov(x = x, file = out_filepaths_merged[i],
                 split.by.strand = split.by.strand && rna_strand_resolve(in_file) != "unstranded",
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
                      split.by.strand = split.by.strand && rna_strand_resolve(in_file) != "unstranded",
                      seq_info = seq_info, is_pre_collapsed = is_pre_collapsed)
      }
      if (verbose) message(", All readlengths merged")
      export.bigWig(x = x, file = out_file,
                   split.by.strand = split.by.strand && rna_strand_resolve(in_file) != "unstranded",
                   seq_info = seq_info, is_pre_collapsed = is_pre_collapsed)
    } else {
      export.bigWig(x = fimport(in_file), file = out_file,
                    split.by.strand = split.by.strand && rna_strand_resolve(in_file) != "unstranded",
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

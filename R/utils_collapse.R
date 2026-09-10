#' Convert a GRanges Object to 1 width reads
#'
#' There are 5 ways of doing this\cr
#' 1. Take 5' ends, reduce away rest (5prime)\cr
#' 2. Take 3' ends, reduce away rest (3prime)\cr
#' 3. Tile to 1-mers and include all (tileAll)\cr
#' 4. Take middle point per GRanges (middle)\cr
#' 5. Get original with metacolumns (None)\cr
#' You can also do multiple at a time, then output is GRangesList, where
#' each list group is the operation (5prime is [1], 3prime is [2] etc)\cr
#' Many other ways to do this have their own functions, like startSites and
#' stopSites etc.
#' To retain information on original width, set addSizeColumn to TRUE.
#' To compress data, 1 GRanges object per unique read, set addScoreColumn to
#' TRUE. This will give you a score column with how many duplicated reads there
#' were in the specified region.
#'
#' NOTE: Note: For cigar based ranges (GAlignments),
#' the 5' end is the first non clipped base (neither soft clipped or hard clipped
#' from 5'). This is following the default
#' of bioconductor.
#' For special case of GAlignmentPairs, 5prime will only use left (first)
#' 5' end and read and 3prime will use only right (last) 3' end of read
#' in pair. tileAll and middle can possibly find poinst that are not in the
#' reads since: lets say pair is 1-5 and 10-15, middle is 7, which is not in
#' the read.
#'
#' @param gr GRanges, GAlignment or GAlignmentPairs object to reduce.
#' @param method character, default \code{"5prime"},
#' the method to reduce ranges, see NOTE for more info.
#' @param addScoreColumn logical (FALSE), if TRUE, add a score column that
#'  sums up the hits per unique range. This will make each read unique, so
#'  that each read is 1 time, and score column gives the number of
#'  collapsed hits.
#'  A useful compression. If addSizeColumn is FALSE, it will not differentiate
#'  between reads with same start and stop, but different length. If
#'  addSizeColumn is FALSE, it will remove it. Collapses after conversion.
#' @param addSizeColumn logical (FALSE), if TRUE, add a size column that
#'  for each read, that gives original width of read. Useful if you need
#'  original read lengths. This takes care of soft clips etc.
#'  If collapsing reads, each unique range will be grouped also by size.
#' @param reuse.score.column logical (TRUE), if addScoreColumn is TRUE,
#'  and a score column exists, will sum up the scores to create a new score.
#'  If FALSE, will skip old score column and create new according to number
#'  of replicated reads after conversion.
#'  If addScoreColumn is FALSE, this argument is ignored.
#' @inheritParams readWidths
#' @importFrom GenomicAlignments first
#' @importFrom GenomicAlignments last
#' @return Converted GRanges object
#' @export
#' @family utils
#' @examples
#' gr <- GRanges("chr1", 1:10,"+")
#' # 5 prime ends
#' convertToOneBasedRanges(gr)
#' # is equal to convertToOneBasedRanges(gr, method = "5prime")
#' # 3 prime ends
#' convertToOneBasedRanges(gr, method = "3prime")
#' # With lengths
#' convertToOneBasedRanges(gr, addSizeColumn = TRUE)
#' # With score (# of replicates)
#' gr <- rep(gr, 2)
#' convertToOneBasedRanges(gr, addSizeColumn = TRUE, addScoreColumn = TRUE)
#'
convertToOneBasedRanges <- function(gr, method = "5prime",
                                    addScoreColumn = FALSE,
                                    addSizeColumn = FALSE,
                                    after.softclips = TRUE,
                                    along.reference = FALSE,
                                    reuse.score.column = TRUE) {
  if (addSizeColumn & is.null(mcols(gr)$size)) {
    mcols(gr) <- S4Vectors::DataFrame(mcols(gr),
                                      size = readWidths(gr, after.softclips,
                                                        along.reference))
  }
  # Convert to positions wanted
  if (!is(gr, "GRanges")) gr <- GRanges(gr)
  if (method == "5prime") {
    gr <- resize(gr, width = 1, fix = "start")
  } else if(method == "3prime") {
    gr <- resize(gr, width = 1, fix = "end")
  } else if(method %in% c("None", "none")) {
  } else if(method == "tileAll") {
    gr <- unlist(tile(gr, width = 1), use.names = FALSE)
  } else if (method == "middle") {
    ranges(gr) <- IRanges(start(gr) + ceiling((end(gr) - start(gr)) / 2),
                          width = 1)
  } else stop("invalid type: must be 5prime, 3prime, None, tileAll or middle")
  # Collapse after conversion
  if (addScoreColumn) {
    gr <- collapseDuplicatedReads(gr, addSizeColumn = addSizeColumn,
                                  reuse.score.column = reuse.score.column)
  }
  return(gr)
}

#' Merge reads by sum of existing scores
#'
#' If you have multiple reads a same location but different read lengths,
#' specified in meta column "size", it will sum up the scores
#' (number of replicates) for all reads at that position
#' @param x a GRanges object
#' @return merged GRanges object
#' @keywords internal
#' @examples
#' gr_s1 <- rep(GRanges("chr1", 1:10,"+"), 2)
#' gr_s2 <- GRanges("chr1", 1:12,"+")
#' gr2 <- GRanges("chr1", 21:40,"+")
#' gr <- c(gr_s1, gr_s2, gr2)
#' res <- convertToOneBasedRanges(gr,
#'    addScoreColumn = TRUE, addSizeColumn = TRUE)
#' ORFik:::collapse.by.scores(res)
#'
collapse.by.scores <- function(x) {
  dt <- data.table(seqnames = as.character(seqnames(x)),
                   start = start(ranges(x)),
                   end = end(ranges(x)),
                   strand = as.character(strand(x)),
                   score = mcols(x)$score)
  dt <- dt[, .(score = sum(score)), .(seqnames, start, end, strand)]
  # TODO change makeGRangesFromDataFrame to internal fast function
  return(makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE, seqinfo = seqinfo(x)))
}

#' Merge multiple ofst file
#'
#' Collapses and sums the score column of each ofst file
#' It is required that each file is of same ofst type.
#' That is if one file has cigar information, all must have it.
#' @param file_paths Full path to .ofst files wanted to merge
#' @param lib_names character, the name to give the resulting score columns.
#' Default: \code{sub(pattern = "\\.ofst$", replacement = "", basename(file_paths))}
#' @param keep_all_scores logical, default TRUE, keep all library scores in the merged file. These
#' score columns are named the libraries full name from \code{bamVarName(df)}.
#' @param sort logical, default TRUE. Sort the ranges. Will make the file smaller and
#' faster to load, but some additional merging time is added.
#' @param keepCigar logical, default TRUE. If CIGAR is defined, keep column. Setting
#' to FALSE compresses the file much more usually.
#' @param max_splits integer, default 20. If number of rows to merge > 2^31,
#' how many times can you allow split merging to try to "rescue" the merging
#' process?
#' @param dt_max_index_size 2^31, the number of rows data.table support, set lower
#' to merge split on lower counts in the ordinary lossless path. This groups
#' whole input files and is independent of the final output target.
#' @param allow_filtering logical, default TRUE. If the distinct merged rows
#' exceed \code{filter_target_rows}, discard whole chromosome/start positions,
#' lowest pooled score first. FALSE aborts without returning a filtered result.
#' @param filter_target_rows numeric whole number, default \code{2^31 - 2}.
#' Maximum final alignment rows, strictly below data.table's 32-bit ceiling.
#' Zero is allowed. Dropping whole positions can retain fewer rows than this.
#' @param filter_seed integer, default 1. Seed for deterministic pseudo-random
#' priorities among equally scored positions across all chromosomes. Does not
#' change R's RNG state; independent of input order and scratch partitioning.
#' @param max_filter_score numeric, default Inf. Only positions whose pooled
#' score is at most this ceiling may be removed. Abort if too few are eligible.
#' @param filter_chunk_rows numeric whole number, default 5e6. Maximum raw rows
#' per input batch and rows per compacted partition. Batches are collapsed
#' before writing, then compacted in balanced rounds while inputs are read.
#' Pooled-score output collapses across studies immediately. Linux available
#' RAM and sampled uncompressed sizes may reduce this budget with an eightfold
#' workspace allowance; this estimate is not a hard byte limit. Temporary
#' combines contain up to two partitions. A single position's distinct rows
#' (per input when keeping library scores) must fit the effective budget.
#' Disk usage depends on distinct merged data, not just the batch size.
#' @param filter_tmpdir character, default \code{tempdir()} for \code{ofst_merge}
#' and \code{out_dir} for \code{mergeLibs}. Existing writable
#' directory on a disk with space for scratch partitions and summaries. Only
#' a newly created private subdirectory is removed on completion or error.
#' @param filter_fallback_dir character or NULL. An existing alternate directory
#' used for one automatic restart after a scratch-storage failure. Default NULL
#' for \code{ofst_merge}; \code{mergeLibs} defaults to its \code{out_dir}. The
#' input files are read again with unchanged filtering settings. NULL disables
#' fallback. Invalid inputs, memory errors and user interrupts are not retried.
#' A small ownership cache here also records primary scratch for crash cleanup.
#' @param filter_input_summary logical, default FALSE. For pooled-score output,
#' reread inputs after selection to reconstruct exact per-input removal counts.
#' This adds disk I/O and bounded per-input compaction. When FALSE these counts
#' are NA, not zero; raw input counts and exact global/chromosome totals remain
#' available. When keeping all library scores, per-input counts are available
#' without another pass regardless of this flag. Only used if filtering occurs.
#' @details Chromosome and strand columns may be character or factor, including
#' different factor levels in different files. Factors are retained internally
#' to reduce memory use. Coordinates must be integer-valued and fit in a 32-bit
#' integer. Missing scores contribute zero, as in the per-library score mode.
#' All alignment columns except an explicitly discarded CIGAR remain merge keys.
#' If raw input rows exceed the final target, files (including a single large
#' file) are merged in bounded batches before counting distinct
#' merged alignment keys. No filtering occurs if duplicate collapse is enough.
#' Positional filtering requires single-end-style \code{seqnames} and
#' \code{start}; missing positions and negative/infinite scores are rejected
#' in this rescue path. It sums scores over all original files, both strands
#' and all CIGARs at each position. Retained alignment keys are not collapsed
#' to positions. This is a lossy size rescue, not a biological QC filter.
#' Equal pooled scores are ordered by seeded 53-bit hash priorities, with
#' chromosome/start used only for actual priority collisions. Whole groups
#' are removed until the target is met; unequal group sizes mean this is not
#' an equal-probability sample of individual alignment rows. Integer count
#' sums are exact up to double precision's integer limit; fractional sums may
#' exhibit ordinary floating-point rounding effects across partitionings.
#' Scratch processing bounds intermediate partitions, but the final table
#' and its assembly still need to fit in RAM. Both score modes use the same
#' removal decision. Filtering is not performed separately on merge chunks.
#' Ordinary R errors and interrupts clean this run's private caches. SIGKILL,
#' OOM kills and machine crashes cannot execute R cleanup handlers. On a later
#' run, positively identified abandoned caches in the configured directories
#' can be removed, including primary scratch registered by the output-side
#' cache. Linux process identity checks include host, user, boot, PID namespace
#' and process start time; active, foreign or unverifiable owners are left alone.
#' Unmarked folders and symlinks are never claimed for automatic cleanup.
#' Final outputs and input files are never deleted by cache cleanup. Recovery
#' files from interrupted output publication are retained with a message.
#' @return a data.table of merged result, it is merged on all columns except "score".
#' The returned file will contain the scores of each file + the aggregate sum score.
#' Only when positions were actually removed, \code{attr(result, "removal_summary")}
#' is an \code{ofst_removal_summary} list with parameters, score/priority cutoffs,
#' before/after/removed alignment rows, positions and total score, and
#' \code{by_chromosome} / \code{by_input} data.tables. Per-input rows count
#' distinct original alignment keys within each file, so their sum can exceed
#' global merged rows. In pooled mode these counts require
#' \code{filter_input_summary = TRUE}; otherwise they are NA and
#' \code{input_summary_computed} is FALSE. Raw input row counts are always
#' recorded separately. Schema version 2 also records the batch strategy. The
#' summary is absent for unfiltered results. FST does not preserve attributes;
#' \code{mergeLibs} saves this object in a companion RDS file.
#' @importFrom data.table setnames
ofst_merge <- function(file_paths,
                       lib_names = sub("\\.ofst$", "", basename(file_paths)),
                       keep_all_scores = TRUE, keepCigar = TRUE, sort = TRUE,
                       max_splits = 20L, dt_max_index_size = 2^31,
                       allow_filtering = TRUE, filter_target_rows = 2^31 - 2,
                       filter_seed = 1L, max_filter_score = Inf,
                       filter_chunk_rows = 5e6, filter_tmpdir = tempdir(),
                       filter_fallback_dir = NULL, filter_input_summary = FALSE) {
  restore_rng <- .ofst_rng_restore()
  on.exit(restore_rng(), add = TRUE)
  if (!is.character(file_paths) || !is.null(dim(file_paths)) || !length(file_paths) ||
      anyNA(file_paths) || any(!nzchar(file_paths)))
    stop("file_paths must be a non-empty character vector without missing or empty paths.")
  missing_paths <- file_paths[!file.exists(file_paths) | dir.exists(file_paths)]
  if (length(missing_paths)) .ofst_abort("input files do not exist or are directories: ", paste(missing_paths, collapse = ", "))
  if (length(file_paths) != length(lib_names))
    .ofst_abort("lib_names must have one name per input file (", length(file_paths), ").")
  if (!is.character(lib_names) || !is.null(dim(lib_names)) ||
      anyNA(lib_names) || any(!nzchar(lib_names)))
    stop("lib_names must contain non-missing, non-empty character names.")
  for (arg in c("keep_all_scores", "keepCigar", "sort")) {
    value <- get(arg)
    if (!is.logical(value) || !is.null(dim(value)) || length(value) != 1L || is.na(value))
      stop(arg, " must be TRUE or FALSE.")
  }
  .ofst_filter_controls(allow_filtering, filter_target_rows, filter_seed,
                        max_filter_score, filter_chunk_rows, filter_tmpdir, filter_fallback_dir, filter_input_summary)
  .plan_splits(0, limit = dt_max_index_size, max_splits = max_splits)
  columns <- .validate_schema(file_paths)
  if (keep_all_scores && (anyDuplicated(lib_names) || any(lib_names %in% columns)))
    stop("lib_names must be unique and must not collide with OFST column names.")

  # Plan splits vs data.table's index limit
  row_numbers <- vapply(file_paths, function(x) as.double(fst::metadata_fst(x)$nrOfRows), 0)
  if (sum(row_numbers) > filter_target_rows) {
    controls <- list(allow_filtering = allow_filtering, filter_target_rows = filter_target_rows,
                     filter_seed = filter_seed, max_filter_score = max_filter_score,
                     filter_chunk_rows = filter_chunk_rows, filter_tmpdir = filter_tmpdir,
                     filter_fallback_dir = filter_fallback_dir, filter_input_summary = filter_input_summary)
    # Every original input participates in one global decision. Never filter
    # first-round chunks independently, even when their raw row sum is large.
    keys <- setdiff(fst::metadata_fst(file_paths[1L])$columnNames,
                    c("score", if (!keepCigar) "cigar"))
    dt <- .ofst_merge_filtered(file_paths, lib_names, row_numbers, keys,
                               keep_all_scores, controls)
    if (sort) .ofst_sort(dt, setdiff(names(dt), "score"))
    return(.ofst_finalize(dt))
  }
  plan <- .plan_splits(row_numbers, limit = dt_max_index_size,
                       max_splits = max_splits)

  # First-round: per-chunk merge
  file_paths_split <- split(file_paths, plan$split_vector)
  lib_names_split  <- split(lib_names,  plan$split_vector)

  is_single_chunked <- plan$is_single_chunked
  if (!is_single_chunked) message(plan$note_message)
  message(plan$message)
  sort_on_chunking <- is_single_chunked & sort
  merge_chunk <- function(g) {
    message("- Merging chunk ", g, "/", length(file_paths_split))
    .ofst_merge_lossless(
      .read_fst_list(file_paths_split[[g]]),
      lib_names = lib_names_split[[g]],
      keep_all_scores = keep_all_scores,
      keepCigar = keepCigar,
      sort = sort_on_chunking
    )
  }

  # Second-round: combine chunks
  dt <- if (is_single_chunked) {
    merge_chunk(1L)
  } else {
    # Return the list directly to the reducer: retaining chunk_list in this
    # frame would keep consumed inputs alive throughout the second round.
    merge_chunks <- function() {
      chunks <- lapply(seq_along(file_paths_split), merge_chunk)
      message("Split round 2")
      chunks
    }
    .ofst_merge_lossless(
      merge_chunks(),
      lib_names = lib_names,
      keep_all_scores = keep_all_scores,
      keepCigar = keepCigar,
      sort = sort, chunkified = TRUE
    )
  }
  .ofst_finalize(dt)
}

.ofst_finalize <- function(dt) {
  # Make seqnames and strand factor
  # Match the public output's factor levels without expanding N character
  # pointers. tabulate only allocates a level-sized count vector.
  seq_levels <- levels(dt$seqnames)
  used <- tabulate(dt$seqnames, nbins = length(seq_levels)) > 0L
  data.table::set(dt, j = "seqnames",
                  value = .ofst_relevel(dt$seqnames, sort(seq_levels[used])))
  data.table::set(dt, j = "strand",
                  value = .ofst_relevel(dt$strand, c("+", "-", "*")))
  message("Done merging")
  dt[]
}

.ofst_factor <- function(x, column) {
  if (is.character(x)) return(factor(x))
  if (!is.factor(x)) stop("OFST column '", column, "' must be character or factor.")
  lv <- levels(x)
  # A range check avoids allocating several N-row logical masks on every
  # ordinary factor. Only malformed inputs need the slower repair path.
  bounds <- suppressWarnings(range(as.integer(x), na.rm = TRUE))
  if (bounds[1L] < 1L || bounds[2L] > length(lv) || anyNA(lv) || anyDuplicated(lv)) {
    codes <- as.integer(x)
    codes[is.na(codes) | codes < 1L | codes > length(lv)] <- NA_integer_
    return(factor(lv[codes]))
  }
  if (is.ordered(x)) return(structure(as.integer(x), levels = lv, class = "factor"))
  x
}

.ofst_relevel <- function(x, new_levels) {
  if (identical(levels(x), new_levels)) return(x)
  structure(match(levels(x), new_levels)[as.integer(x)],
            levels = new_levels, class = "factor")
}

.ofst_sort <- function(dt, columns) {
  # data.table sorts factors by codes, but OFST historically sorted their
  # character labels. Recode only the factor keys, using a small level table
  # to preserve data.table's character ordering independently of input levels.
  for (col in columns) {
    if (!is.factor(dt[[col]])) next
    labels <- data.table::data.table(label = levels(dt[[col]]))
    data.table::setorderv(labels, "label")
    data.table::set(dt, j = col, value = .ofst_relevel(dt[[col]], labels$label))
  }
  data.table::setorderv(dt, columns)
  dt
}

.validate_ofst_table <- function(d) {
  if (!data.table::is.data.table(d)) stop("OFST input must be a data.table.")
  if (anyDuplicated(names(d))) stop("OFST column names must be unique.")
  coordinates <- if (any(c("cigar1", "cigar2", "start1", "start2") %in% names(d)))
    c("start1", "start2", "cigar1", "cigar2") else "start"
  missing <- setdiff(c("seqnames", coordinates, "strand", "score"), names(d))
  if (length(missing)) stop("Missing required OFST columns: ", paste(missing, collapse = ", "))
  for (col in intersect(c("seqnames", "strand", "cigar", "cigar1", "cigar2"), names(d))) {
    if ((!is.factor(d[[col]]) && !is.character(d[[col]])) || !is.null(dim(d[[col]])))
      stop("OFST column '", col, "' must be character or factor.")
  }
  for (col in intersect(c("start", "start1", "start2", "end", "size", "width", "score"), names(d))) {
    if (!is.numeric(d[[col]]) || is.object(d[[col]]) || !is.null(dim(d[[col]])))
      stop("OFST column '", col, "' must be integer or numeric.")
  }
  invisible(d)
}

.normalize_dt <- function(d) {
  .validate_ofst_table(d)
  for (col in c("seqnames", "strand"))
    data.table::set(d, j = col, value = .ofst_factor(d[[col]], col))
  strand_levels <- levels(d$strand)
  used <- tabulate(d$strand, nbins = length(strand_levels)) > 0L
  if (any(!strand_levels[used] %in% c("+", "-", "*")))
    stop("OFST strand values must be '+', '-', '*' or missing.")
  for (col in intersect(c("cigar", "cigar1", "cigar2"), names(d))) {
    if (is.factor(d[[col]]))
      data.table::set(d, j = col, value = as.character(.ofst_factor(d[[col]], col)))
  }
  for (col in intersect(c("start", "start1", "start2", "end", "size", "width"), names(d))) {
    x <- d[[col]]
    if (!is.double(x)) next
    bounds <- suppressWarnings(range(x, na.rm = TRUE))
    if (bounds[1L] < -.int32_max || bounds[2L] > .int32_max ||
        any(x != trunc(x), na.rm = TRUE))
      stop("OFST coordinate column '", col, "' must contain 32-bit integer-valued coordinates.")
    data.table::set(d, j = col, value = as.integer(x))
  }
  d
}

.read_fst_list <- function(paths) {
  lapply(paths, function(x) tryCatch(
    .normalize_dt(fst::read_fst(x, as.data.table = TRUE)),
    error = function(e) stop("Invalid OFST file '", x, "': ", conditionMessage(e), call. = FALSE)))
}

.validate_schema <- function(paths) {
  schemas <- lapply(paths, function(x) tryCatch(fst::metadata_fst(x)$columnNames,
    error = function(e) .ofst_abort("cannot inspect OFST file '", x, "': ", conditionMessage(e))))
  for (i in seq_along(schemas)) {
    if (anyDuplicated(schemas[[i]])) .ofst_abort("duplicate column names in '", paths[i], "'.")
    if (!setequal(schemas[[1L]], schemas[[i]]))
      .ofst_abort("only ofst files with identical columns can be merged. File '", paths[i],
                  "' differs from '", paths[1L], "'; missing: ",
                  paste(setdiff(schemas[[1L]], schemas[[i]]), collapse = ", "),
                  "; extra: ", paste(setdiff(schemas[[i]], schemas[[1L]]), collapse = ", "), ".")
  }
  invisible(schemas[[1L]])
}

.plan_splits <- function(row_counts, limit = 2^31, max_splits = 20L) {
  if (!is.numeric(limit) || length(limit) != 1L || is.na(limit) ||
      !is.finite(limit) || limit <= 0 || limit > 2^31)
    stop("dt_max_index_size must be positive and at most 2^31.")
  if (!is.numeric(max_splits) || length(max_splits) != 1L || is.na(max_splits) ||
      !is.finite(max_splits) || max_splits < 1 || max_splits != trunc(max_splits) ||
      max_splits > .Machine$integer.max)
    stop("max_splits must be a positive integer.")
  total <- sum(row_counts)
  if (total < limit) {
    return(list(must_split = FALSE, is_single_chunked = TRUE, splits = 1L,
                split_vector = rep(1L, length(row_counts)),
                message = "Merging all libraries without splitting required.."))
  }
  # More groups than files cannot improve a file-based split, and a huge
  # max_splits must not itself allocate a huge planning vector.
  for (s in seq_len(min(max_splits, length(row_counts)))[-1L]) {
    tmp_vec <- ceiling(seq_len(length(row_counts)) / (length(row_counts) / s))
    group_sums <- vapply(split(row_counts, tmp_vec), sum, numeric(1))
    if (all(group_sums < limit)) {
      return(list(must_split = TRUE, is_single_chunked = FALSE,
                  splits = s, split_vector = tmp_vec, message = "Split round 1",
                  note_message = paste0("Total rows exceed limit (", limit,
                                        ") - splitting into ", s, " chunk(s).")))
    }
  }
  .ofst_abort("max_splits = ", max_splits, " is not enough for safe chunking at dt_max_index_size=", limit,
              ". Largest input has ", .ofst_number(max(row_counts)),
              " rows; this splitter groups whole files. Increase the split limit/count or use the scratch-backed rescue by lowering filter_target_rows.")
}

.merge_by_keys_reduce <- function(dt_list, by, sort = FALSE) {
  Reduce(function(x, y) merge.data.table(x, y, by = by, all = TRUE, sort = sort), dt_list)
}

.keys_nonlib_nonscore <- function(dt, lib_names) {
  setdiff(names(dt), c("score", lib_names))
}

.recompute_total_score <- function(dt, lib_names) {
  present_libs <- intersect(lib_names, names(dt))
  if (!length(present_libs)) stop("No per-library columns present; cannot recompute total score.")
  # Avoid rowSums(.SD)'s N-row by library-count dense matrix.
  total <- numeric(nrow(dt))
  for (column in present_libs) {
    x <- as.double(dt[[column]])
    if (anyNA(x)) x[is.na(x)] <- 0
    total <- total + x
  }
  data.table::set(dt, j = "score", value = total)
  dt
}

.int32_max <- 2^31 - 1
.int32_min <- -.int32_max

.is_safe_int32 <- function(x) {
  if (!is.double(x)) return(FALSE)                 # only handle double -> int
  if (!length(x)) return(TRUE)
  if (anyNA(x) || any(!is.finite(x))) return(FALSE)
  # exact integers within int32 bounds
  if (any(x != round(x))) return(FALSE)
  r <- range(x)
  (r[1] >= .int32_min) && (r[2] <= .int32_max)
}

.downcast_cols_if_safe <- function(dt, cols) {
  cols <- intersect(cols, names(dt))
  for (col in cols) {
    x <- dt[[col]]
    if (is.integer(x)) next
    if (.is_safe_int32(x)) data.table::set(dt, j = col, value = as.integer(round(x)))
  }
  invisible(dt)
}

.downcast_score_if_safe <- function(dt) .downcast_cols_if_safe(dt, "score")

.ofst_merge_lossless <- function(dt_list, lib_names, keep_all_scores = TRUE,
                                keepCigar = TRUE, sort = TRUE,
                                chunkified = FALSE) {
  if (!is.list(dt_list) || !length(dt_list) ||
      !all(vapply(dt_list, data.table::is.data.table, logical(1))))
    .ofst_abort("dt_list must contain at least one data.table, with no non-table elements.")
  if (!chunkified && length(dt_list) != length(lib_names))
    .ofst_abort("lib_names must have one name per input table.")

  colnames <- names(dt_list[[1L]])
  non_key_columns <- "score"
  if (chunkified && keep_all_scores) non_key_columns <- c(non_key_columns, lib_names)
  merge_keys <- setdiff(colnames, non_key_columns)
  if (!keepCigar) merge_keys <- setdiff(merge_keys, "cigar")

  if (keep_all_scores) {
    if (!chunkified) {
      # Normalize & collapse per-library BEFORE merging
      dt_list <- lapply(dt_list, function(d) {
        keys <- setdiff(names(d), "score")
        if (!keepCigar && "cigar" %in% keys) keys <- setdiff(keys, "cigar")
        d[, .(score = sum(score, na.rm = TRUE)), by = c(keys)]
      })

      # Rename each library's 'score' to its lib name
      for (i in seq_along(dt_list))
        data.table::setnames(dt_list[[i]], "score", lib_names[i])
    } else dt_list <- lapply(dt_list, function(d) {d[, score := NULL]})

    # Outer-merge libraries on keys, then recompute total
    dt <- .merge_by_keys_reduce(dt_list, by = merge_keys, sort = FALSE)
    # Downcast per-lib cols if safe, then total
    dt <- .downcast_cols_if_safe(dt, lib_names)
    dt <- .recompute_total_score(dt, lib_names)
  } else {
    if (!chunkified) {
      if (!keepCigar && "cigar" %in% colnames) {
        dt_list <- lapply(dt_list, function(d) { d[, cigar := NULL]; d })
      }
    }
    # A singleton still needs collapsing (it may contain duplicate keys), but
    # must not enter the reduction loop or return an uninitialized result.
    if (length(dt_list) == 1L) {
      dt <- dt_list[[1L]][, .(score = sum(score, na.rm = TRUE)), by = c(merge_keys)]
      dt_list <- NULL
    }
    while (length(dt_list) > 1) {
      nrows <- vapply(dt_list, nrow, numeric(1))
      sums <- cumsum(nrows)

      indices <- which(sums < .int32_max)
      message("-- Merging sub-chunks: ", paste(indices, collapse = ", ")," (of total: ", length(dt_list), ")")
      if (length(indices) < 2) {
        message("No pair of chunks has total rows smaller than '.int32_max' : ",
             .int32_max, ", can not merge using BOOST multithreading, using slow join.")

        score_names <- paste0(".ofst_score_", seq_along(dt_list))
        while (any(score_names %in% merge_keys)) score_names <- paste0("_", score_names)
        for (i in seq_along(dt_list)) {
          # Raw files may contain duplicate keys. Joining them directly can
          # multiply counts; second-round chunks are already aggregated.
          if (!chunkified)
            dt_list[[i]] <- dt_list[[i]][, .(score = sum(score, na.rm = TRUE)), by = c(merge_keys)]
          data.table::setnames(dt_list[[i]], "score", score_names[i])
        }
        # Outer-merge libraries on keys, then recompute total
        dt <- tryCatch(.merge_by_keys_reduce(dt_list, by = merge_keys, sort = FALSE),
          error = function(e) .ofst_abort("lossless chunk join failed. Retry ofst_merge on the original files with a smaller filter_target_rows to use global scratch-backed rescue. Cause: ", conditionMessage(e)))
        dt <- .downcast_cols_if_safe(dt, score_names)
        dt <- .recompute_total_score(dt, score_names)
        dt[, (score_names) := NULL]
        dt_list <- list()
      } else {
        # Drop both the previous accumulator alias and the consumed input
        # tables before grouping allocates its workspace and result.
        dt <- NULL
        combined <- data.table::rbindlist(dt_list[indices], use.names = TRUE, fill = TRUE)
        dt_list[indices] <- NULL
        gc()
        dt <- combined[, .(score = sum(score, na.rm = TRUE)), by = c(merge_keys)]
        combined <- NULL
        dt_list <- c(list(dt), dt_list)
      }
      gc()
    }
  }
  # Downcast total score if safe
  dt <- .downcast_score_if_safe(dt)
  if (sort) dt <- .ofst_sort(dt, setdiff(names(dt), "score"))
  dt
}

#' Collapse duplicated reads
#'
#' For every GRanges, GAlignments read, with the same:
#' seqname, start, (cigar) / width and strand, collapse and give a new
#' meta column called "score", which contains the number of duplicates
#' of that read. If score column already exists, will return input object!
#' @param x a GRanges, GAlignments or GAlignmentPairs object
#' @param addScoreColumn logical, default: (TRUE), if FALSE,
#' only collapse and not keep score column of counts for collapsed reads.
#' Returns directly without collapsing if reuse.score.column is FALSE and
#' score is already defined.
#' @param ... alternative arguments for class instances. For example, see:
#' \code{?'collapseDuplicatedReads,GRanges-method'}
#' @return a GRanges, GAlignments, GAlignmentPairs or data.table object,
#'  same as input
#' @export
#' @examples
#' gr <- rep(GRanges("chr1", 1:10,"+"), 2)
#' collapseDuplicatedReads(gr)
setGeneric("collapseDuplicatedReads", function(x, addScoreColumn = TRUE, ...) standardGeneric("collapseDuplicatedReads"))

#' @inherit collapseDuplicatedReads
#' @inheritParams convertToOneBasedRanges
setMethod("collapseDuplicatedReads", "GRanges",
          function(x, addScoreColumn = TRUE, addSizeColumn = FALSE,
                   reuse.score.column = TRUE) {
  if (addSizeColumn) {
    if (!("size" %in% colnames(mcols(x))))
      stop("addSizeColumn is TRUE, and no size column found!")
  }

  dt <- data.table(seqnames = as.character(seqnames(x)),
                   start = start(ranges(x)),
                   end = end(ranges(x)),
                   strand = as.character(strand(x)))

  if (reuse.score.column & ("score" %in% colnames(mcols(x)))) { # reuse
    dt[, score := mcols(x)$score]
    if (addSizeColumn) {
      dt[, size := mcols(x)$size]
      dt <- dt[, .(score = sum(score)), .(seqnames, start, end, strand, size)]
    } else {
      dt <- dt[, .(score = sum(score)), .(seqnames, start, end, strand)]
    }
  } else { # Do not reuse or "score" does not exist
    if (addSizeColumn) {
      dt[, size := mcols(x)$size]
      dt <- dt[, .(score = .N), .(seqnames, start, end, strand, size)]
    } else {
      dt <- dt[, .(score = .N), .(seqnames, start, end, strand)]
    }
  }
  if (!addScoreColumn) dt$score <- NULL
  # TODO change makeGRangesFromDataFrame to internal fast function
  return(makeGRangesFromDataFrame(dt, keep.extra.columns = TRUE, seqinfo = seqinfo(x)))
})

#' @inherit collapseDuplicatedReads
#' @inheritParams convertToOneBasedRanges
setMethod("collapseDuplicatedReads", "GAlignments",
          function(x, addScoreColumn = TRUE, reuse.score.column = TRUE) {
  if (("score" %in% colnames(mcols(x))) & !reuse.score.column) return(x)

  dt <- data.table(seqnames = factor(seqnames(x)),
                   start = start(ranges(x)),
                   cigar = cigar(x),
                   strand = factor(strand(x)), levels = c("+", "-", "*"))
  if (reuse.score.column & ("score" %in% colnames(mcols(x)))) { # reuse
    dt[, score := mcols(x)$score]
    dt <- dt[, .(score = sum(score)), .(seqnames, start, cigar, strand)]
  } else { # Do not reuse or "score" does not exist
    dt <- dt[, .(score = .N), .(seqnames, start, cigar, strand)]
  }

  if (!addScoreColumn) dt$score <- NULL
  return(getGAlignments(dt))
})

#' @inherit collapseDuplicatedReads
setMethod("collapseDuplicatedReads", "GAlignmentPairs",
          function(x, addScoreColumn = TRUE) {
  if ("score" %in% colnames(mcols(x))) return(x)

  dt <- data.table(seqnames = factor(x@first@seqnames),
                   start1 = x@first@start,
                   start2 = x@last@start,
                   cigar1 = factor(x@first@cigar),
                   cigar2 = factor(x@last@cigar),
                   strand = factor(x@first@strand, levels = c("+", "-", "*")))
  dt <- dt[, .(score = .N), .(seqnames, start1, start2,
                              cigar1, cigar2, strand)]
  if (!addScoreColumn) dt$score <- NULL
  return(getGAlignmentsPairs(dt))
})

#' @inherit collapseDuplicatedReads
#' @param keepCigar logical, default FALSE. Keep the cigar information
#' @inheritParams convertToOneBasedRanges
setMethod("collapseDuplicatedReads", "data.table",
          function(x, addScoreColumn = TRUE, addSizeColumn = FALSE,
                   reuse.score.column = TRUE, keepCigar = FALSE) {
  is_GAlignmentPair <- "cigar1" %in% colnames(x)
  if (is_GAlignmentPair)
    stop("Paired end collapse on data.table not supported yet!")
  required_columns <- c("seqnames", "start", "strand", "score")
  stopifnot(all(required_columns %in% colnames(x)))
  size_exists <- "size" %in% colnames(x)
  cigar_exists <- "cigar" %in% colnames(x)
  grouping <- c("seqnames", "start", "strand")

  if (addSizeColumn) {
    if (size_exists) {
      grouping <- c(grouping, "size")
    } else
      warning("addSizeColumn is TRUE, and no size column found!")
  }

  if (keepCigar & cigar_exists) grouping <- c(grouping, "cigar")

  if (reuse.score.column & ("score" %in% colnames(x))) { # reuse
    x <- x[, .(score = sum(score)), by = grouping]
  } else { # Do not reuse or "score" does not exist
    x <- x[, .(score = .N), by = grouping]
  }
  if (!addScoreColumn) x$score <- NULL
  # TODO change makeGRangesFromDataFrame to internal fast function
  return(x)
})

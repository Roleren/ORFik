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
#' to merge split on lower counts
#' @return a data.table of merged result, it is merged on all columns except "score".
#' The returned file will contain the scores of each file + the aggregate sum score.
#' @importFrom data.table setnames
ofst_merge <- function(file_paths,
                       lib_names = sub("\\.ofst$", "", basename(file_paths)),
                       keep_all_scores = TRUE, keepCigar = TRUE, sort = TRUE,
                       max_splits = 20L, dt_max_index_size = 2^31) {
  stopifnot(length(file_paths) > 0 && is.character(file_paths))
  stopifnot(all(file.exists(file_paths)))
  stopifnot(length(file_paths) == length(lib_names))
  if (keep_all_scores) stopifnot(length(lib_names) == length(unique(lib_names)))
  .validate_schema(file_paths)

  # Plan splits vs data.table's index limit
  row_numbers <- unlist(lapply(file_paths, function(x)
    data.table::data.table(fst::metadata_fst(x)$nrOfRows)))
  plan <- .plan_splits(row_numbers, limit = dt_max_index_size,
                       max_splits = max_splits)

  # First-round: per-chunk merge
  file_paths_split <- split(file_paths, plan$split_vector)
  lib_names_split  <- split(lib_names,  plan$split_vector)

  is_single_chunked <- plan$is_single_chunked
  if (!is_single_chunked) message(plan$note_message)
  message(plan$message)
  sort_on_chunking <- is_single_chunked & sort
  chunk_list <- lapply(seq_along(file_paths_split), function(g) {
    message("- Merging chunk ", g, "/", length(file_paths_split))
    ofst_merge_internal(
      .read_fst_list(file_paths_split[[g]]),
      lib_names = lib_names_split[[g]],
      keep_all_scores = keep_all_scores,
      keepCigar = keepCigar,
      sort = sort_on_chunking
    )
  })

  # Second-round: combine chunks
  dt <- if (is_single_chunked) {
    chunk_list[[1L]]
  } else {
    message("Split round 2")
    ofst_merge_internal(
      chunk_list,
      lib_names = lib_names,
      keep_all_scores = keep_all_scores,
      keepCigar = keepCigar,
      sort = sort, chunkified = TRUE
    )
  }
  # Make seqnames and strand factor
  dt[, seqnames := as.factor(seqnames)]
  dt[, strand := factor(strand, levels = c("+", "-", "*"))]
  message("Done merging")
  return(dt[])
}

.normalize_dt <- function(d) {
  # Robustly coerce possible factors (even malformed) to character
  for (col in intersect(c("seqnames","strand","cigar"), names(d))) {
    v <- d[[col]]
    if (is.factor(v)) {
      # Decode factor codes via levels[] indexing to survive malformed factors
      lv <- levels(v)
      codes <- as.integer(v)
      d[[col]] <- as.character(lv[pmax.int(pmin.int(codes, length(lv)), 1L)])
      # Any out-of-range codes turn into NA via indexing outside [1, nlevels]
      d[[col]][is.na(codes) | codes < 1L | codes > length(lv)] <- NA_character_
    }
    if (!is.character(d[[col]])) d[[col]] <- as.character(d[[col]])
  }
  # Ensure integer-like coords are integer (optional, keeps joins precise)
  for (col in intersect(c("start","end","size","width"), names(d))) {
    if (is.double(d[[col]])) d[[col]] <- as.integer(round(d[[col]]))
  }
  return(d)
}

.read_fst_list <- function(paths) {
  lapply(paths, function(x) .normalize_dt(fst::read_fst(x, as.data.table = TRUE)))
}

.validate_schema <- function(paths) {
  meta <- table(unlist(lapply(paths, function(x)
    data.table::data.table(fst::metadata_fst(x)$columnNames))))
  if (!all(meta == length(paths))) {
    print(meta)
    stop("Some libraries had columns not in others! ",
         "Only ofst files with identical columns can be merged.")
  }
  invisible(TRUE)
}

.plan_splits <- function(row_counts, limit = 2^31, max_splits = 20L) {
  total <- sum(row_counts)
  if (total < limit) {
    return(list(must_split = FALSE, is_single_chunked = TRUE, splits = 1L,
                split_vector = rep(1L, length(row_counts)),
                message = "Merging all libraries without splitting required.."))
  }
  for (s in seq_len(max_splits)[-1L]) {
    tmp_vec <- ceiling(seq_len(length(row_counts)) / (length(row_counts) / s))
    group_sums <- vapply(split(row_counts, tmp_vec), sum, numeric(1))
    if (all(group_sums < limit)) {
      return(list(must_split = TRUE, is_single_chunked = FALSE,
                  splits = s, split_vector = tmp_vec, message = "Split round 1",
                  note_message = paste0("Total rows exceed limit (", limit,
                                        ") - splitting into ", s, " chunk(s).")))
    }
  }
  stop("max_splits = ", max_splits, " is not enough for safe chunking.")
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
  dt[, score := rowSums(.SD, na.rm = TRUE), .SDcols = present_libs]
  dt
}

.int32_max <- 2^31 - 1
.int32_min <- -.int32_max

.is_safe_int32 <- function(x) {
  if (!is.double(x)) return(FALSE)                 # only handle double -> int
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

ofst_merge_internal <- function(dt_list, lib_names, keep_all_scores = TRUE,
                                keepCigar = TRUE, sort = TRUE,
                                chunkified = FALSE) {
  if (!chunkified) stopifnot(length(dt_list) == length(lib_names))
  colnames <- names(dt_list[[1L]])
  non_key_columns <- "score"
  if (chunkified) non_key_columns <- c(non_key_columns, lib_names)

  if (keep_all_scores) {

    if (!chunkified) {
      # Normalize & collapse per-library BEFORE merging
      dt_list <- lapply(dt_list, function(d) {
        keys <- setdiff(names(d), "score")
        if (!keepCigar && "cigar" %in% keys) keys <- setdiff(keys, "cigar")
        d[, .(score = sum(score, na.rm = TRUE)), by = keys]
      })

      # Rename each library's 'score' to its lib name
      for (i in seq_along(dt_list))
        data.table::setnames(dt_list[[i]], "score", lib_names[i])
    } else dt_list <- lapply(dt_list, function(d) {d[, score := NULL]})

    merge_keys <- setdiff(colnames, non_key_columns)

    # 4) Outer-merge libraries on keys, then recompute total
    dt <- .merge_by_keys_reduce(dt_list, by = merge_keys, sort = FALSE)
    # Optional: downcast per-lib cols if safe, then total
    dt <- .downcast_cols_if_safe(dt, lib_names)
    dt <- .recompute_total_score(dt, lib_names)
  } else {
    if (!chunkified) {
      if (!keepCigar && "cigar" %in% colnames) {
        dt_list <- lapply(dt_list, function(d) { d[, cigar := NULL]; d })
      }
    }
    dt <- collapseDuplicatedReads(
      data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE),
      addSizeColumn = TRUE,
      keepCigar = keepCigar
    )
  }
  # Downcast total score if safe
  dt <- .downcast_score_if_safe(dt)
  if (sort) data.table::setorderv(dt, setdiff(names(dt), "score"))
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
                   strand = factor(strand(x)))
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
                   strand = factor(x@first@strand))
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

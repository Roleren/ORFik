#' Converts bed style data.frame to GRanges
#'
#' For info on columns, see:
#' https://www.ensembl.org/info/website/upload/bed.html
#' @param x A \code{\link{data.frame}} from imported bed-file,
#'  to convert to GRanges
#' @param skip.name default (TRUE), skip name column (column 4)
#' @return a \code{\link{GRanges}} object from bed
#' @family utils
#' @keywords internal
bedToGR <- function(x, skip.name = TRUE) {
  if (ncol(x) < 3) stop("bed file must have a least 3 columns!")
  starts <- x[, 2L] + 1L
  ends <- x[, 3L]
  gr <- GRanges(x[, 1L], IRanges(starts, ends))

  bed4 <- (ncol(x) >= 4) & !skip.name #name column
  if (bed4) mcols(gr) <- S4Vectors::DataFrame(mcols(gr), name = x[, 4L])
  bed5 <- ncol(x) >= 5 # score column
  if (bed5) mcols(gr) <- S4Vectors::DataFrame(mcols(gr), score = x[, 5L])
  bed6 <- ncol(x) >= 6
  if (bed6) strand(gr) <- x[, 6L]
  if (ncol(x) > 6L) mcols(gr) <- x[, 7L:ncol(x)]
  return(gr)
}

#' Internal GRanges loader from fst data.frame
#' @param df a data.frame/data.table with columns minimum 4 columns:
#' seqnames, start, strand\cr
#' Additional specific columns are:\cr
#' - width (if not set, width is set to 1 for all reads)\cr
#' Additional columns will be assigned as meta columns
#' @param keep.extra.columns logical, default TRUE, keep meta cols.
#' @inheritParams import.ofst
#' @return GRanges object
#' @importFrom S4Vectors new2
#' @keywords internal
getGRanges <- function(df, keep.extra.columns = TRUE, seqinfo = NULL) {
  stopifnot(is(df, "data.frame"))
  if (!isTRUEorFALSE(keep.extra.columns))
    stop("'keep.extra.columns' must be TRUE or FALSE")
  if (!all(c("seqnames", "start", "strand") %in% colnames(df)))
    stop("df must at minimum have 4 columns named: seqnames, start, width and strand")
  if (!is(df, "data.table")) setDT(df)

  widths <- if ("width" %in% colnames(df)) {
    df$width
  } else {
    if ("end" %in% colnames(df)) {
      df$end - df$start  + 1L
    } else rep.int(1L, nrow(df))
  }

  ranges <- new2("IRanges", start = df$start,
                 width = widths,
                 NAMES = df$NAMES,
                 elementMetadata = NULL,
                 check = FALSE)
  if (is.null(levels(df$seqnames))) {
    df[, seqnames := factor(seqnames, levels = unique(seqnames))]
  }
  if (is.null(levels(df$strand))) {
    df[, strand := factor(strand, levels = c("+", "-", "*"))]
  }
  if (!is.null(seqinfo)) {
    stopifnot(is(seqinfo, "Seqinfo"))
  } else seqinfo <- Seqinfo(levels(df$seqnames))

  if (!is.null(df$NAMES)) df$NAMES <- NULL


  mcols <- NULL
  if (keep.extra.columns) {
    mcols_hits <- !(colnames(df) %in% c("seqnames", "start", "strand", "width", "end"))
    has_mcols <- any(mcols_hits)
    if (has_mcols) {
      mcols <- df[, which(mcols_hits), with = FALSE]
    }
  }
  mcols <- S4Vectors:::normarg_mcols(mcols, "GRanges", nrow(df))

  new2("GRanges", seqnames = Rle(df$seqnames), ranges = ranges, strand = Rle(df$strand),
       elementMetadata = mcols, seqinfo = seqinfo, check = FALSE)
}

#' Faster version (also less safe) of makeGRangesFromDataFrame
#'
#' @inheritParams getGRanges
#' @return GRanges object
#' @export
#' @examples
#' df <- data.frame(start = rep(1L, 1e5), end = 10L, strand = "+", seqnames = "1")
#'
#' system.time(res <- makeGRangesFromDataFrame(df))
#' system.time(res_fast <- makeGRangesFromDataFrameFast(df))
#' identical(res, res_fast)
#'
#' # Use width instead of end, does not work in original
#' df2 <- data.frame(start = rep(1L, 1e5), width = 10L, strand = "+", seqnames = "1")
#' system.time(makeGRangesFromDataFrameFast(df2))
#' df_small <- data.frame(start = 1L, end = 10L, strand = "+", seqnames = "1")
#' system.time(res <- makeGRangesFromDataFrame(df_small))
#' system.time(res_fast <- makeGRangesFromDataFrameFast(df_small))
#' identical(res, res_fast)
makeGRangesFromDataFrameFast <- function(df, keep.extra.columns = TRUE,
                                         seqinfo = NULL) {
  return(getGRanges(df, keep.extra.columns, seqinfo))
}

#' Internal GAlignments loader from fst data.frame
#' @param df a data.frame/data.table with columns minimum 4 columns:
#' seqnames, start ("pos" in final GA object), cigar and strand.\cr
#' Additional columns will be assigned as meta columns
#' @inheritParams import.ofst
#' @return GAlignments object
#' @importFrom S4Vectors new2
#' @keywords internal
getGAlignments <- function(df, seqinfo = NULL) {
  stopifnot(is(df, "data.frame"))
  if (!all(c("seqnames", "start", "cigar", "strand") %in% colnames(df)))
    stop("df must at minimum have 4 columns named: seqnames, start, cigar and strand")
  if (nrow(df) == 0) return(GenomicAlignments::GAlignments())
  if (!is(df, "data.table")) setDT(df)

  if (is.null(levels(df$seqnames))) {
    df[, seqnames := factor(df$seqnames, levels = unique(df$seqnames))]
  }
  if (!is.null(seqinfo)) {
    stopifnot(is(seqinfo, "Seqinfo"))
  } else seqinfo <- Seqinfo(levels(df$seqnames))

  names <- df$NAMES
  if (!is.null(df$NAMES)) df$NAMES <- NULL
  if (ncol(df) == 4){
    mcols <- NULL
  } else {
    mcols <- df[,5:ncol(df)]
    if (ncol(df) == 5) { # Hm... Is this safe ? What if a score is there ?
      mcols <- data.frame(mcols)
      names(mcols) <- names(df)[5]
    }
  }
  if (!is(df$strand, "factor") && identical(levels(df$strand), c("+", "-", "*"))){
    df[, strand := factor(strand, levels = c("+", "-", "*"))]
  }

  mcols <- S4Vectors:::normarg_mcols(mcols, "GRanges", nrow(df))
  new2("GAlignments", NAMES = names, seqnames = Rle(df$seqnames), start = df$start,
       cigar = as.character(df$cigar), strand = Rle(df$strand), elementMetadata = mcols,
       seqinfo = seqinfo, check = FALSE)

}

#' Internal GAlignmentPairs loader from fst data.frame
#' @param df a data.frame with columns minimum 6 columns:
#' seqnames, start1/start2 (integers), cigar1/cigar2 and strand\cr
#' Additional columns will be assigned as meta columns
#' @inheritParams import.ofst
#' @return GAlignmentPairs object
#' @importFrom S4Vectors new2
#' @keywords internal
getGAlignmentsPairs <- function(df, strandMode = 0, seqinfo = NULL) {
  if (nrow(df) == 0) {
    return(GenomicAlignments::GAlignmentPairs(first = GAlignments(),
                                              last = GAlignments(),
                                              strandMode = strandMode))
  }
  if (is.null(levels(df$seqnames))) {
    df$seqnames <- factor(df$seqnames, levels = unique(df$seqnames))
  }
  if (!is.null(seqinfo)) {
    stopifnot(is(seqinfo, "Seqinfo"))
  } else seqinfo <- Seqinfo(levels(df$seqnames))
  names <- df$NAMES
  if (!is.null(df$NAMES)) df$NAMES <- NULL
  if (ncol(df) == 6){
    mcols <- NULL
  } else {
    mcols <- df[,7:ncol(df)]
    if (ncol(df) == 7) { # Hm... Is this safe ? What if a score is there ?
      mcols <- data.frame(mcols)
      names(mcols) <- names(df)[7]
    }
  }
  mcols <- S4Vectors:::normarg_mcols(mcols, "GRanges", nrow(df))
  # reverse strand for last
  levels = c("+", "-", "*")
  strand2 <- strandTemp <- df$strand <-
    factor(df$strand, levels = levels)
  strandTemp[strand2 == "+"] <- "-"
  strandTemp[strand2 == "-"] <- "+"
  strand2 <- strandTemp
  new2("GAlignmentPairs",
       first = new2("GAlignments", NAMES = names, seqnames = Rle(df$seqnames), start = df$start1,
        cigar = as.character(df$cigar1), strand = Rle(df$strand),
        seqinfo = seqinfo, check = FALSE,
        elementMetadata = DataFrame(data.frame(matrix(nrow = nrow(df), ncol = 0)))),
       last = new2("GAlignments", NAMES = names, seqnames = Rle(df$seqnames), start = df$start2,
             cigar = as.character(df$cigar2), strand = Rle(factor(strand2, levels)),
             seqinfo = seqinfo, check = FALSE,
             elementMetadata = DataFrame(data.frame(matrix(nrow = nrow(df), ncol = 0)))),
       isProperPair = rep(TRUE, nrow(df)), strandMode = as.integer(strandMode),
       elementMetadata = mcols, check = FALSE)
}

#' Find pair of forward and reverse strand wig / bed files and
#' paired end bam files split in two
#'
#' @param paths a character path at least one .wig / .bed file
#' @param f Default (c("forward", "fwd")
#' a character vector for forward direction regex.
#' @param r Default (c("reverse", "rev")
#' a character vector for reverse direction regex.
#' @param format default "wig", for bed do "bed". Also searches
#' compressions of these variants.
#' @return if not all are paired, return original list,
#' if they are all paired, return a data.table with matches as 2 columns
#' @keywords internal
findNGSPairs <- function(paths, f = c("forward", "fwd"),
                         r = c("reverse", "rev"), format = "wig") {
  if (length(paths) == 0) return(paths)

  f <- paste0(f, "\\.", format, "*", collapse = "|")
  r <- paste0(r, "\\.", format, "*", collapse = "|")
  forwardPath <- grep(f, paths)
  reversePath <- grep(r, paths)

  if ((length(forwardPath) != length(reversePath)) |
      length(forwardPath) == 0 | length(reversePath) == 0) return(paths)
  dt <- data.table(forward = paths[forwardPath],
                   reverse = paths[reversePath], match = FALSE)
  dt.sub <- dt[, .(forward = gsub(pattern = f, x = forward, ""),
                   reverse = gsub(pattern = r, x = reverse, ""))]
  matches <- chmatch(dt.sub$forward, dt.sub$reverse)
  if (!anyNA(matches)) {
    dt <- dt[, .(forward, reverse = reverse[matches], match = TRUE)]
  }
  if (all(dt$match)) {
    return(dt)
  }
  return(paths)
}

#' Remove file extension of path
#'
#' Allows removal of compression
#' @param path character path (allows multiple paths)
#' @param basename relative path (TRUE) or full path (FALSE)? (default: FALSE)
#' @importFrom tools file_ext
#' @return character path without file extension
#' @keywords internal
remove.file_ext <- function(path, basename = FALSE) {
  fext <- file_ext(path)
  compressions <- c("gzip", "gz", "bgz", "zip")
  areCompressed <- fext %in% compressions
  regex <- paste0("*\\.",fext,"$")
  if (any(areCompressed)) {
    ext <- file_ext(file_path_sans_ext(path[areCompressed], compression = FALSE))
    regex[areCompressed] <- paste0("*\\.",ext,"\\.", fext[areCompressed],"$")
  }
  out <- c()
  for (p in seq_along(path)) {
    new <- gsub(pattern = regex[p], "", path[p])
    out <- c(out, new)
  }
  return(if(basename){basename(out)}else out)
}

#' A wrapper for seqlevelsStyle
#'
#' To make sure chromosome naming is correct (chr1 vs 1 vs I etc)
#' @param range a ranged object, (GRanges, GAlignment etc)
#' @param chrStyle a GRanges object, TxDb, FaFile,
#' , a \code{\link{seqlevelsStyle}} or \code{\link{Seqinfo}}.
#' (Default: NULL) to get seqlevelsStyle from. In addition if it
#' is a Seqinfo object, seqinfo will be updated.
#' Example of seqlevelsStyle update:
#' Is chromosome 1 called chr1 or 1,
#'  is mitocondrial chromosome called MT or chrM etc.
#' Will use 1st seqlevel-style if more are present.
#' Like: c("NCBI", "UCSC") -> pick "NCBI"
#' @return a GAlignment/GRanges object depending on input.
#' @keywords internal
matchSeqStyle <- function(range, chrStyle = NULL) {
  # TODO: Check if this makes it safer:
  # if (tryCatch(seqlevelsStyle(cage) <- seqlevelsStyle(fiveUTRs),
  #              error = function(e) {TRUE}) == TRUE) {
  #   warning("seqlevels of CAGE/fiveUTRs are not standardized, check them.")
  # } else {
  #   seqlevelsStyle(cage) <- seqlevelsStyle(fiveUTRs)
  # }
  if (!is.null(chrStyle)) {
    valid_seq_style <- seqlevelsStyleSafe(range)
    if (!is.null(valid_seq_style)) {
      if (is.character(chrStyle)) {
        if (!any(seqlevelsStyle(range) %in% chrStyle)) {
          try(seqlevelsStyle(range) <- chrStyle[1], silent = TRUE)
        }
      } else if (is.gr_or_grl(chrStyle) | is(chrStyle, "TxDb") |
                 is(chrStyle, "FaFile") | is(chrStyle, "Seqinfo")) {
        style_is_different <- !any(seqlevelsStyle(range) %in% seqlevelsStyle(chrStyle)[1])
        if (style_is_different) {
          if (!any(seqlevelsStyle(range) %in% chrStyle)) {
            try(seqlevelsStyle(range) <- seqlevelsStyle(chrStyle)[1], silent = TRUE)
          }
        }
      } else stop("chrStyle must be valid GRanges object,",
                  "or a valid chr style!")
    }
  }
  if (is(chrStyle, "Seqinfo")) {
    seqlevels(range) <- seqlevels(chrStyle)
    seqinfo(range) <- chrStyle
  }
  return(range)
}

seqlevelsStyleTest <- function(x) {
  try(seqlevelsStyle(x), silent = TRUE)
}

seqlevelsStyleSafe <- function(x) {
  res <- seqlevelsStyleTest(x)
  if (!is(res, "try-error")) return(res)
  return(NULL)
}

#' Find optimized subset of valid reads
#'
#' Keep only the ones that overlap within the grl ranges.
#' Also sort them in the end
#' @inheritParams validSeqlevels
#' @return the reads as GRanges,  GAlignment or GAlignmentPairs
#' @importFrom GenomicAlignments first
#' @family utils
#' @keywords internal
optimizeReads <- function(grl, reads) {
  if (is(reads, "covRle")) return(reads)
  seqMatch <- validSeqlevels(grl, reads)
  reads <- keepSeqlevels(reads, seqMatch, pruning.mode = "coarse")

  reads <- reads[countOverlaps(reads, grl, type = "within") > 0]
  reads <- if (is(reads, "GAlignmentPairs")) {
    message("For GAlignmentPairs using 'first' read of pair only,",
            " Make sure it is correct!")
    reads <- reads[order(GenomicAlignments::first(reads))]
    } else reads <- sort(reads)

  if (length(reads) == 0) warning("No reads left in 'reads' after",
                                    "optimisation!")

  return(reads)
}

#' Transform object
#'
#' Similar to normal transform like log2 or log10.
#' But keep 0 values as 0, to avoid Inf values and negtive values
#' are made as -scale(abs(x)), to avoid NaN values.
#' @param x a numeric vector or data.frame/data.table of numeric columns
#' @param scale a function, default log2, which function to transform with.
#' @param by.reference logical, FALSE. if TRUE, update object by reference
#' if it is data.table.
#' @return same object class as x, with transformed values
#' @keywords internal
pseudo.transform <- function(x, scale = log2, by.reference = FALSE) {

  if (!by.reference)
    x <- copy(x)
  xl <- x
  if (is.numeric(x)) {
    x[(x > 0)] <- scale(x[(x > 0)])
    x[(xl < 0)] <- - scale(abs(x[(xl < 0)]))
  } else if (is.data.frame(x)) {
    for (i in names(x)) {
      b <- x[, get(i)]
      bigger <- b > 0
      smaller <- (b < 0)
      b[bigger] <- scale(b[bigger])
      b[smaller] <- -scale(b[abs(smaller)])
      x[, (i):= b]
    }
  } else stop("x must be numeric or data.frame/data.table/matrix")

  return(x)
}

#' Create all unique combinations pairs possible
#'
#' Given a character vector, get all unique combinations of 2.
#' @param x a character vector, will unique elements for you.
#' @return a list of character vector pairs
#' @importFrom utils combn
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' ORFik:::combn.pairs(df[, "libtype"])
combn.pairs <- function(x) {
  pairs <- list() # creating compairisons :list of pairs
  comparisons.design <- unique(x)
  my_comparison <- combn(unique(comparisons.design), 2)
  pairs <- list()
  for (i in seq(ncol(my_comparison))) {
    pairs[[i]] <- c(my_comparison[1, i], my_comparison[2, i])
  }
  return(pairs)
}

#' Read RDS or QS format file
#'
#' @param file path to file with "rds" or "qs" file extension
#' @param nthread numeric, number of threads for qs::qread
#' @return R object loaded from file
#' @importFrom qs qread
#' @export
#' @examples
#' df <- ORFik::ORFik.template.experiment()
#' path <- ORFik:::countTablePath(df)
#' read_RDSQS(path)
read_RDSQS <- function(file, nthread = 5) {
  format <- file_ext(file)
  stopifnot(format %in% c("qs", "rds", "covqs", "covrds"))
  if (format %in% c("rds", "covrds")) {
    readRDS(file)
  } else qs::qread(file, nthread = nthread)
}

#' Read RDS or QS format file
#'
#' @param object the object to save
#' @param file path to file with "rds" or "qs" file extension
#' @param nthread numeric, number of threads for qs::qread
#' @return R object loaded from file
#' @importFrom qs qread
#' @export
#' @examples
#' path <- tempfile(fileext = ".qs")
#' # Simple numeric save
#' x <- 1
#' save_RDSQS(x, path)
#' read_RDSQS(path)
#' # Save a list
#' x <- list(a = 1, b = c(1,2,3))
#' save_RDSQS(x, path)
#' read_RDSQS(path)
save_RDSQS <- function(object, file, nthread = 5) {
  stopifnot(is(file, "character"))
  format <- file_ext(file)
  stopifnot(format %in% c("qs", "rds", "covqs", "covrds"))
  if (format %in% c("rds", "covrds")) {
    saveRDS(object, file)
  } else qs::qsave(object, file, nthread = nthread)
}

#' A fast ftp directory check
#'
#' Check if ftp directory exists
#' @param url.dir character, url to a ftp directory.
#' @param report.error logical, FALSE. If TRUE,
#' stop and report error.
#' @return logical, TRUE if url directory exists
#' @importFrom httr http_error
#' @keywords internal
exists.ftp.dir.fast <- function(url.dir, report.error = FALSE) {
  url.dir.not.exists <- tryCatch(
    expr = {
      httr::http_error(url.dir)
    },
    error = function(e){
      e
    })
  if (is.logical(url.dir.not.exists)) { # Directory exists
    return(TRUE)
  } else { # Directory does not exist
    if (report.error) {
      msg <- url.dir.not.exists$message
      if (length(grep("Timeout was reached.", msg)) == 1 |
          length(grep("Recv failure:", msg)) == 1) {
        stop("FTP port 21 is not connecting to || ", url.dir,
             " || wait 5 minutes, if it still does not work",
             " check your firewall and FTP settings!")
      } else if (length(grep("Could not resolve host", msg)) == 1){
        stop("Server is not up:,", url.dir, "check again in 5 minutes!")
      }  else if (length(grep("Server denied you to change", msg)) == 1){
        stop("Directory does not exist: ", url.dir)
      }
    }
  }
  # Else return that file does not exist
  return(FALSE)
}

#' A fast ftp file check
#'
#' Check if ftp file exists
#' @param url character, url to a ftp file
#' @param report.error logical, FALSE. If TRUE,
#' stop and report error.
#' @return logical, TRUE if file exists
#' @importFrom httr http_error
#' @keywords internal
exists.ftp.file.fast <- function(url, report.error = FALSE) {
  # Stop early if directory does not exist
  if (!exists.ftp.dir.fast(paste0(dirname(url), "/"), report.error))
    return(FALSE)
  # Else check if file exists
  url.not.exists <- tryCatch(
    expr = {
      httr::http_error(url)
    },
    error = function(e){
      e
    })

  if (is.logical(url.not.exists)) { # File exists
    return(TRUE)
  } else { # File does not exist
    if (report.error) {
      msg <- url.not.exists$message
      if (length(grep("Timeout was reached.", msg)) == 1 |
          length(grep("Recv failure:", msg)) == 1) {
        stop("FTP port 21 is not connecting to || ", url,
             " || wait 5 minutes, if it still does not work",
             " check your internet connection, firewall and FTP settings!")
      } else if (length(grep("Given file does not exist", msg)) == 1){
        stop("File does not exist,", url)
      }
    }
  }
  return(FALSE)
}

local_DTthreads <- function(n = 1, env = parent.frame()) {
  old <- data.table::getDTthreads()
  withr::defer(data.table::setDTthreads(old), envir = env)
  data.table::setDTthreads(n)
}

#' Detects the mounted drive based on a mounted path
#' @param ref_path = path.expand(config()["ref"])
#' @return character, name of FileSystem drive of mounted path,
#' NA_character_ if not found
#' @export
#' @examples
#' detect_drive(tempdir())
detect_drive <- function(ref_path = path.expand(config()["ref"])) {
  if (.Platform$OS.type != "unix") return(NA_character_)
  ref_path <- path.expand(ref_path)
  # Get disk usage information
  drive_info <- system("df -h", intern = TRUE)[-1]

  # Iterate through the lines (skip the first row, which is the header)
  candidate_drives <- c()
  for (line in drive_info) {
    parts <- strsplit(line, " +")[[1]]  # Split by spaces

    # Ensure we have enough columns (Filesystem, Size, Used, Avail, Use%, Mounted on)
    if (length(parts) >= 6 && grepl(paste0("^", parts[6]), ref_path)) {
      candidate_drives <- c(candidate_drives, parts[1])  # Store filesystem path
    }
  }

  # Filter out root "/" if another drive is available
  if (length(candidate_drives) > 1) {
    candidate_drives <- setdiff(candidate_drives, system("df -h / | awk 'NR==2{print $1}'", intern = TRUE))
  }

  if (length(candidate_drives) != 1) {
    warning("Could not determine the correct unique Filesystem drive for the reference path:\n  ",
            ref_path)
    return(NA_character_)
  }

  return(candidate_drives)
}


#' System usage for Linux (Auto-detects correct drive if not provided)
#' @param drive path, the Filesystem drive !(Not the mounted name),
#' use drive = ORFik:::detect_drive("My_directory_inside this mount_name") to
#' get custom drive
#' @param one_liner Logical, default FALSE. Instead return a length 1 character string
#' with all the info.
#' @return A list with system info, if one_liner is TRUE, then a length 1 character string.
#' @export
#' @examples
#' get_system_usage()
get_system_usage <- function(drive = detect_drive(), one_liner = FALSE) {
  is_windows <- .Platform$OS.type != "unix"
  if (is_windows) return(list())

  sysname <- Sys.info()[["sysname"]]
  is_linux <- sysname == "Linux"
  is_macos <- sysname == "Darwin"

  # ---- CPU usage ----
  cpu_call <- if (is_linux) {
    "top -bn1 | grep 'Cpu(s)' | awk '{print $2 + $4}'"
  } else if (is_macos) {
    "top -l 1 | grep 'CPU usage' | awk -F'[:,]' '{print $2}' | awk -F'%' '{print 100 - $1}'"
  } else {
    return(list())
  }
  cpu_usage <- as.numeric(system(cpu_call, intern = TRUE))

  # ---- Memory usage (in GB) ----
  if (is_linux) {
    mem_info <- system("free -g | awk 'NR==2{print $3, $2}'", intern = TRUE)
    mem_vals <- as.numeric(strsplit(trimws(mem_info), "\\s+")[[1]])
    mem_usage <- mem_vals[1]
    mem_total <- mem_vals[2]
  } else if (is_macos) {
    page_size <- as.numeric(system("sysctl -n hw.pagesize", intern = TRUE))
    vm_stats <- system("vm_stat", intern = TRUE)
    get_pages <- function(label) {
      as.numeric(gsub("[^0-9]", "", grep(label, vm_stats, value = TRUE)))
    }
    free_pages <- get_pages("Pages free:")
    active_pages <- get_pages("Pages active:")
    speculative_pages <- get_pages("Pages speculative:")
    inactive_pages <- get_pages("Pages inactive:")
    wired_pages <- get_pages("Pages wired down:")

    used_pages <- active_pages + inactive_pages + wired_pages + speculative_pages
    total_pages <- used_pages + free_pages

    mem_usage <- round((used_pages * page_size) / 1e9, 2)  # in GB
    mem_total <- round((total_pages * page_size) / 1e9, 2)
  } else {
    return(list())
  }

  mem_percent <- round((mem_usage / mem_total) * 100, 2)

  # ---- Drive usage ----
  if (!is.na(drive)) {
    drive_line <- suppressWarnings(system(paste0("df -h | grep '",
                              drive, "'", " | tail -1"), intern = TRUE))
    if (length(attr(drive_line, "status")) == 1 || length(drive_line) == 0) {
      drive_vals <- as.character(rep(NA, 6))
    } else {
      drive_vals <- strsplit(trimws(drive_line), " +")[[1]]
    }
  } else {
    drive_vals <- as.character(rep(NA, 6))
  }

  drive_total <- drive_vals[2]
  drive_used <- drive_vals[3]
  drive_free <- drive_vals[4]
  drive_percent <- drive_vals[5]

  # ---- Output ----
  usage <- list(
    CPU_Usage_Percent = cpu_usage,
    Memory_Usage_GB = mem_usage,
    Memory_Total_GB = mem_total,
    Memory_Usage_Percent = mem_percent,
    Drive = drive,
    Drive_Total = drive_total,
    Drive_Used = drive_used,
    Drive_Free = drive_free,
    Drive_Usage_Percent = drive_percent
  )

  if (one_liner) {
    usage <- get_system_usage_one_liner(usage)
  }

  return(usage)
}


get_system_usage_one_liner <- function(usage) {
  cat(paste0("CPU (", usage$CPU_Usage_Percent, "%),",
             " Memory (", usage$Memory_Usage_Percent, "%),",
             " Drive ", usage$Drive, " (", usage$Drive_Usage_Percent, "%)\n"))
}

sufficient_memory_to_run_this_check <- function(to_run_GB, step = "indexing", wiggle_room = 1,
                                                ref_path = path.expand(config()["ref"])) {
  system_info <- get_system_usage(detect_drive(ref_path))
  cat("-- System resource usage:\n")
  get_system_usage_one_liner(system_info)
  memory_on_computer <- system_info$Memory_Total_GB
  if ((system_info$Memory_Total_GB - wiggle_room) < to_run_GB) {
    message("Your ",step, " might fail, you specified ", to_run_GB,
            "GB max ram, and you only have ", memory_on_computer, "GB ram available.")
  }
}


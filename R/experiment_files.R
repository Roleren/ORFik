#' Get all library files in folder/folders of given types
#'
#' Will try to guess paired / unpaired wig, bed, bam files.
#'
#' Set pairedEndBam if you have paired end reads as a single bam file.
#' @inheritParams create.experiment
#' @importFrom tools file_ext
#' @return (data.table) All files found from types parameter.
#' With 2 extra column (logical), is it wig pairs, and paired bam files.
#' @keywords internal
findLibrariesInFolder <- function(dir, types, pairedEndBam = FALSE) {
  notDir <- !all(dir.exists(dir))
  if (notDir) stop(paste(dir[!dir.exists(dir)], "is not a existing directory!"))

  regex <- paste("\\.", types, collapse = "|", sep = "")
  # Find files in multiple dirs in correct order
  files <- pasteDir(unlist(lapply(dir, function(d)
    list.files(d, pattern = regex, full.names = TRUE))))

  files <- clean_file_extension(files, types)
  pairs <- find_pairs_from_folder(files)

  return(format_file_pairs_result(pairs, files, pairedEndBam) )
}

clean_file_extension <- function(files, types) {
  # Remove .bai bam index files etc
  fext <- file_ext(files)
  if (!(all(fext %in% types))) {
    files <- files[fext != "bai"]
    fext <- fext[fext != "bai"]

    compressed = fext %in% c("gzip", "gz", "bgz", "zip")
    if (any(compressed)) {
      fext[compressed] <-file_ext(file_path_sans_ext(files[compressed],
                                                     compression = FALSE))
    }
  }
  names(files) <- fext
  files <- files[fext %in% types]
  if (length(files) == 0) stop("Found no valid files in folder")
  return(files)
}

find_pairs_from_folder <- function(files, fext = names(files)) {
  stopifnot(length(files) == length(fext))
  # BigWig pairs
  bwg_files <- findNGSPairs(files[fext == "bigWig"], format = "bigWig")
  # Wig pairs
  wig_files <- findNGSPairs(files[fext == "wig"], format = "wig")
  # Bed pairs
  bed_pairs <- findNGSPairs(files[fext == "bed"], format = "bed")
  # Paired end double bam
  bam_pairs <- findNGSPairs(files[fext == "bam"],
                            f = c("_R1_00", "_F", "_Forward", "_forward"),
                            r = c("_R2_00", "_R", "_Reverse", "_reverse"),
                            format = "bam")
  pairs <- data.table()
  if (is(bwg_files, "data.table")) pairs <- rbind(pairs, bwg_files)
  if (is(wig_files, "data.table")) pairs <- rbind(pairs, wig_files)
  if (is(bed_pairs, "data.table")) pairs <- rbind(pairs, bed_pairs)
  if (is(bam_pairs, "data.table")) pairs <- rbind(pairs, bam_pairs)
  return(pairs)
}

format_file_pairs_result <- function(pairs, files, pairedEndBam) {
  pair_exists <- nrow(pairs) > 0
  if (pair_exists) {
    all_files_are_pairs <- nrow(pairs)*2 != length(files)
    if (all_files_are_pairs) { # if more than just matched pairs
      others <- files[!(files %in% c(pairs$forward, pairs$reverse))]
      file_dt <- data.table(forward = others, reverse = "", match = FALSE)
      files <- rbind(pairs, file_dt)
      if (any(pairedEndBam)) {
        if (length(pairedEndBam) == 1) {
          files[files$reverse == "",] <- "paired-end"
        } else if (length(pairedEndBam) == nrow(files)) {
          files[(files$reverse == "") & pairedEndBam,] <- "paired-end"
        } else stop("pairedEndBam must be either length 1 or
                    number of files in experiment")
      }
    } else files <- pairs
    if (nrow(files) == 0) stop("Found no valid files in folder")
  } else if (any(pairedEndBam)) {
    files <- data.table(forward = files, reverse = "", match = TRUE)
    files[pairedEndBam == TRUE,]$reverse <- "paired-end"
  }
  return(files)
}

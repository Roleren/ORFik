#' Create count table info for QC report
#'
#' The better the annotation / gtf used, the more results you get.
#' @inheritParams outputLibs
#' @inheritParams QCreport
#' @inheritParams countTable_regions
#' @return a data.table of the count info
#' @keywords internal
QC_count_tables <- function(df, out.dir, type = "ofst",
                            use_simplified_reads = TRUE,
                            force = TRUE, forceRemake = FALSE,
                            library.names = bamVarName(df),
                            BPPARAM = bpparam()) {
  stopifnot(is.logical(use_simplified_reads))

  outputLibs(df, chrStyle = findFa(df), type = type, force = force,
             library.names = library.names, BPPARAM = BPPARAM)
  # TODO: test if needed
  if (use_simplified_reads) {
    suppressMessages(convertLibs(df, out.dir = NULL,
                                 library.names = library.names, force = force))
  }


  # Make count tables
  message("--------------------------")
  dt_list <- countTable_regions(df, geneOrTxNames = "tx",
                                longestPerGene = FALSE,
                                out.dir = out.dir, lib.type = type,
                                library.names = library.names,
                                forceRemake = forceRemake,
                                BPPARAM = BPPARAM)
  return(invisible(NULL))
}

#' Create alignment feature statistcs
#'
#' Among others how much reads are in mRNA, introns, intergenic,
#' and check of reads from rRNA and other ncRNAs.
#' The better the annotation / gtf used, the more results you get.
#' @inheritParams outputLibs
#' @inheritParams QCreport
#' @return a data.table of the statistcs
#' @keywords internal
alignmentFeatureStatistics <- function(df, type = "ofst", force = TRUE,
                                       library.names = bamVarName(df),
                                       BPPARAM = bpparam()) {
  message("--------------------------")
  message("Making alignment statistics for lib:")
  message("- Loading annotation regions..")
  # Special regions rRNA etc..
  types <- c()
  # TODO: Check if there is a way to get this from txdb directly
  txdb <- loadTxdb(df)
  fa <- findFa(df)
  gff.df <- importGtfFromTxdb(txdb, stop.error = FALSE)
  if (is.null(gff.df)) warnings("No biotypes defined in GTF,",
                                " skiping biotype analysis!")
  types <- unique(gff.df$transcript_biotype)
  # The ncRNAs regions to check
  types <-types[types %in% c("Mt_rRNA", "snRNA", "snoRNA", "lincRNA", "miRNA",
                             "rRNA", "Mt_rRNA", "ribozyme", "Mt_tRNA")]

  # Helper function, sum countOverlaps with weight
  sCo <- function(region, lib) {
    weight <- "score"
    if (!(weight %in% colnames(mcols(lib))))
      weight <- NULL
    return(sum(countOverlapsW(region, lib, weight = weight)))
  }
  tx <- loadRegion(txdb, "tx")
  mrna <- loadRegion(txdb, "mrna")
  cds <- loadRegion(txdb, "cds")
  leaders <- loadRegion(txdb, "leaders")
  trailers <- loadRegion(txdb, "trailers")
  introns <- loadRegion(txdb, "introns")

  outputLibs(df, chrStyle = fa, type = type, force = force,
             library.names = library.names, BPPARAM = BPPARAM)
  libs <- library.names
  message("- Calculating alignment statistics..")
  finals <- bplapply(libs, function(s, sCo, tx, gff.df, libs, env) {
    message(s)
    lib <- get(s, envir = env)
    # Raw stats
    aligned_reads <- ifelse(!is.null(mcols(lib)$score),
                            sum(mcols(lib)$score), length(lib))
    res <- data.frame(Sample = s, Raw_reads = as.numeric(NA),
                      Trimmed_reads = as.numeric(NA),
                      Aligned_reads = aligned_reads)
    res$percentage_aligned_raw = 100 * (res$Aligned_reads / res$Raw_reads)

    # mRNA region stats
    tx_reads <- findOverlaps(lib, tx)
    tx_reads_indices <- unique(from(tx_reads))
    tx_reads_total <- sum(mcols(lib)$score[tx_reads_indices])
    mRNA_reads <- findOverlaps(lib, mrna)
    mRNA_reads_indices <- unique(from(mRNA_reads))
    mRNA_reads_total <- sum(mcols(lib)$score[mRNA_reads_indices])
    mRNA_proportion_covered <- length(unique(to(mRNA_reads)))
    # ncRNA, Intronic and intergenic are only using count not aligning to a mrna / transcript
    ncRNA_reads <- findOverlaps(lib, tx[!(names(tx) %in% names(mrna))])
    ncRNA_reads_total <- sum(mcols(lib)$score[unique(from(ncRNA_reads[!(from(ncRNA_reads) %in% mRNA_reads_indices)]))])

    intronic_reads <- findOverlaps(lib, introns)
    intronic_reads_total <- sum(mcols(lib)$score[unique(from(intronic_reads[!(from(intronic_reads) %in% tx_reads_indices)]))])
    intergenic_reads_total <- sum(mcols(lib[-unique(c(tx_reads_indices, unique(from(intronic_reads))))])$score)

    res_mrna <- data.table(mRNA = mRNA_reads_total,
                           LEADERS = sum(mcols(lib)$score[unique(from(findOverlaps(lib, leaders)))]),
                           CDS = sum(mcols(lib)$score[unique(from(findOverlaps(lib, cds)))]),
                           TRAILERs = sum(mcols(lib)$score[unique(from(findOverlaps(lib, trailers)))]),
                           Transcript = tx_reads_total,
                           ncRNA = ncRNA_reads_total,
                           Introns = intronic_reads_total,
                           Intergenic = intergenic_reads_total)

    res_mrna[,percentage_mrna_aligned := round(100* (mRNA / res$Aligned_reads), 2)]
    res_mrna[,percentage_tx_aligned := round(100* (Transcript / res$Aligned_reads), 2)]
    res_mrna[,ratio_cds_mrna := round(CDS / mRNA, 2)]
    res_mrna[, ratio_cds_leader := round(CDS / LEADERS, 2)]

    # Special region stats
    numbers <- c()
    for (t in types) {
      valids <- gff.df[grep(x = gff.df$transcript_biotype, pattern = t)]
      numbers <- c(numbers, sCo(tx[unique(valids$transcript_id)], lib))
    }

    # Lib width distribution, after soft.clip
    widths <- round(summary(readWidths(lib)))
    res_widths <- data.frame(matrix(widths, nrow = 1))
    colnames(res_widths) <- paste(names(widths), "read length")
    res_final <- cbind(res, res_widths, res_mrna)
    if (length(numbers) > 0) {
      res_extra <- data.frame(matrix(numbers, nrow = 1))
      colnames(res_extra) <- c(types)
      res_final <- cbind(res_final, res_extra)
    }
    return(res_final)
  }, sCo = sCo, tx = tx, gff.df = gff.df,
  libs = libs, env = envExp(df), BPPARAM = BPPARAM)

  return(rbindlist(finals))
}

#' Add trimming info to QC report
#'
#' Only works if alignment was done using ORFik with STAR.
#' @inheritParams QCreport
#' @param finals a data.table with current output from QCreport
#' @param alignment_folder character, default: \code{libFolder(df, "unique")}.
#' All unique folders. trim_folders should then be relative as:
#' file.path(alignment_folder, "..", "trim/")
#' @return a data.table of the update finals object with trim info
#' @keywords internal
trim_detection <- function (df, finals,
                            alignment_folder = libFolder(df, "unique")) {
  found_data <- FALSE
  trim_folders <- file.path(alignment_folder, "..", "trim/")
  if (all(dir.exists(trim_folders))) {
    message("Create raw read counts")
    raw_data <- try(rbindlist(lapply(trim_folders, trimming.table)), silent = TRUE)
    if (!is(raw_data, "try-error")) found_data <- TRUE
  }
  if (found_data) {
    matches <- unlist(lapply(X = df$filepath, function(x) {
      match <- sapply(raw_data$raw_library, function(p) grep(pattern = p, x, fixed = TRUE))
      if (isEmpty(match)) NA else names(unlist(match))[1]
    }))
    matches <- match(matches, raw_data$raw_library)
    if (nrow(finals) < nrow(raw_data)) {
      message("A subset of raw data will be used.")
    }
    if (any(is.na(matches))) {
      message("Could not find raw read count for some data, setting to NA.")
      trim_detection_message(df, raw_data, finals, matches)
    }
    raw_data <- raw_data[matches, ]
    class(finals$Raw_reads) <- "numeric"
    class(finals$Trimmed_reads) <- "numeric"
    finals$Raw_reads <- raw_data$raw_reads
    finals$Trimmed_reads <- raw_data$trim_reads
    finals$percentage_aligned_raw = round(100 * (finals$Aligned_reads/finals$Raw_reads),
                                          4)
  } else {
    message("Could not find raw read counts of data, setting to NA")
    message(paste0("No existing fastp json files in folder:", trim_folders))
  }
  return(finals)
}

trim_detection_message <- function(df, raw_data, finals, order) {
  message(paste("Alignment folder had", nrow(finals), "libraries"))
  message(paste("ORFik experiment has", nrow(df), "libraries"))
  message(paste("Trim folder had", nrow(raw_data), "libraries"))
  print(paste(c("Matches in the order:", order), collapse = " "))
  print("Trim data:")
  print(raw_data)
  print("ORFik experiment paths:")
  print(df$filepath)
  message("unexpected behavior, did you delete or merge any files?",
          "else report this bug on ORFik issues page!")
  message("Could not find raw read counts of data, setting to NA")
  return(invisible(NULL))
}

#' Load ORFik QC Statistics report
#'
#' Loads the pre / post alignment statistcs made in ORFik.
#'
#' The ORFik QC uses the aligned files (usually bam files),
#' fastp and STAR log files
#' combined with annotation to create relevant statistics.
#' @inheritParams QCreport
#' @param path path to QC statistics report, default:
#' file.path(dirname(df$filepath[1]), "/QC_STATS/STATS.csv")
#' @family QC report
#' @return data.table of QC report or NULL if not exists
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' ## First make QC report
#' # QCreport(df)
#' # stats <- QCstats(df)
QCstats <- function(df, path = file.path(QCfolder(df), "STATS.csv")) {
  if (!file.exists(path)) {
    message("No QC report made, run QCreport. Or wrong path given.")
    return(invisible(NULL))
  }
  return(fread(path, header = TRUE))
}

#' Make table of readlengths
#'
#' Summarizing all libraries in experiment,
#' make a table of proportion of read lengths.
#' @param output.dir NULL or character path, default: NULL, plot not saved to disc.
#' If defined saves plot to that directory with the name "./readLengths.csv".
#' @inheritParams heatMapRegion
#' @inheritParams QC_count_tables
#' @return a data.table object of the the read length data with columns:
#' \code{c("sample", "sample_id", "read length", "counts",
#'  "counts_per_sample", "perc_of_counts_per_sample")}
#' @keywords internal
readLengthTable <- function(df, output.dir = NULL, type = "ofst",
                            force = TRUE,
                            library.names = bamVarName(df),
                            BPPARAM = bpparam()) {
  file.name <- file.path(output.dir, "readLengths.csv")
  if (file.exists(file.name)) {
    message("Using previously stored readlengths in QC_STATS folder!")
    return(fread(file.name, header = TRUE))
  }

  outputLibs(df, type = type, force = force, library.names = library.names, BPPARAM = BPPARAM)

  dt_read_lengths <- rbindlist(lapply(seq_along(library.names), function(sample_id) {
    lib <- library.names[sample_id]
    data.table(sample = lib, sample_id, table(readWidths(get(lib, envir = envExp(df)))))
  }))


  colnames(dt_read_lengths) <- c("sample", "sample_id", "read length", "counts")
  dt_read_lengths[, counts_per_sample := sum(counts), by = sample_id]
  dt_read_lengths[, perc_of_counts_per_sample :=
                    (counts / counts_per_sample)*100,
                  by = sample_id]

  if (!is.null(output.dir)) {
    fwrite(dt_read_lengths, file.name)
  }
  return(dt_read_lengths)
}

#' Find shifted Ribo-seq frame distributions
#'
#' Per library: get coverage over CDS per frame per readlength
#' Return as data.datable with information and best frame found.
#' Can be used to automize re-shifting of read lengths (find read lengths
#' where frame 0 is not the best frame over the entire cds)
#' @inheritParams RiboQC.plot
#' @param orfs GRangesList, default loadRegion(df, part = "cds")
#' @param libraries a list of loaded libraries, default:
#'  outputLibs(df, type = type, output.mode = "envirlist")
#' @return data.table with columns: fraction (library) frame (0, 1, 2)
#' score (coverage) length (read length)
#'    percent (coverage percentage of library)
#'    percent_length (coverage percentage of library and length)
#'    best_frame (TRUE/FALSE, is this the best frame per length)
#' @export
#' @examples
#' df <- ORFik.template.experiment()[9,]
#' dt <- orfFrameDistributions(df, BPPARAM = BiocParallel::SerialParam())
#' ## Check that frame 0 is best frame for all
#' all(dt[frame == 0,]$best_frame)
orfFrameDistributions <- function(df, type = "pshifted", weight = "score",
                                  orfs = loadRegion(df, part = "cds"),
                                  libraries = outputLibs(df, type = type, output.mode = "envirlist"),
                                  BPPARAM = BiocParallel::bpparam()) {
  stopifnot(is(libraries, "list"))
  cds <- orfs
  # Frame distribution over all
  frame_sum_per1 <- bplapply(libraries, FUN = function(lib, cds, weight) {
    total <- regionPerReadLength(cds, lib,
                                 withFrames = TRUE, scoring = "frameSumPerL",
                                 weight = weight, drop.zero.dt = TRUE)
    total[, length := fraction]
    #hits <- get(lib, mode = "S4")[countOverlaps(get(lib, mode = "S4"), cds) > 0]
    name <- attr(lib, "name_short")
    total[, fraction := rep(name, nrow(total))]
  }, cds = cds, weight = weight, BPPARAM = BPPARAM)
  frame_sum_per <- rbindlist(frame_sum_per1)

  frame_sum_per$frame <- as.factor(frame_sum_per$frame)
  frame_sum_per$fraction <- as.factor(frame_sum_per$fraction)
  frame_sum_per[, percent := (score / sum(score))*100, by = fraction]
  frame_sum_per[, percent_length := (score / sum(score))*100, by = .(fraction, length)]
  frame_sum_per[, best_frame := (percent_length / max(percent_length)) == 1, by = .(fraction, length)]
  frame_sum_per[, fraction := factor(fraction, levels = names(libraries),
                                     labels = gsub("^RFP_", "", names(libraries)), ordered = TRUE)]

  frame_sum_per[, fraction := factor(fraction, levels = unique(fraction), ordered = TRUE)]
  frame_sum_per[]
  return(frame_sum_per)
}

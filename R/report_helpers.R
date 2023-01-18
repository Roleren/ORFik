#' Create count table info for QC report
#'
#' The better the annotation / gtf used, the more results you get.
#' @inheritParams outputLibs
#' @inheritParams QCreport
#' @return a data.table of the count info
#' @keywords internal
QC_count_tables <- function(df, out.dir, type = "ofst",
                            BPPARAM = bpparam()) {
  outputLibs(df, findFa(df), type = type, BPPARAM = BPPARAM)
  # TODO: test if needed
  suppressMessages(convertLibs(df, NULL)) # Speedup by reducing unwanted information

  # Make count tables
  message("--------------------------")
  dt_list <- countTable_regions(df, geneOrTxNames = "tx",
                                longestPerGene = FALSE,
                                out.dir = out.dir, lib.type = type,
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
alignmentFeatureStatistics <- function(df, type = "ofst",
                                       BPPARAM = bpparam()) {
  message("--------------------------")
  message("Making alignment statistics for lib:")
  # Special regions rRNA etc..
  types <- c()
  # TODO: Check if there is a way to get this from txdb directly
  txdb <- loadTxdb(df)
  fa <- findFa(df)
  outputLibs(df, fa, type = type, BPPARAM = BPPARAM)
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
  libs <- bamVarName(df)
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
#' @return a data.table of the update finals object with trim info
#' @keywords internal

trim_detection <- function (df, finals) {
  out.dirs <- unique(dirname(df$filepath))
  trim_folders <- file.path(out.dirs, "..", "trim/")
  if (all(dir.exists(trim_folders))) {
    message("Create raw read counts")
    raw_data <- rbindlist(lapply(trim_folders, trimming.table))
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
    message(paste0("No folder called:", trim_folder))
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
                            BPPARAM = bpparam()) {
  file.name <- file.path(output.dir, "readLengths.csv")
  if (file.exists(file.name)) {
    message("Using previously stored readlengths in QC_STATS folder!")
    return(fread(file.name, header = TRUE))
  }

  outputLibs(df, type = type, BPPARAM = BPPARAM)
  dt_read_lengths <- data.table(); sample_id <- 1
  for(lib in bamVarName(df)) {
    dt_read_lengths <- rbind(dt_read_lengths,
                             data.table(sample = lib, sample_id,
                              table(readWidths(get(lib, envir = envExp(df))))))
    sample_id <- sample_id + 1
  }

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
#' @return data.table with columns: fraction (library) frame (0, 1, 2)
#' score (coverage) length (read length)
#'    percent (coverage percentage of library)
#'    percent_length (coverage percentage of library and length)
#'    best_frame (TRUE/FALSE, is this the best frame per length)
#' @export
#' @examples
#' df <- ORFik.template.experiment()[3,]
#' dt <- orfFrameDistributions(df, BPPARAM = BiocParallel::SerialParam())
#' ## Check that frame 0 is best frame for all
#' all(dt[frame == 0,]$best_frame)
orfFrameDistributions <- function(df, type = "pshifted", weight = "score",
                                  BPPARAM = BiocParallel::bpparam()) {
  outputLibs(df, type = type)
  cds <- loadRegion(df, part = "cds")
  libs <- bamVarName(df)

  # Frame distribution over all
  frame_sum_per1 <- bplapply(libs, FUN = function(lib, cds, weight, env) {
    total <- regionPerReadLength(cds, get(lib, mode = "S4", envir = env),
                                 withFrames = TRUE, scoring = "frameSumPerL",
                                 weight = weight, drop.zero.dt = TRUE)
    total[, length := fraction]
    #hits <- get(lib, mode = "S4")[countOverlaps(get(lib, mode = "S4"), cds) > 0]
    total[, fraction := rep(lib, nrow(total))]
  }, cds = cds, weight = weight, env = envExp(df),BPPARAM = BPPARAM)
  frame_sum_per <- rbindlist(frame_sum_per1)

  frame_sum_per$frame <- as.factor(frame_sum_per$frame)
  frame_sum_per$fraction <- as.factor(frame_sum_per$fraction)
  frame_sum_per[, percent := (score / sum(score))*100, by = fraction]
  frame_sum_per[, percent_length := (score / sum(score))*100, by = .(fraction, length)]
  frame_sum_per[, best_frame := (percent_length / max(percent_length)) == 1, by = .(fraction, length)]
  frame_sum_per[, fraction := factor(fraction, levels = libs,
                                     labels = bamVarName(df, skip.libtype = TRUE), ordered = TRUE)]

  frame_sum_per[, fraction := factor(fraction, levels = unique(fraction), ordered = TRUE)]
  frame_sum_per[]
  return(frame_sum_per)
}

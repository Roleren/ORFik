#' Create count table info for QC report
#'
#' The better the annotation / gtf used, the more results you get.
#' @inheritParams QCreport
#' @return a data.table of the count info
QC_count_tables <- function(df, out.dir, BPPARAM = bpparam()) {
  txdb <- loadTxdb(df)
  loadRegions(txdb, parts = c("mrna", "leaders", "cds", "trailers", "tx"))
  outputLibs(df, leaders, type = "ofst", BPPARAM = BPPARAM)
  libs <- bamVarName(df)
  # Update this to use correct
  convertLibs(df, NULL) # Speedup by reducing unwanted information

  # Make count tables
  dt_list <- countTable_regions(df, geneOrTxNames = "tx",
                                longestPerGene = FALSE,
                                out.dir = out.dir,
                                BPPARAM = BPPARAM)
  # Special regions rRNA etc..
  gff.df <- importGtfFromTxdb(txdb)
  types <- unique(gff.df$transcript_biotype)
  types <-types[types %in% c("Mt_rRNA", "snRNA", "snoRNA", "lincRNA", "miRNA",
                             "rRNA", "Mt_rRNA", "ribozyme", "Mt_tRNA")]
  # Put into csv, the standard stats
  message("Making alignment statistics for lib:")
  sCo <- function(region, lib) {
    weight <- "score"
    if (!(weight %in% colnames(mcols(lib))))
      weight <- NULL
    return(sum(countOverlapsW(region, lib, weight = weight)))
  }
  finals <- bplapply(libs, function(s, dt_list, sCo, tx, gff.df, libs) {
    message(s)
    lib <- get(s)
    # Raw stats
    res <- data.frame(Sample = s, Raw_reads = as.numeric(NA),
                      Trimmed_reads = as.numeric(NA),
                      Aligned_reads = length(lib))
    res$percentage_aligned_raw = 100 * (res$Aligned_reads / res$Raw_reads)

    # mRNA region stats
    index <- which(s == libs)
    res_mrna <- data.table(mRNA = colSums(assay(dt_list[["mrna"]]))[index],
                           LEADERS = colSums(assay(dt_list[["leaders"]]))[index],
                           CDS = colSums(assay(dt_list[["cds"]]))[index],
                           TRAILERs = colSums(assay(dt_list[["trailers"]]))[index])
    res_mrna[,percentage_mrna_aligned := round(100* (mRNA / res$Aligned_reads), 6)]
    res_mrna[,ratio_cds_mrna := round(CDS / mRNA, 6)]
    res_mrna[, ratio_cds_leader := round(CDS / LEADERS, 6)]

    # Special region stats
    numbers <- sCo(tx, lib)
    for (t in types) {
      valids <- gff.df[grep(x = gff.df$transcript_biotype, pattern = t)]
      numbers <- c(numbers, sCo(tx[unique(valids$transcript_id)], lib))
    }

    res_extra <- data.frame(matrix(numbers, nrow = 1))
    colnames(res_extra) <- c("All_tx_types", types)

    # Lib width distribution, after soft.clip
    widths <- round(summary(readWidths(lib)))
    res_widths <- data.frame(matrix(widths, nrow = 1))
    colnames(res_widths) <- paste(names(widths), "read length")
    cbind(res, res_widths, res_mrna, res_extra)
  }, dt_list = dt_list, sCo = sCo, tx = tx, gff.df = gff.df,
     libs = libs, BPPARAM = BPPARAM)

  return(rbindlist(finals))
}

#' Add trimming info to QC report
#'
#' Only works if alignment was done using ORFik with STAR.
#' @inheritParams QCreport
#' @param finals a data.table with current output from QCreport
#' @return a data.table of the update finals object with trim info
trim_detection <- function(df, finals, out.dir) {
  # Update raw reads to real number
  # Needs a folder called trim
  trim_folder <- file.path(out.dir, "..", "trim/")
  if (dir.exists(trim_folder)) {
    message("Create raw read counts")

    raw_data <- trimming.table(trim_folder)
    # Verify rows have equal order as experiment
    order <- unlist(lapply(X = raw_data$raw_library,
                           function(p) grep(pattern = p, x = df$filepath, fixed = TRUE)))
    order <- unique(order)
    notMatch <-
      !all(seq(nrow(df)) %in% order) | length(seq(nrow(df))) != length(order)
    if (length(order) != nrow(raw_data)) {
      message(paste("ORFik experiment has", nrow(df), "libraries"))
      message(paste("Trim folder had", nrow(raw_data), "libraries"))
      print(paste(c("Matches in the order:", order), collapse = " "))
      print(raw_data)
      print("Input file paths:")
      print(df$filepath)
      message("unexpected behavior, did you delete or merge any files?",
           "else report this bug on ORFik github page!")
      message("Could not find raw read counts of data, setting to NA")

    } else if (notMatch) { # did not match all
      message("Could not find raw read count of your data, setting to NA")
    } else {
      class(finals$Raw_reads) <- "numeric"
      class(finals$Trimmed_reads) <- "numeric"
      finals[order,]$Raw_reads <- raw_data$raw_reads
      finals[order,]$Trimmed_reads <- raw_data$trim_reads
      finals$percentage_aligned_raw = round(100* (finals$Aligned_reads /
                                         finals$Raw_reads), 4)
    }
  } else {
    message("Could not find raw read counts of data, setting to NA")
    message(paste0("No folder called:", trim_folder))
  }
  return(finals)
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
#' df <- ORFik.template.experiment()[3,]
#' ## First make QC report
#' # QCreport(df)
#' # stats <- QCstats(df)
QCstats <- function(df, path = file.path(dirname(df$filepath[1]),
                                         "/QC_STATS/STATS.csv")) {
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
#' @param stats path to ORFik QC stats .csv file, or the experiment object.
#' @param output.dir NULL or character path, default: NULL, plot not saved to disc.
#' If defined saves plot to that directory with the name "./readLengths.csv".
#' @return a data.table object of the the read length data with columns:
#' \code{c("sample", "sample_id", "read length", "counts",
#'  "counts_per_sample", "perc_of_counts_per_sample")}
readLengthTable <- function(df, output.dir = NULL, type = "ofst",
                            BPPARAM = bpparam()) {
  file.name <- file.path(output.dir, "readLengths.csv")
  if (file.exists(file.name)) return(fread(file.name, header = TRUE))

  outputLibs(df, type = type, BPPARAM = BPPARAM)
  dt_read_lengths <- data.table(); sample_id <- 1
  for(lib in bamVarName(df)) {
    dt_read_lengths <- rbind(dt_read_lengths, data.table(sample = lib, sample_id, table(readWidths(get(lib)))))
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

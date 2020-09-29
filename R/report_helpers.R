#' Create count table info for QC report
#'
#' The better the annotation / gtf used, the more results you get.
#' @inheritParams QCreport
#' @return a data.table of the count info
QC_count_tables <- function(df, out.dir) {
  txdb <- loadTxdb(df)
  loadRegions(txdb, parts = c("mrna", "leaders", "cds", "trailers", "tx"))
  outputLibs(df, leaders, type = "ofst")
  libs <- bamVarName(df)
  # Update this to use correct
  convertLibs(df, NULL) # Speedup by reducing unwanted information

  # Make count tables
  dt_list <- countTable_regions(df, geneOrTxNames = "tx",
                                longestPerGene = FALSE,
                                out.dir = out.dir)
  # Special regions rRNA etc..
  gff.df <- importGtfFromTxdb(txdb)
  types <- unique(gff.df$transcript_biotype)
  types <-types[types %in% c("Mt_rRNA", "snRNA", "snoRNA", "lincRNA", "miRNA",
                             "rRNA", "Mt_rRNA", "ribozyme", "Mt_tRNA")]
  # Put into csv, the standard stats
  message("Making summary counts for lib:")
  sCo <- function(region, lib) {
    weight <- "score"
    if (!(weight %in% colnames(mcols(lib))))
      weight <- NULL
    return(sum(countOverlapsW(region, lib, weight = weight)))
  }
  for (s in libs) { # For each library
    message(s)
    lib <- get(s)
    # Raw stats
    res <- data.frame(Sample = s, Raw_reads = NA,
                      Aligned_reads = length(lib))
    res$ratio_aligned_raw = res$Aligned_reads / res$Raw_reads

    # mRNA region stats
    res_mrna <- data.table(mRNA = colSums(assay(dt_list[["mrna"]]))[s],
                           LEADERS = colSums(assay(dt_list[["leaders"]]))[s],
                           CDS = colSums(assay(dt_list[["cds"]]))[s],
                           TRAILERs = colSums(assay(dt_list[["trailers"]]))[s])
    res_mrna[,ratio_mrna_aligned := round(mRNA / res$Aligned_reads, 6)]
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

    final <- cbind(res, res_widths, res_mrna, res_extra)

    if (s == libs[1]) finals <- final
    else finals <- rbind(finals, final)
  }
  return(finals)
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
  if (dir.exists(paste0(out.dir, "/../trim/"))) {
    message("Create raw read counts")
    oldDir <- getwd()
    setwd(paste0(out.dir, "/../trim/"))
    raw_library <- system('ls *.json', intern = TRUE)
    lib_string <- 'grep -m 1 -h "total_reads" *.json | grep -Eo "[0-9]*"'
    raw_reads <- as.numeric(system(lib_string, intern = TRUE))
    raw_data <- data.table(raw_library, raw_reads)
    raw_data$raw_library <- gsub("report_",
                                 x = raw_data$raw_library, replacement = "")
    raw_data$raw_library <- gsub(".json",
                                 x = raw_data$raw_library, replacement = "")
    order <- unlist(lapply(X = raw_data$raw_library,
                           function(p) grep(pattern = p, x = df$filepath)))
    order <- unique(order)
    notMatch <-
      !all(seq(nrow(df)) %in% order) | length(seq(nrow(df))) != length(order)
    if (length(order) != nrow(raw_data)) {
      message(paste("ORFik experiment has", nrow(df), "libraries"))
      message(paste("Trim folder had", nrow(raw_data), "libraries"))
      print(paste(c("Matches in the order:", order), collapse = " "))
      print(raw_data)
      print(df$filepath)
      stop("unexpected behavior, report this bug on github page!")
    } else if (notMatch) { # did not match all
      message("Could not find raw read count of your data, setting to 20 M")
    } else {
      finals[order,]$Raw_reads <- raw_data$raw_reads
      finals$ratio_aligned_raw = round(finals$Aligned_reads /
                                         finals$Raw_reads, 4)
    }
    setwd(oldDir)
  } else {
    message("Could not find raw read counts of data, setting to NA")
    message(paste0("No folder called:", paste0(out.dir, "/../trim/")))
  }
  return(finals)
}

#' Make plot of ORFik QCreport
#'
#' From post-alignment QC relative to annotation, make a plot for all samples.
#' Will contain things like aligned_reads, read lengths, reads overlapping leaders,
#' cds, trailers, rRNA, tRNA etc.
#' @param stats path to ORFik QC stats .csv file, or the experiment object.
#' @param output.dir NULL or character path, default: NULL, plot not saved to disc.
#' If defined saves plot to that directory with the name "/STATS_plot.png".
#' @return ggplot object of the the statistics data
#' @importFrom data.table melt
#' @export
#' @examples
#' df <- ORFik.template.experiment()[3,]
#' ## First make QC report
#' # QCreport(df)
#' ## Now you can get plot
#' # QCstats.plot(df)
QCstats.plot <- function(stats, output.dir = NULL) {
  if (is(stats, "experiment")) {
    stats <- QCstats(stats)
    if (is.null(stats))
      stop("No QC report made for experiment, run ORFik QCreport")
  } else {
    stats <- fread(stats)
  }
  if (colnames(stats)[1] == "V1") colnames(stats)[1] <- "sample_id"

  dt_plot <- melt(stats, id.vars = c("Sample", "sample_id"))
  gg_STAT <- ggplot(dt_plot, aes(x=sample_id, y = value, group = Sample, fill = Sample)) +
    geom_bar(aes(color = Sample), stat="identity", position=position_dodge())+
    scale_fill_brewer(palette="Paired")+
    ylab("Annotation & Alignment feature, value") +
    xlab("Samples") +
    facet_wrap(  ~ variable, scales = "free") +
    theme_minimal()
  if (!is.null(output.dir)) {
    ggsave(paste0(output.dir, "/STATS_plot.png"), gg_STAT, width = 13, height = 8)
  }
  return(gg_STAT)
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
#' paste0(dirname(df$filepath[1]), "/QC_STATS/STATS.csv")
#' @family QC report
#' @return data.table of QC report or NULL if not exists
#' @export
#' @examples
#' df <- ORFik.template.experiment()[3,]
#' ## First make QC report
#' # QCreport(df)
#' # stats <- QCstats(df)
QCstats <- function(df, path = paste0(dirname(df$filepath[1]),
                                      "/QC_STATS/STATS.csv")) {
  if (!file.exists(path)) {
    message("No QC report made, run QCreport. Or wrong path given.")
    return(invisible(NULL))
  }
  return(fread(path, header = TRUE))
}

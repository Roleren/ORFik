#' A post Alignment quality control of reads
#'
#' From this report you will get a summary csv table, with distribution of
#' aligned reads, overlap of reads to transcript regions, like
#' leader, cds, trailer, tRNAs, rRNAs, snoRNAs etc.
#' It will also make you some correlation plots and meta coverage plots,
#' so you get a good understanding of how good the quality of your NGS
#' data production + aligner step were.
#' You will also get count tables over mrna, leader, cds and trailer
#' separately.
#'
#'
#' Everything will be outputed in the directory of your NGS data,
#' inside the folder QC_STATS/, relative to data location in 'df'
#'
#' To make a ORFik experiment, see ?ORFik::experiment
#' @param df an ORFik \code{\link{experiment}}
#' @param out.dir optional output directory, default: dirname(df$filepath[1])
#' @return NULL (objects stored to disc)
#' @export
#' @examples
#' # 1. Pick directory
#' dir <- system.file("extdata", "", package = "ORFik")
#' # 2. Pick an experiment name
#' exper <- "ORFik"
#' # 3. Pick .gff/.gtf location
#' txdb <- system.file("extdata", "annotations.gtf", package = "ORFik")
#' template <- create.experiment(dir = dir, exper, txdb = txdb,
#'                               viewTemplate = FALSE)
#' template$X5[6] <- "heart" # <- fix non unique row
#' # read experiment
#' df <- read.experiment(template)
#' # Save with: save.experiment(df, file = "path/to/save/experiment.csv")
#'
#' # ORFikQC(df)
ORFikQC <- function(df, out.dir = dirname(df$filepath[1])) {
  # When experiment is ready, everything down from here is automatic
  message("Started ORFik QC report:")
  validateExperiments(df)

  exp_dir <- out.dir
  stats_folder <- pasteDir(exp_dir, "/QC_STATS/")
  if (!dir.create(stats_folder)) {
    if (!dir.exists(stats_folder)) stop("Could not create directory!")
  }

  txdb <- loadTxdb(df)
  loadRegions(txdb, parts = c("mrna", "leaders", "cds", "trailers", "tx"))
  outputLibs(df, leaders)

  # Special regions rRNA etc..
  if (!(file_ext(df@txdb) %in% c("gtf", "gff", "gff3", "gff2"))) {
    # Try to use reference of sqlite
    if (!(file_ext(metadata(txdb)[3,2]) %in%
          c("gtf", "gff", "gff3", "gff2"))) {
      stop("Could not find valid gtf / gff file, only data base object!")
    }
    gff.df <- import(metadata(txdb)[3,2])
  } else gff.df <- import(df@txdb)
  types <- unique(gff.df$transcript_biotype)
  types <-types[types %in% c("Mt_rRNA", "snRNA", "snoRNA", "lincRNA", "miRNA",
                             "rRNA", "ribozyme", "Mt_tRNA")]

  # Make count tables
  message("Making count tables for region:")
  countDir <- paste0(stats_folder, "countTable_")
  for (region in c("mrna", "leaders", "cds", "trailers")) {
    message(region)
    path <- paste0(countDir,region)
    dt <- makeSummarizedExperimentFromBam(df, region = region,
                                          geneOrTxNames = "tx",
                                          longestPerGene = FALSE,
                                          saveName = path)
    assign(paste0("ct_", region), colSums(assay(dt)))
  }

  # Put into csv, the standard stats
  libs <- bamVarName(df)
  message("Making summary counts for lib:")
  for (s in libs) { # For each library
    message(s)
    lib <- get(s)
    # Raw stats
    res <- data.frame(Sample = s, Raw_reads = 2e7,
                      Aligned_reads = length(lib))
    res$ratio_aligned_raw = res$Aligned_reads / res$Raw_reads

    # mRNA region stats
    sCo <- function(region, lib) {
      return(sum(countOverlaps(region, lib)))
    }
    res_mrna <- data.table(mRNA = get("ct_mrna")[lib],
                           LEADERS = get("ct_leaders")[lib],
                           CDS = get("ct_cds")[lib],
                           TRAILERs = get("ct_trailers")[lib])
    res_mrna[,ratio_mrna_aligned := mRNA / res$Aligned_reads]
    res_mrna[,ratio_cds_mrna := CDS / mRNA]
    res_mrna[, ratio_cds_leader := CDS / LEADERS]

    # Special region stats
    numbers <- sCo(tx, lib)
    for (t in types) {
      valids <- gff.df[grep(x = gff.df$transcript_biotype, pattern = t)]
      numbers <- c(numbers, sCo(tx[unique(valids$transcript_id)], lib))
    }

    res_extra <- data.frame(matrix(numbers, nrow = 1))
    colnames(res_extra) <- c("All_tx_types", types)

    # Lib width distribution, after soft.clip
    widths <- summary(readWidths(lib))
    res_widths <- data.frame(matrix(widths, nrow = 1))
    colnames(res_widths) <- names(widths)

    final <- cbind(res, res_widths, res_mrna, res_extra)

    if (s == libs[1]) finals <- final
    else finals <- rbind(finals, final)
  }

  # Update raw reads to real number
  # Needs a folder called trim
  if (dir.exists(paste0(exp_dir, "../trim/"))) {
    message("Create raw read counts")
    oldDir <- getwd()
    setwd(paste0(exp_dir, "../trim/"))
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
    notMatch <-
      !all(seq(nrow(df)) %in% order) | length(seq(nrow(df))) != length(order)
    if (notMatch) { # did not match all
      message("Could not find raw read count of your data, setting to 20 M")
    } else {
      finals[order,]$Raw_reads <- raw_data$raw_reads
      finals$ratio_aligned_raw = finals$Aligned_reads / finals$Raw_reads
    }
    setwd(oldDir)
  } else {
    message("Could not find raw read counts of data, setting to 20 M")
    message(paste0("No folder called:", paste0(exp_dir, "../trim/")))
  }

  write.csv(finals, file = pasteDir(stats_folder, "STATS.csv"))

  QCplots(df, mrna, stats_folder)

  message(paste("Everything done, saved QC to:", stats_folder))
}

#' Correlation and coverage plots for ORFikQC
#'
#' Output will be stored in same folder as the
#' libraries in df.
#' @param df an ORFik \code{\link{experiment}}
#' @param region (default: mrna), make raw count matrices of
#' whole mrnas or one of (leaders, cds, trailers)
#' @param stats_folder directory to save
#' @return NULL (objects stored to disc)
#' @importFrom GGally ggpairs
#' @importFrom AnnotationDbi metadata
QCplots <- function(df, region, stats_folder) {
  message("Making QC plots:")
  message("- Correlation plots")
  # Load fpkm values
  data_for_pairs <- makeSummarizedExperimentFromBam(df, region = region,
                                                    geneOrTxNames = "tx",
                                                    longestPerGene = FALSE,
                                                    type = "fpkm")
  paired_plot <- ggpairs(as.data.frame(data_for_pairs),
                         columns = 1:ncol(data_for_pairs))
  ggsave(pasteDir(stats_folder, "cor_plot.png"), paired_plot,
         height=400, width=400, units = 'mm', dpi=300)
  paired_plot <- ggpairs(as.data.frame(log2(data_for_pairs)),
                         columns = 1:ncol(data_for_pairs))
  ggsave(pasteDir(stats_folder, "cor_plot_log2.png"), paired_plot,
         height=400, width=400, units = 'mm', dpi=300)

  # window coverage over mRNA regions
  message("- Meta coverage plots")
  txNames <- filterTranscripts(df, 100, 100, 100, longestPerGene = FALSE)
  transcriptWindow(leaders[txNames], get("cds", mode = "S4")[txNames],
                   trailers[txNames], df = df, outdir = stats_folder,
                   allTogether = TRUE)
}

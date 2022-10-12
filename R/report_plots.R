#' Make plot of ORFik QCreport
#'
#' From post-alignment QC relative to annotation, make a plot for all samples.
#' Will contain among others read lengths, reads overlapping leaders,
#' cds, trailers, mRNA / rRNA etc.
#' @param stats the experiment object or path to custom ORFik QC folder where a file
#' called "STATS.csv" is located.
#' @param output.dir NULL or character path, default: NULL, plot not saved to disc.
#' If defined saves plot to that directory with the name "/STATS_plot.pdf".
#' @param plot.ext character, default: ".pdf". Alternatives: ".png" or ".jpg".
#' @param as_gg_list logical, default FALSE. Return as a list of ggplot objects
#'  instead of as a grob. Gives you the ability to modify plots more directly.
#' @return the plot object, a grob of ggplot objects of the the statistics data
#' @importFrom data.table melt
#' @importFrom gridExtra grid.arrange
#' @export
#' @examples
#' df <- ORFik.template.experiment()[3,]
#' ## First make QC report
#' # QCreport(df)
#' ## Now you can get plot
#' # QCstats.plot(df)
QCstats.plot <- function(stats, output.dir = NULL, plot.ext = ".pdf",
                         as_gg_list = FALSE) {
  if (is(stats, "experiment")) {
    path <- QCfolder(stats)
    stats <- QCstats(stats)
    if (is.null(stats))
      stop("No QC report made for experiment, run ORFik QCreport")
  } else { # From path to QC folder
    path <- stats
    stats <- fread(file.path(stats, "STATS.csv"))
  }
  if (colnames(stats)[1] == "V1") colnames(stats)[1] <- "sample_id"

  temp_theme <-  theme(legend.text=element_text(size=8),
                       legend.key.size = unit(0.3, "cm"),
                       plot.title = element_text(size=11),
                       strip.text.x = element_blank(),
                       panel.grid.minor = element_blank())

  stats$sample_id <-  factor(stats$Sample,
                             labels = as.character(seq(length(stats$Sample))),
                             levels = stats$Sample)
  stats$Sample <-  factor(stats$Sample, levels = stats$Sample)
  colnames(stats) <- gsub("percentage", "%", colnames(stats))
  # Update all values to numeric
  stats[,(seq.int(ncol(stats))[-c(1,2)]):= lapply(.SD, as.numeric),
        .SDcols = seq.int(ncol(stats))[-c(1,2)]]
  dt_plot <- melt(stats, id.vars = c("Sample", "sample_id"))

  step_counts <- c("mRNA", "ncRNA", "Introns", "Intergenic")
  stat_regions <- colnames(stats)[c(which(colnames(stats) %in% step_counts))]
  dt_STAT <- dt_plot[(variable %in% stat_regions),]
  dt_STAT_normalized <- copy(dt_STAT)
  dt_STAT_normalized[, sample_total := sum(value), by = .(sample_id)]
  dt_STAT_normalized[, percentage := (value / sample_total)*100]

  gg_STAT <- ggplot(data = dt_STAT_normalized, aes(x = sample_id, y = percentage)) +
    geom_bar(aes(fill = variable), stat="identity", position = "stack")+
    theme_minimal() +
    ylab("% content") +
    scale_y_continuous(breaks = c(50, 100)) +
    temp_theme

  # Read lengths
  dt_read_lengths <- readLengthTable(NULL, output.dir = path)
  dt_read_lengths <- dt_read_lengths[perc_of_counts_per_sample > 1, ]
  dt_read_lengths[, sample_id := as.factor(sample_id)]
  gg_read_lengths <- ggplot(dt_read_lengths, aes(y = perc_of_counts_per_sample, x = `read length`, fill = sample_id)) +
    geom_bar(stat="identity", position = "dodge")+
    ylab("% counts") +
    facet_wrap(  ~ sample_id, nrow = length(levels(dt_read_lengths$sample_id))) +
    scale_y_continuous(breaks = c(15, 30)) +
    theme_minimal() +
    temp_theme

  # mRNA regions
  mRNA_regions <- colnames(stats)[colnames(stats) %in% c("LEADERS", "CDS", "TRAILERs")]
  dt_mRNA_regions <- dt_plot[(variable %in% mRNA_regions),]

  gg_mRNA_regions <- ggplot(dt_mRNA_regions, aes(y = value, x = sample_id)) +
    geom_bar(aes(fill = variable), stat="identity", position = "fill")+
    ylab("ratio") +
    scale_y_continuous(breaks = c(0.5, 1.0)) +
    theme_minimal() +
    temp_theme

  lay <- rbind(c(2),
               c(2),
               c(2),
               c(1,3),
               c(1,3))
  plot_list <- list(gg_STAT, gg_read_lengths, gg_mRNA_regions)
  final <- gridExtra::grid.arrange(grobs = plot_list, layout_matrix = lay)

  if (!is.null(output.dir)) {
    ggsave(file.path(output.dir, paste0("STATS_plot", plot.ext)), final, width = 13,
           height = 8, dpi = 300)
  }
  if (as_gg_list) return(plot_list)
  return(final)
}

#' Correlation plots between all samples
#'
#' Get correlation plot of raw counts and/or log2(count + 1) over
#' selected region in: c("mrna", "leaders", "cds", "trailers")\cr\cr
#' Note on correlation: Pearson correlation, using pairwise observations
#' to fill in NA values for the covariance matrix.
#' @inheritParams QCstats.plot
#' @inheritParams QCplots
#' @inheritParams cor_table
#' @inheritParams cor_plot
#' @param output.dir directory to save to, named : cor_plot,
#' cor_plot_log2 and/or cor_plot_simple with either .pdf or .png
#' @param type which value to use, "fpkm", alternative "counts".
#' @param height numeric, default 400 (in mm)
#' @param width numeric, default 400 (in mm)
#' @param size numeric, size of dots, default 0.15. Deprecated.
#' @param data_for_pairs a data.table from ORFik::countTable of counts wanted.
#' Default is fpkm of all mRNA counts over all libraries.
#' @return invisible(NULL) / if as_gg_list is TRUE, return a list of raw plots.
correlation.plots <- function(df, output.dir,
                              region = "mrna", type = "fpkm",
                              height = 400, width = 400, size = 0.15, plot.ext = ".pdf",
                              complex.correlation.plots = TRUE,
                              data_for_pairs = countTable(df, region, type = type),
                              as_gg_list = FALSE, text_size = 4,
                              method = c("pearson", "spearman")[1]) {
  message("- Correlation plots")
  if (nrow(df) == 1) { # Avoid error from ggplot2 backend
    message("-  Skipping correlation plots (only 1 sample)")
    return(invisible(NULL))
  }
  message("  - raw scaled fpkm (simple)")
  dt_cor <- cor_table(as.data.table(data_for_pairs), method = method)
  cor_plot1 <- cor_plot(dt_cor, text_size = text_size,
                        label_name = paste0(method, "\ncorrelation"))
  if (!is.null(output.dir)) {
    ggsave(pasteDir(output.dir, paste0("cor_plot_simple", plot.ext)), cor_plot1,
           height = height, width = width, units = 'mm', dpi = 300)
  }
  plot_list <- list(cor_plot1)
  if (complex.correlation.plots) warning("Complex correlation plots are for now deprecated!")

  if (as_gg_list) return(plot_list)
  return(invisible(NULL))
}

#' Get correlation between columns
#' @param dt a data.table
#' @param method c("pearson", "spearman")[1]
#' @param upper_triangle logical, default TRUE. Make lower triangle values NA.
#' @param decimals numeric, default 2. How many decimals for correlation
#' @param melt logical, default TRUE.
#' @param na.rm.melt logical, default TRUE. Remove NA values from melted table.
#' @return a data.table with 3 columns, Var1, Var2 and Cor
#' @importFrom data.table melt.data.table setcolorder
cor_table <- function(dt, method = c("pearson", "spearman")[1],
                      upper_triangle = TRUE, decimals = 2, melt = TRUE,
                      na.rm.melt = TRUE) {
  stopifnot(is(dt, "data.table"))
  cor <- round(cor(dt, method = method), decimals)
  # Get lower triangle of the correlation matrix
  get_lower_tri<-function(cormat){
    cormat[upper.tri(cormat)] <- NA
    return(cormat)
  }
  # Get upper triangle of the correlation matrix
  get_upper_tri <- function(cormat){
    cormat[lower.tri(cormat)]<- NA
    return(cormat)
  }
  if (upper_triangle) cor <- get_upper_tri(cor)
  if (melt) {
    cor <- suppressWarnings(data.table::melt.data.table(as.data.table(cor),
                                       value.name = "Cor",
                                       variable.name = "Var2", na.rm = FALSE))

    cor[, Var1 := rep(unique(Var2), times = length(unique(Var2)))]
    data.table::setcolorder(cor, neworder = c(3,1,2))
    if (na.rm.melt) cor <- cor[!is.na(Cor),]
  }

  return(as.data.table(cor))
}

#' Get correlation between columns
#' @param dt_cor a data.table, with column Cor
#' @param col colors c(low = "blue", high = "red", mid = "white", na.value = "white")
#' @param limit default (-1, 1), defined by:
#'  \code{c(ifelse(min(dt_cor$Cor, na.rm = TRUE) < 0, -1, 0), 1)}
#' @param midpoint midpoint of correlation values in label coloring.
#' @param label_name name of correlation method, default
#' \code{"Pearson Correlation"} with newline after Pearson.
#' @param legend.position default c(0.4, 0.7), other: "top", "right",..
#' @param text_size size of correlation numbers
#' @param legend.direction default "horizontal", or "vertical"
#' @return a ggplot (heatmap)
cor_plot <- function(dt_cor, col = c(low = "blue", high = "red", mid = "white", na.value = "white"),
                     limit = c(ifelse(min(dt_cor$Cor, na.rm = TRUE) < 0, -1, 0), 1),
                     midpoint = mean(limit), label_name = "Pearson\nCorrelation",
                     text_size = 4, legend.position = c(0.4, 0.7),
                     legend.direction = "horizontal") {
  stopifnot(is(dt_cor, "data.table"))
  if (!all(c("low", "high", "mid", "na.value") %in% names(col)))
    stop("'col' must have names low, high, mid and na.value")

  ggheatmap <- ggplot(dt_cor, aes(Var2, Var1, fill = Cor)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = col["low"], high = col["high"], mid = col["mid"], na.value = col["na.value"],
                         midpoint = midpoint, limit = limit, space = "Lab",
                         name = label_name) +
    theme_minimal()+ # minimal theme
    theme(axis.text.x = element_blank()) +
    coord_fixed()

  wh <- if (legend.direction == "horizontal") {c(7,1)} else c(1,7)
  # Add text boxes
  ggheatmap <- ggheatmap +
    geom_text(aes(Var2, Var1, label = Cor), color = "black", size = text_size) +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.ticks = element_blank(),
      legend.justification = c(1, 0),
      legend.position = legend.position,
      legend.direction = legend.direction)+
    guides(fill = guide_colorbar(barwidth = wh[1], barheight = wh[2],
                                 title.position = "top", title.hjust = 0.5))
  return(ggheatmap)
}

#' Simple PCA analysis
#'
#' Detect outlier libraries with PCA analysis.
#' Will output PCA plot of PCA component 1 (x-axis) vs
#' PCA component 2 (y-axis) for each library (colored by library),
#' shape by replicate. Will be extended to allow batch correction
#' in the future.
#' @inheritParams QCplots
#' @param output.dir default NULL, else character path to directory.
#' File saved as "PCAplot_(experiment name)(plot.ext)"
#' @param table data.table, default countTable(df, "cds", type = "fpkm"),
#' a data.table of counts per column (default normalized fpkm values).
#' @param title character, default "CDS fpkm".
#' @param subtitle character, default: \code{paste("Numer of genes:", nrow(table))}
#' @param color.by.group logical, default TRUE. Colors in PCA plot represent
#' unique library groups, if FALSE. Color each sample in seperate color
#' (harder to distinguish for > 10 samples)
#' @param return.data logical, default FALSE. Return data instead of plot
#' @return ggplot or invisible(NULL) if output.dir is defined or < 3 samples.
#' Returns data.table with PCA analysis if return.data is TRUE.
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' # Select only Ribo-seq and RNA-seq
#' pcaExperiment(df[df$libtype %in% c("RNA", "RFP"),])
pcaExperiment <- function(df, output.dir = NULL,
                          table = countTable(df, "cds", type = "fpkm"),
                          title = "PCA analysis by CDS fpkm",
                          subtitle = paste("Numer of genes/regions:", nrow(table)),
                          plot.ext = ".pdf",
                          return.data = FALSE,
                          color.by.group = TRUE) {
  if (nrow(df) < 3) {
    message("-  Skipping PCA analysis (< 3 samples)")
    return(invisible(NULL))
  }
  pca <- stats::prcomp(table, scale = FALSE)
  dt <- data.table(pca$rotation, keep.rownames = TRUE)
  if (color.by.group) {
    dt$sample <- bamVarName(df, skip.replicate = TRUE)
  } else dt$sample <- dt$rn

  if (any(df$rep > 1, na.rm = TRUE)) {
    dt$replicate <- df$rep
    dt$replicate[is.na(dt$replicate)] <- 1
    dt$replicate <- as.factor(dt$replicate)
  } else dt$replicate <- as.factor("1")

  plot <- ggplot(data = dt,
                 aes(x = PC1, y = PC2)) +
    geom_point(aes(shape=replicate, color = sample),
               size = 3, alpha = 0.8) +
    scale_fill_brewer() +
    ggtitle(title, subtitle) +
    theme_bw() +
    theme(legend.text=element_text(size=7))
  if(!is.null(output.dir)) {
    if (output.dir == "auto") {
      path <- file.path(dirname(df$filepath[1]), "QC_STATS",
                        paste0("PCAplot_", df@experiment, plot.ext))
    } else {
      path <- file.path(output.dir,
                        paste0("PCAplot_", df@experiment, plot.ext))
    }
    ggsave(path, plot,
           height = 4 + (nrow(df)*0.1), width = 5 + (nrow(df)*0.1))
    if (return.data) return(dt)
    return(invisible(NULL))
  }
  if (return.data) return(dt)
  return(plot)
}

#' Quality control for pshifted Ribo-seq data
#'
#' Combines several statistics from the pshifted reads into a plot:\cr
#' -1 Coding frame distribution per read length\cr
#' -2 Alignment statistics\cr
#' -3 Biotype of non-exonic pshifted reads\cr
#' -4 mRNA localization of pshifted reads\cr
#' @param type type of library loaded, default pshifted,
#'  warning if not pshifted might crash if too many read lengths!
#' @param width width of plot, default 6.6 (in inches)
#' @param height height of plot, default 4.5 (in inches)
#' @param weight which column in reads describe duplicates, default "score".
#' @param bar.position character, default "dodge". Should Ribo-seq frames
#' per read length be positioned as "dodge" or "stack" (on top of each other).
#' @inheritParams QCstats.plot
#' @inheritParams outputLibs
#' @return the plot object, a grob of ggplot objects of the the data
#' @importFrom ggplot2 theme
#' @export
#' @examples
#' df <- ORFik.template.experiment()
#' df <- df[3,] #lets only p-shift RFP sample at index 3
#' #shiftFootprintsByExperiment(df)
#' #RiboQC.plot(df)
RiboQC.plot <- function(df, output.dir = QCfolder(df),
                        width = 6.6, height = 4.5, plot.ext = ".pdf",
                        type = "pshifted", weight = "score", bar.position = "dodge",
                        as_gg_list = FALSE,
                        BPPARAM = BiocParallel::SerialParam(progressbar = TRUE)) {
  stats <- QCstats(df)
  stopifnot(bar.position %in% c("stack", "dodge"))

  if (colnames(stats)[1] == "V1") colnames(stats)[1] <- "sample_id"

  stats$sample_id <-  factor(stats$Sample,
                             labels = as.character(seq(length(stats$Sample))),
                             levels = stats$Sample)
  stats$Sample <-  factor(stats$Sample, levels = stats$Sample)
  colnames(stats) <- gsub("percentage", "%", colnames(stats))
  stats[,(seq.int(ncol(stats))[-c(1,2)]):= lapply(.SD, as.numeric),
        .SDcols = seq.int(ncol(stats))[-c(1,2)]]
  dt_plot <- melt(stats, id.vars = c("Sample", "sample_id"))

  step_counts <- c("mRNA", "rRNA")
  stat_regions <- colnames(stats)[c(which(colnames(stats) %in% step_counts))]
  dt_STAT <- dt_plot[(variable %in% stat_regions),]
  dt_STAT_normalized <- copy(dt_STAT)
  dt_STAT_normalized[, sample_total := sum(value), by = .(sample_id)]
  dt_STAT_normalized[, percentage := (value / sample_total)*100]

  # Assign themes
  temp_theme <-  theme(legend.text=element_text(size=8),
                       legend.key.size = unit(0.3, "cm"),
                       plot.title = element_text(size=11),
                       strip.text.x = element_blank(),
                       panel.grid.minor = element_blank())
  rm.minor <- theme(legend.text=element_text(size=8),
                    legend.key.size = unit(0.3, "cm"),
                    plot.title = element_text(size=11),
                    panel.grid.minor = element_blank())

  # frame distributions
  frame_sum_per <- orfFrameDistributions(df, type = type, weight = weight,
                                         BPPARAM = BPPARAM)
  gg_frame_per_stack <- ggplot(frame_sum_per, aes(x = length, y = percent)) +
    geom_bar(aes(fill = frame), stat="identity", position = bar.position)+
    scale_x_continuous(breaks = unique(frame_sum_per$length)) +
    theme_minimal() +
    rm.minor+
    xlab("read length") +
    ylab("percent") +
    facet_wrap(  ~ fraction, ncol = 1, scales = "free_y") +
    scale_y_continuous(breaks = c(15, 35))
  gg_frame_per_stack

  # content: all_tx_types > 1%
  all_tx_types <- which(colnames(stats) == "ratio_cds_leader") + 1
  all_tx_regions <- colnames(stats)[c(all_tx_types:length(colnames(stats)))]
  all_tx_regions <- c("mRNA", all_tx_regions)
  dt_all_tx_regions<- dt_plot[(variable %in% all_tx_regions),]
  dt_all_tx_regions[, sample_total := sum(value), by = .(sample_id)]
  dt_all_tx_regions[, percentage := (value / sample_total)*100]
  dt_all_tx_other <- dt_all_tx_regions[percentage < 1,]
  dt_all_tx_regions[percentage < 1, variable := "other"]
  gg_all_tx_regions <-
    ggplot(dt_all_tx_regions, aes(y = percentage, x = sample_id)) +
    geom_bar(aes(fill = variable), stat="identity", position = "stack")+
    ylab("percent") +
    xlab("sample id") +
    theme_minimal() +
    temp_theme +
    labs(fill = "tx. type")  +
    scale_y_continuous(breaks = c(50, 100))
  gg_all_tx_regions
  # content: all_tx_types <= 1%
  gg_all_tx_other <-
    ggplot(dt_all_tx_other, aes(y = percentage, x = sample_id)) +
    geom_bar(aes(fill = variable), stat="identity", position = "stack")+
    ylab("percent") +
    xlab("sample id") +
    theme_minimal() +
    temp_theme +
    labs(fill = "tx. other") +
    scale_fill_grey()  +
    scale_y_continuous(n.breaks = 3)
  gg_all_tx_other

  # Aligned reads
  dt_aligned <- dt_plot[(variable %in% "Aligned_reads"),]
  gg_all_mrna <- ggplot(dt_aligned, aes(y = value, x = sample_id)) +
    geom_bar(stat="identity", position = "stack")+
    ylab("alignments") +
    xlab("sample id") +
    theme_minimal() +
    temp_theme +
    labs(fill = "region") +
    scale_y_continuous(n.breaks = 3, labels = function(x) format(x, scientific = TRUE))
  gg_all_mrna
  # 5' UTR, CDS & 3' UTR
  mRNA_regions <- colnames(stats)[colnames(stats) %in% c("LEADERS", "CDS", "TRAILERs")]
  dt_mRNA_regions <- dt_plot[(variable %in% mRNA_regions),]
  dt_mRNA_regions[, variable := as.character(variable)]
  dt_mRNA_regions[variable == "LEADERS", variable := "5'UTR"]
  dt_mRNA_regions[variable == "TRAILERs", variable := "3'UTR"]
  dt_mRNA_regions[, variable := factor(as.character(dt_mRNA_regions$variable),
                                       levels = c(unique(dt_mRNA_regions$variable)),
                                       ordered = TRUE)]
  dt_mRNA_regions[, sample_total := sum(value), by = .(sample_id)]
  dt_mRNA_regions[, percentage := (value / sample_total)*100]
  gg_mRNA_regions <- ggplot(dt_mRNA_regions, aes(y = percentage, x = sample_id)) +
    geom_bar(aes(fill = variable), stat="identity", position = "stack")+
    ylab("percent") +
    xlab("sample id") +
    scale_y_continuous(breaks = c(50, 100)) +
    theme_minimal() +
    temp_theme +
    labs(fill = "region")

  lay <- rbind(c(2, 1),
               c(2, 3),
               c(2, 4),
               c(2, 5))

  plot_list <- list(gg_all_mrna, gg_frame_per_stack, gg_all_tx_regions, gg_all_tx_other, gg_mRNA_regions)
  final <- gridExtra::grid.arrange(grobs = plot_list, layout_matrix = lay)
  if (!is.null(output.dir)) {
    ggsave(file.path(output.dir, paste0("STATS_plot_Ribo-seq_Check", plot.ext)), final,
           width = width, height = height, dpi = 300)
  }
  if (as_gg_list) return(plot_list)
  return(final)
}


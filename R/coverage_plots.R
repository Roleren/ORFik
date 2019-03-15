#' Plot area around TIS for p-shifted reads
#'
#' Usefull to validate p-shifting is correct
#' Can be used for any coverage of region around a point, like TIS, TSS,
#' stop site etc.
#'
#'
#' @param hitMap a data.frame/data.table, given from metaWindow
#' (must have columns: position, (score or count) and frame)
#' @param length an integer (29), which length is this for?
#' @param region a character (start), either "start or "stop"
#' @param output character string (NULL), if set, saves the plot as pdf or png
#' to path given. If no format is given, is save as pdf.
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
#' @importFrom data.table setDF
#' @family coveragePlot
#' @examples
#' # An ORF
#' grl <- GRangesList(tx1 = GRanges("1", IRanges(1, 6), "+"))
#' # Ribo-seq reads
#' range <- IRanges(c(rep(1, 3), 2, 3, rep(4, 2), 5, 6), width = 1 )
#' reads <- GRanges("1", range, "+")
#' coverage <- coveragePerTiling(grl, reads, TRUE, as.data.table = TRUE,
#'                               withFrames = TRUE)
#' ORFik:::pSitePlot(coverage)
#'
#' # See vignette for more examples
#'
pSitePlot <- function(hitMap, length = 29, region = "start", output = NULL) {
  min <- min(hitMap$position)
  max <- max(hitMap$position)

  if (nrow(hitMap) > 30) {
    min <- min + 2
    max <- max - 1
  }
  by <- ifelse(length(hitMap$position) > 80, 3, 1)

  if (is.null(hitMap$score)) hitMap$score <- hitMap$count
  if (is(hitMap, "data.table")) setDF(hitMap)
  plot <- ggplot(hitMap, aes(x = factor(position), y = score,
                             fill = factor(frame))) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = paste("Length", length, "over", region, "of canonical CDS")) +
    xlab(paste("\nshift from first", region, "nucleotide [bp]")) +
    ylab("Averaged counts") +
    scale_x_discrete(breaks = seq(min, max, by)) +
    guides(fill = FALSE)

  return(return(savePlot(plot, output)))
}

#' Get coverage window plot of reads
#'
#' Spanning a region like a transcripts, plot how the reads distribute.
#'
#' If you return this function without assigning it and output is NULL,
#' it will automaticly plot the figure in your session. If output is assigned,
#' no plot will be shown in session.
#' @param coverage a data.table, output of scaledWindowCoverage
#' @param output character string (NULL), if set, saves the plot as pdf or png
#' to path given. If no format is given, is save as pdf.
#' @param scoring character vector (zscore), either of zScore,
#' transcriptNormalized, sum, mean, median, NULL. Set NULL if already scored.
#' @param colors character vector colors to use in plot
#' @param title a character (metaplot) (what is the title of plot?)
#' @param type a character (transcript), what should legends say is
#' the whole region? Transcript, gene, non coding rna etc.
#' @import ggplot2
#' @importFrom data.table copy
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
#' @family coveragePlot
#' @export
#' @examples
#' library(data.table)
#' coverage <- data.table(position = seq(20), score = cumsum(seq(20)))
#' windowCoveragePlot(coverage)
#' # See vignette for a more practical example
#'
windowCoveragePlot <- function(coverage, output = NULL, scoring = "zscore",
                               colors = c('skyblue4', 'orange'),
                               title = "Coverage metaplot",
                               type = "transcript") {
  cov <- copy(coverage)
  if (is.null(cov$feature))
    cov[, feature := rep("meta", nrow(cov))]
  if (is.null(cov$fraction)) {
    cov[, fraction := rep("range", nrow(cov))]
  }
  colors <- colors[seq(length(unique(cov$fraction)))]

  cov$feature  <- factor(cov$feature, levels = unique(cov$feature),
                        labels = unique(cov$feature))
  cov$fraction <- factor(cov$fraction, levels = unique(cov$fraction),
                         labels = unique(cov$fraction))

  coverage_score <- coverageScorings(cov, scoring)

  coverage_score[, `:=` (fraction_min=min(score)), by = list(fraction)]

  plot <- ggplot(data = as.data.frame(coverage_score),
                 aes(x = position, ymax = score, ymin = fraction_min,
                     y = score, colour = as.factor(fraction))) +
    geom_ribbon(stat = "identity", position = "identity",
                aes(fill = as.factor(fraction), alpha = 0.5)) +
    geom_line() +
    theme_bw() + theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    ggtitle(label = title,
            subtitle = paste0("Genes n=", length(unique(coverage$genes)))) +
    xlab(paste("Scaled position in", type)) +
    ylab(paste0(scoring, " over ", type)) +
    theme(legend.position = "none") +
    facet_grid(fraction ~ feature, scales = "free")

  return(savePlot(plot, output))
}

#' Create a heatmap of coverage
#'
#' Coverage rows in heat map is fraction
#' Coverage column in heat map is score, default zscore of counts
#'
#' See vignette for example
#' @inheritParams windowCoveragePlot
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
#' @import ggplot2
#' @family coveragePlot
#' @examples
#' # An ORF
#' grl <- GRangesList(tx1 = GRanges("1", IRanges(1, 6), "+"))
#' # Ribo-seq reads
#' range <- IRanges(c(rep(1, 3), 2, 3, rep(4, 2), 5, 6), width = 1 )
#' reads <- GRanges("1", range, "+")
#' reads$size <- c(rep(28, 5), rep(29, 4)) # read size
#' coverage <- ORFik:::windowPerReadLength(grl, reads = reads, upstream = 0,
#'                                         downstream = 5)
#'
#'
#' ORFik:::coverageHeatMap(coverage)
#'
#' # See vignette for more examples
#'
coverageHeatMap <- function(coverage, output = NULL, scoring = "zscore") {
  coverage$fraction <- factor(coverage$fraction,
                              levels = unique(coverage$fraction),
                              labels = unique(coverage$fraction))

  min <- min(coverage$position)
  max <- max(coverage$position)
  by <- ifelse(length(coverage$position) > 80, 3, 1)

  plot <- ggplot(as.data.frame(coverage) ,
                 aes(x = position, y = fraction, fill = score)) +
    geom_tile()  +
    scale_fill_gradientn(colours = c("white", "yellow2", "yellow3",
                                     "lightblue", "blue", "navy"),
                         name = scoring) +
    xlab("Position relative to start site") +
    ylab("Protected fragment length") +
    scale_x_continuous(breaks = seq(min, max, by)) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +
    theme(text = element_text(size = 12))

  return(savePlot(plot, output))
}

#' Helper function for writing plots to disc
#' @param plot the ggplot to save
#' @param output character string (NULL), if set, saves the plot as pdf or png
#' to path given. If no format is given, is save as pdf.
#' @param width width of output in mm
#' @param height height of output in mm
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
#' @family coveragePlot
savePlot <- function(plot, output = NULL, width = 200, height = 150) {
  if (!is.null(output)) {
    if (is.character(output) && dir.exists(dirname(output))) {
      ext <- tools::file_ext(output)
      if (ext != "pdf" & ext != "png") output <- paste0(output, ".pdf")
      ggsave(output, plot = plot, width = width, height = height, units = "mm",
             dpi = 150, limitsize = FALSE)
    } else {
      stop("output does not name a valid directory")
    }
    return(NULL)
  } else return(plot)
}
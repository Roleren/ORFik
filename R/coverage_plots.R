#' Plot area around TIS for p-shifted reads
#'
#' Usefull to validate p-shifting is correct
#' Can be used for any coverage of region around a point, like TIS, TSS,
#' stop site etc.
#'
#' Remember if you want to change anything like colors, just return the
#' ggplot object, and reassign like: obj + scale_color_brewer() etc.
#' @param hitMap a data.frame/data.table, given from metaWindow
#' (must have columns: position, (score or count) and frame)
#' @param length an integer (29), which length is this for?
#' @param region a character (start), either "start or "stop"
#' @param output character (NULL), if set, saves the plot as pdf or png
#' to path given. If no format is given, is save as pdf.
#' @param type character (canonical CDS), type for plot
#' @param scoring character (Average sum) which scoring did you use ?
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
pSitePlot <- function(hitMap, length = 29, region = "start", output = NULL,
                      type = "canonical CDS", scoring = "Averaged counts") {
  hitMap <- setDT(copy(hitMap))
  if (is.null(hitMap$score)) hitMap[, score := count]
  if (is.null(hitMap$frame)) stop("frame must be included in hitMap for now!")
  plot <- ggplot(hitMap, aes(x = factor(position), y = score,
                             fill = factor(frame))) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = paste("Length", length, "over", region, "of", type)) +
    xlab(paste("\nshift from first", region, "nucleotide [bp]")) +
    ylab(scoring) +
    scale_x_discrete(breaks = xAxisScaler(hitMap$position)) +
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
#'
#' Remember if you want to change anything like colors, just return the
#' ggplot object, and reassign like: obj + scale_color_brewer() etc.
#' @param coverage a data.table, e.g. output of scaledWindowCoverage
#' @param output character string (NULL), if set, saves the plot as pdf or png
#' to path given. If no format is given, is save as pdf.
#' @param scoring character vector (zscore), either of zScore,
#' transcriptNormalized, sum, mean, median, NULL. Set NULL if already scored.
#' @param colors character vector colors to use in plot, will fix automaticly,
#' using binary splits with colors c('skyblue4', 'orange').
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
  cov <- setDT(copy(coverage))
  if (is.null(cov$feature))
    cov[, feature := rep("meta", nrow(cov))]
  if (is.null(cov$fraction)) {
    cov[, fraction := rep("range", nrow(cov))]
  }
  cov$feature  <- factor(cov$feature, levels = unique(cov$feature),
                         labels = unique(cov$feature))
  cov$fraction <- factor(cov$fraction, levels = unique(cov$fraction),
                         labels = unique(cov$fraction))

  coverage_score <- coverageScorings(cov, scoring)

  coverage_score[, `:=` (fraction_min=min(score)), by = fraction]
  nGenes <- getNGenesCoverage(coverage)
  subTitle <- ifelse(any(nGenes > 0), paste0("Genes n=", nGenes), "")
  colors <- matchColors(cov, colors)

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
    ggtitle(label = title, subtitle = subTitle) +
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
#'
#' Remember if you want to change anything like colors, just return the
#' ggplot object, and reassign like: obj + scale_color_brewer() etc.
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
#' ORFik:::coverageHeatMap(coverage)
#'
#' # See vignette for more examples
#'
coverageHeatMap <- function(coverage, output = NULL, scoring = "zscore") {
  coverage$fraction <- factor(coverage$fraction,
                              levels = unique(coverage$fraction),
                              labels = unique(coverage$fraction))

  plot <- ggplot(coverage, aes(x = position, y = fraction, fill = score)) +
    geom_tile()  +
    scale_fill_gradientn(colours = c("white", "yellow2", "yellow3",
                                     "lightblue", "blue", "navy"),
                         name = scoring) +
    xlab("Position relative to start site") +
    ylab("Protected fragment length") +
    scale_x_continuous(breaks = xAxisScaler(coverage$position)) +
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
#' @param dpi (300) dpi of plot
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
#' @family coveragePlot
savePlot <- function(plot, output = NULL, width = 200, height = 150,
                     dpi = 300) {
  if (!is.null(output)) {
    if (is.character(output) && dir.exists(dirname(output))) {
      ext <- tools::file_ext(output)
      if (ext != "pdf" & ext != "png") output <- paste0(output, ".pdf")
      ggsave(output, plot = plot, width = width, height = height, units = "mm",
             dpi = dpi, limitsize = FALSE)
    } else {
      stop("output does not name a valid directory")
    }
    return(NULL)
  } else return(plot)
}

#' Scale x axis correctly
#'
#' Works for all coverage plots, that need 0 position aligning
#'
#' @param covPos a numeric vector of positions in coverage
#' @return a numeric vector from the seq() function, aligned to 0.
xAxisScaler <- function(covPos) {
  pos <- length(covPos)
  min <- min(covPos)
  max <- max(covPos)
  by <- ifelse(pos > 80, ifelse(pos > 150, 6, 3), 1)

  return(seq(min, max, by) - (min %% by))
}

#' Match coloring of coverage plot
#'
#' Check that colors match with the number of fractions.
#' @param coverage a data.table with coverage
#' @param colors a character vector of colors
#' @return number of genes in coverage
matchColors <- function(coverage, colors) {
  nFractions <- length(unique(coverage$fraction))
  nColors <- length(colors)
  if (nColors == 0 || nFractions == 0)
    stop("did not define fraction or colors")

  if (nColors < nFractions) {
    return(rep(colors, nFractions)[seq(nFractions)])
  }
  return(colors[seq(nFractions)])
}

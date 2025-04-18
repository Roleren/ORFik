#' Plot area around TIS as histogram
#'
#' Usefull to validate p-shifting is correct
#' Can be used for any coverage of region around a point, like TIS, TSS,
#' stop site etc.
#'
#' The region is represented as a histogram with different colors for the
#' 3 frames. To make it easy to see patterns in the reads.
#' Remember if you want to change anything like colors, just return the
#' ggplot object, and reassign like: obj + scale_color_brewer() etc.
#' @param hitMap a data.frame/data.table, given from metaWindow
#' (must have columns: position, (score or count) and frame)
#' @param length an integer (29), which read length is this for?
#' @param region a character (start), either "start or "stop"
#' @param output character (NULL), if set, saves the plot as pdf or png
#' to path given. If no format is given, is save as pdf.
#' @param type character (canonical CDS), type for plot
#' @param scoring character, default: (Averaged counts),
#' which scoring did you use ?
#' see ?coverageScorings for info and more alternatives.
#' @param forHeatmap a logical (FALSE), should the plot be part of
#' a heatmap? It will scale it differently. Removing title, x and y labels, and
#' truncate spaces between bars.
#' @param title character, title of plot. Default "auto", will make it:
#' paste("Length", length, "over", region, "of", type).
#' Else set your own (set to NULL to remove all together).
#' @param facet logical, default FALSE. If you input multiple read lengths,
#' specified by fraction column of hitMap, it will split the plots for
#' each read length, putting them under each other. Ignored if forHeatmap
#' is TRUE.
#' @param frameSum logical default FALSE. If TRUE, add an addition plot
#' to the right, sum per frame over all positions per length.
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
#' @importFrom data.table setDT
#' @importFrom cowplot plot_grid
#' @family coveragePlot
#' @export
#' @examples
#' # An ORF
#' grl <- GRangesList(tx1 = GRanges("1", IRanges(1, 6), "+"))
#' # Ribo-seq reads
#' range <- IRanges(c(rep(1, 3), 2, 3, rep(4, 2), 5, 6), width = 1 )
#' reads <- GRanges("1", range, "+")
#' coverage <- coveragePerTiling(grl, reads, TRUE, as.data.table = TRUE,
#'                               withFrames = TRUE)
#' pSitePlot(coverage)
#'
#' # See vignette for more examples
#'
pSitePlot <- function(hitMap, length = unique(hitMap$fraction),
                      region = "start", output = NULL,
                      type = "canonical CDS", scoring = "Averaged counts",
                      forHeatmap = FALSE, title = "auto",
                      facet = FALSE, frameSum = FALSE) {
  hitMap <- setDT(copy(hitMap))
  if (is.null(hitMap$score)) hitMap[, score := count]
  if (is.null(hitMap$frame)) {
    hitMap[, frame := rep("1", nrow(hitMap))]
  }
  plot <- ggplot(hitMap, aes(x = position, y = score,
                             fill = factor(frame))) +
    guides(fill = "none")
  if (forHeatmap) {
    plot <- plot + ggtitle("") +
      geom_bar(stat = "identity", width = 1) +
      xlab("") +
      ylab("") +
      theme(axis.ticks = element_blank(),
            axis.text = element_blank(),
            plot.margin = unit(c(0.1,0.3,-0.5,0.8), "cm")) +
      theme(panel.background=element_rect(fill="white", colour="gray")) +
      scale_fill_grey()
  } else {
    if (!is.null(title)) {
      if (title == "auto") {
        title <- paste("Length", length, "over", region, "of", type)
      }
    }

    plot <- plot +
      geom_bar(stat = "identity") +
      labs(title = title) +
      xlab(paste("\nshift from first", region, "nucleotide [bp]")) +
      ylab(prettyScoring(scoring)) +
      scale_x_continuous(breaks = xAxisScaler(hitMap$position)) +
      scale_y_continuous(n.breaks = 3)
    if (facet) {
      plot <- plot +
        facet_wrap(~ fraction, ncol = 1, scales = "free_y", strip.position = "right")
    }
    if (frameSum) {
      hitMap2 <- setDT(copy(hitMap))
      hitMap2[, position := frame]
      plot2 <- pSitePlot(hitMap2, facet = facet, title = "Sum per Frame")
      plot2 <- plot2 + ylab(label = "") + xlab("Frame")
      plot <- cowplot::plot_grid(plot, plot2, ncol = 2,
                                 rel_widths = c(3,1), align = "v")
    }
  }

  return(savePlot(plot, output))
}

#' Get meta coverage plot of reads
#'
#' Spanning a region like a transcripts, plot how the reads distribute.
#'
#' If coverage has a column called feature, this can be used to subdivide the
#' meta coverage into parts as (5' UTRs, cds, 3' UTRs) These are the columns
#' in the plot.
#' The fraction column divide sequence libraries. Like ribo-seq and rna-seq.
#' These are the rows of the plot.
#' If you return this function without assigning it and output is NULL,
#' it will automaticly plot the figure in your session. If output is assigned,
#' no plot will be shown in session. NULL is returned and object is saved to
#' output.
#'
#' Colors:
#' Remember if you want to change anything like colors, just return the
#' ggplot object, and reassign like: obj + scale_color_brewer() etc.
#' @param coverage a data.table, e.g. output of scaledWindowCoverage
#' @param output character string (NULL), if set, saves the plot as pdf or png
#' to path given. If no format is given, is save as pdf.
#' @param scoring character vector, default "zscore", either of zscore,
#' transcriptNormalized, sum, mean, median, .. or  NULL. Set NULL if already scored.
#' see ?coverageScorings for info and more alternatives.
#' @param colors character vector colors to use in plot, will fix automaticly,
#' using binary splits with colors c('skyblue4', 'orange').
#' @param title a character (metaplot) (what is the title of plot?)
#' @param type a character (transcripts), what should legends say is
#' the whole region? Transcripts, genes, non coding rnas etc.
#' @param scaleEqual a logical (FALSE), should all fractions (rows), have same
#'  max value, for easy comparison of max values if needed.
#' @param setMinToZero a logical (FALSE), should minimum y-value be 0 (TRUE).
#' With FALSE minimum value is minimum score at any position. This parameter
#' overrides scaleEqual.
#' @import ggplot2
#' @importFrom data.table copy
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
#' @family coveragePlot
#' @export
#' @examples
#' library(data.table)
#' coverage <- data.table(position = seq(20),
#'                        score = sample(seq(20), 20, replace = TRUE))
#' windowCoveragePlot(coverage)
#'
#' #Multiple plots in one frame:
#' coverage2 <- copy(coverage)
#' coverage$fraction <- "Ribo-seq"
#' coverage2$fraction <- "RNA-seq"
#' dt <- rbindlist(list(coverage, coverage2))
#' windowCoveragePlot(dt, scoring = "log10sum")
#'
#' # See vignette for a more practical example
#'
windowCoveragePlot <- function(coverage, output = NULL, scoring = "zscore",
                               colors = c('skyblue4', 'orange'),
                               title = "Coverage metaplot",
                               type = "transcripts", scaleEqual = FALSE,
                               setMinToZero = FALSE) {
  cov <- setDT(copy(coverage))
  if (is.null(cov$feature))
    cov[, feature := rep("meta", nrow(cov))]
  if (is.null(cov$fraction)) {
    cov[, fraction := rep("range", nrow(cov))]
  }
  if (!is(cov$feature, "factor")) {
    cov$feature  <- factor(cov$feature, levels = unique(cov$feature),
                           labels = unique(cov$feature))
  }
  cov$fraction <- factor(cov$fraction, levels = unique(cov$fraction),
                         labels = unique(cov$fraction))

  coverage_score <- coverageScorings(cov, scoring)

  # Decide y-min value
  if (setMinToZero) {
    coverage_score[, `:=` (fraction_min = 0)]
  } else coverage_score[, `:=` (fraction_min=min(score)), by = fraction]
  if (scaleEqual) coverage_score[, `:=` (fraction_min = min(fraction_min))]

  nGenes <- getNGenesCoverage(coverage)
  subTitle <- ifelse(any(nGenes > 0), paste0("Genes n=", nGenes), "")
  colors <- matchColors(cov, colors)
  # Y-axis facet name, scaling
  max_name <- max(nchar(unique(as.character(coverage_score$fraction))))
  fraction_text <- 12 - min(6, (max_name / 40) *12)

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
    ylab(paste0(prettyScoring(scoring), " over ", type)) +
    theme(legend.position = "none") +
    facet_grid(fraction ~ feature, scales = ifelse(scaleEqual, "free_x", "free")) +
    theme(strip.text.y = element_text(size = fraction_text),
          strip.text.x = element_text(size = 12))

  return(savePlot(plot, output))
}

#' Helper function for writing plots to disc
#' @param plot the ggplot to save
#' @param output character string (NULL), if set, saves the plot as pdf or png
#' to path given. If no format is given,
#' is save as specified by plot.ext argument.
#' @param width width of output in mm
#' @param height height of output in mm
#' @param dpi (300) dpi of plot
#' @param plot.ext character, default: ".pdf". Alternatives: ".png" or ".jpg".
#' @param limitsize logical, default FALSE. If TRUE, activate ggplot max size restriction.
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
#' @family coveragePlot
#' @keywords internal
savePlot <- function(plot, output = NULL, width = 200, height = 150, plot.ext = ".pdf",
                     dpi = 300, limitsize = FALSE) {
  if (!is.null(output)) {
    if (is.character(output) && dir.exists(dirname(output))) {
      ggsave(image_path_format_append(output, plot.ext), plot = plot, width = width,
             height = height, units = "mm", dpi = dpi, limitsize = limitsize)
    } else {
      stop("output does not name a valid directory")
    }
    return(NULL)
  }
  return(plot)
}

image_path_format_append <- function(output, plot.ext, all_formats = c("pdf", "png", "jpg")) {
  ext <- tools::file_ext(output)
  if (!(ext %in% all_formats)) output <- paste0(output, plot.ext)
  return(output)
}

#' Scale x axis correctly
#'
#' Works for all coverage plots, that need 0 position aligning
#'
#' It basicly bins the x axis on floor(length of x axis / 20) or
#' 1 if x < 20
#' @param covPos a numeric vector of positions in coverage
#' @return a numeric vector from the seq() function, aligned to 0.
#' @keywords internal
xAxisScaler <- function(covPos) {
  pos <- length(unique(covPos))
  min <- min(covPos)
  max <- max(covPos)
  by <- max(floor(pos / 20), 1)

  return(seq(min, max, by) - (min %% by))
}

#' Scale y axis correctly
#'
#' Works for all coverage plots.
#'
#' @param covPos a levels object from a factor of y axis
#' @param increments.y increments of y axis, default "auto".
#' Or a numeric value < max position & > min position.
#' @return a character vector from the seq() function, aligned to 0.
#' @keywords internal
yAxisScaler <- function(covPos, increments.y = "auto") {
  covPos <- as.integer(covPos)
  pos <- length(covPos)
  min <- min(covPos)
  max <- max(covPos)
  if (increments.y == "auto") {
    by <- ifelse(pos > 25, ifelse(pos > 50, ifelse(pos > 70, ifelse(pos > 120,
                                                   ifelse(pos > 300, 100, 50), 20), 10), 4), 1)
  } else if (is.numeric(increments.y)) {
    by <- increments.y
    if ((by < min) | (by > max))
      stop("increments.y must be > min pos and < max pos of y-axis")
  } else stop("increments.y must be auto or a numeric value")

  return(as.character(seq.int(min, max, by)))
}

#' Prettify scoring name
#' @param scoring a character (the scoring)
#' @return a new scoring name or the same if pretty
#' @keywords internal
prettyScoring <- function(scoring) {
  if (is.null(scoring)) return("raw")
  if (scoring == "log2sum") {
    scoring <- "log2(sum)"
  } else if (scoring == "log10sum") {
    scoring <- "log10(sum)"
  } else if (scoring == "transcriptNormalized") {
    scoring <- "Transcript Normalized"
  }
  return(scoring)
}

#' Match coloring of coverage plot
#'
#' Check that colors match with the number of fractions.
#' @param coverage a data.table with coverage
#' @param colors a character vector of colors
#' @return number of genes in coverage
#' @keywords internal
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

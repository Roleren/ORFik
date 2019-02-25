#' Plot area around TIS for p-shifted reads
#'
#' Usefull to validate p-shifting is correct
#' @param hitMap a data.frame, given from metaWindow
#' @param length an integer (29), which length is this for?
#' @param region a character (start), either "start or "stop".
#' @return a ggplot
pSitePlot <- function(hitMap, length = 29, region = "start") {
  plot <- ggplot(hitMap, aes(x = factor(position), y = counts,
                          fill = factor(frame))) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
    labs(title = paste("Length", length, "over", region, "of canonical CDS")) +
    xlab(paste("\nshift from first", region, "nucleotide [bp]")) +
    ylab("Averaged counts") +
    guides(fill = FALSE)

  return(plot)
}

#' Get coverage window plot of reads
#'
#' Spanning a region like a transcripts, plot how the reads distribute.
#'
#' If you return this function without assigning it and output is NULL,
#' it will automaticly plot the figure in your session. If output is assigned,
#' no plot will be shown in session.
#' @param coverage a data.table, output of scaledWindowCoverage
#' @param output character string (NULL), if set, saves the plot as pdf
#' to path given.
#' @param scoring character vector (zscore), either of zScore,
#' transcriptNormalized, sum, mean
#' @param colors character vector colors to use in plot
#' @param title a character (metaplot) (what is the title of plot?)
#' @param type a character (transcript), what should legends say is
#' the whole region? Transcript, gene, non coding rna etc.
#' @import ggplot2
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
windowCoveragePlot <- function(coverage, output = NULL, scoring = "zscore",
                               colors = c('skyblue4', 'orange'),
                               title = "Coverage metaplot",
                               type = "transcript") {
  if(is.null(coverage$feature))
    coverage[, feature := rep("meta", nrow(coverage))]
  if(is.null(coverage$fraction)) {
    coverage[, fraction := rep("range", nrow(coverage))]
    colors <- colors[1]
  }
  coverage$feature  <- factor(coverage$feature,
                              levels = unique(coverage$feature),
                              labels = unique(coverage$feature))
  coverage$fraction <- factor(coverage$fraction,
                              levels = unique(coverage$fraction),
                              labels = unique(coverage$fraction))

  coverage_score <- coverageScorings(coverage, scoring)

  coverage_score[, `:=` (fraction_min=min(score)), by = list(fraction)]

  plot <- ggplot(data=as.data.frame(coverage_score),
                 aes(x=position, ymax=score, ymin=fraction_min,
                     y=score, colour = as.factor(fraction))) +
    geom_ribbon(stat="identity", position = "identity",
                aes(fill= as.factor(fraction), alpha=0.5)) +
    geom_line() +
    theme_bw() +   theme(panel.grid.major = element_blank(),
                         panel.grid.minor = element_blank()) +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values=colors) +
    scale_color_manual(values=colors) +
    ggtitle(label = title,
            subtitle = paste0("Genes n=",
                              length(unique(coverage_score$genes))) ) +
    xlab(paste("Scaled position in", type)) +
    ylab(paste0(scoring, " over ", type)) +
    theme(legend.position="none") +
    facet_grid(fraction ~ feature, scales = "free")

  if (!is.null(output)) {
    if(is.character(output) && dir.exists(dirname(outName))) {
      if (tools::file_ext(output) != "pdf") output <- paste0(output, ".pdf")
      ggsave(output, plot = plot, width = 200, height=150, units = "mm",
             dpi = 100, limitsize = FALSE)
    } else {
      stop("output does not name a valid directory")
    }
    return(NULL)
  } else return(plot)
}

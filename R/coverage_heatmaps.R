#' Create a heatmap of coverage
#'
#' Creates a ggplot representing a heatmap of coverage:\cr
#' \itemize{
#'  \item{Rows            : }{Position in region}
#'  \item{Columns         : }{Read length}
#'  \item{Index intensity : }{(color) coverage scoring per index.}
#' }
#' Coverage rows in heat map is fraction, usually fractions is divided into
#' unique read lengths (standard Illumina is 76 unique widths, with some
#' minimum cutoff like 15.)
#' Coverage column in heat map is score, default zscore of counts. These are
#' the relative positions you are plotting to. Like +/- relative to TIS or TSS.
#'
#' Colors:
#' Remember if you want to change anything like colors, just return the
#' ggplot object, and reassign like: obj + scale_color_brewer() etc.
#' Standard colors are:\cr
#' \itemize{
#'  \item{0 reads in whole readlength :}{gray}
#'  \item{few reads in position       :}{white}
#'  \item{medium reads in position    :}{yellow}
#'  \item{many reads in position      :}{dark blue}
#' }
#'
#' @inheritParams windowCoveragePlot
#' @param legendPos a character, Default "right". Where should the fill legend
#' be ? ("top", "bottom", "right", "left")
#' @param addFracPlot Add margin histogram plot on top of heatmap with
#'  fractions per positions
#' @param xlab the x-axis label, default "Position relative to start site"
#' @param ylab the y-axis label, default "Protected fragment length"
#' @param colors character vector, default: "default", this gives you:
#'  c("white", "yellow2", "yellow3", "lightblue", "blue", "navy"),
#'  do "high" for more high contrasts, or specify your own colors.
#' @param scoring character vector, default "zscore",
#' Which scoring did you use to create? either of zscore,
#' transcriptNormalized, sum, mean, median, ..
#' see ?coverageScorings for info and more alternatives.
#' @param title a character, default NULL (no title),
#' what is the top title of plot?
#' @param gradient.max numeric, defualt: max(coverage$score). What data value
#' should the top color be ? Good to use if you want to compare 2 samples, with the
#' same color intensity, in that case set this value to the max score of the
#' 2 coverage tables.
#' @inheritParams yAxisScaler
#' @return a ggplot object of the coverage plot, NULL if output is set,
#' then the plot will only be saved to location.
#' @import ggplot2
#' @importFrom cowplot plot_grid
#' @family heatmaps
#' @family coveragePlot
#' @export
#' @examples
#' # An ORF
#' grl <- GRangesList(tx1 = GRanges("1", IRanges(1, 6), "+"))
#' # Ribo-seq reads
#' range <- IRanges(c(rep(1, 3), 2, 3, rep(4, 2), 5, 6), width = 1 )
#' reads <- GRanges("1", range, "+")
#' reads$size <- c(rep(28, 5), rep(29, 4)) # read size
#' coverage <- windowPerReadLength(grl, reads = reads, upstream = 0,
#'                                 downstream = 5)
#'
#' coverageHeatMap(coverage)
#'
#' # With top sum bar
#' coverageHeatMap(coverage, addFracPlot = TRUE)
#' # See vignette for more examples
#'
coverageHeatMap <- function(coverage, output = NULL, scoring = "zscore",
                            legendPos = "right", addFracPlot = FALSE,
                            xlab = "Position relative to start site",
                            ylab = "Protected fragment length",
                            colors = "default", title = NULL,
                            increments.y = "auto",
                            gradient.max = max(coverage$score)) {
  coverage$fraction <- factor(coverage$fraction,
                              levels = unique(coverage$fraction),
                              labels = unique(coverage$fraction))
  if (colors == "default") {
    colors = c("white", "yellow2", "yellow3", "lightblue", "blue", "navy")
  } else if (colors == "high") {
    colors <- c("white", "yellow2", "yellow3", "lightblue", "blue", "blue",
                "blue", "navy", "black")
  }

  plot <- ggplot(coverage, aes(x = position, y = fraction, fill = score)) +
    geom_tile()  +
    scale_fill_gradientn(colours = colors,
                         limits = c(min(coverage$score), gradient.max),
                         name = prettyScoring(scoring)) +
    xlab(xlab) +
    ylab(ylab) +
    scale_x_continuous(breaks = xAxisScaler(coverage$position)) +
    theme_bw() + theme(panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank()) +
    scale_y_discrete(breaks = yAxisScaler(levels(coverage$fraction),
                                          increments.y)) +
    theme(legend.position = legendPos)

  if (addFracPlot) {
    if (legendPos != "none") legendPos <- "bottom"
    plot2 <- pSitePlot(coverage, forHeatmap = TRUE) + ggtitle(title)
    plot <- plot_grid(plot2,
                      plot + theme(legend.position = legendPos,
                                   plot.margin = unit(c(0,0.3,0,0.8), "cm")),
                      ncol = 1, rel_heights = c(1,4), align = "v")
  } else plot <- plot + ggtitle(title)
  return(savePlot(plot, output))
}

#' Create coverage heatmaps of specified region
#'
#' Simplified input space for easier abstraction of coverage heatmaps\cr
#' Pick your region and plot \cr
#' Input CAGE file if you use TSS and want improved 5' annotation.
#'
#' @inheritParams heatMapL
#' @param region a character, default "TIS", can be any combination of the
#'  set: c("TSS", "TIS", "TTS"), which are: Transcription start site
#'  (5' end of mrna), Translation initation site (5' end of CDS),
#'  Translation termination site (3' end of CDS)
#' @param outdir a character path, default: "default", saves to:
#' \code{paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/")}, a created
#' folder within the ORFik experiment data folder for plots. Change if you
#' want custom location.
#' @param cage a character path to library file or a \code{\link{GRanges}},
#' \code{\link{GAlignments}} preloaded file of CAGE data. Only used if
#' "TSS" is defined as region, to redefine 5' leaders.
#' @return invisible(NULL), plots are saved
#' @family heatmaps
#' @export
#' @examples
#' # Toy example, will not give logical output, but shows how it works
#' df <- ORFik.template.experiment()[3,] # Only third library
#' #heatMapRegion(df, "TIS", outdir = "default")
#' #
#' # Do also TSS, add cage for specific TSS
#' # heatMapRegion(df, c("TSS", "TIS"), cage = "path/to/cage.bed")
#'
#' # Do on pshifted reads instead of original files
#' remove.experiments(df) # Remove loaded experiment first
#' # heatMapRegion(df, "TIS", type = "pshifted")
heatMapRegion <- function(df, region = "TIS", outdir = "default",
                          scores = c("transcriptNormalized", "sum"),
                          type = "ofst", cage = NULL, format = ".png",
                          acceptedLengths = 21:75, upstream = c(50, 30),
                          downstream = c(29, 69),
                          shifting = c("5prime", "3prime")) {

  if (outdir == "default") outdir <- paste0(dirname(df$filepath[1]), "/QC_STATS/heatmaps/")
  if (!(any(region %in% c("TIS", "TSS", "TTS")))) stop("region must be either TSS, TIS or TTS")
  dir.create(outdir, showWarnings = FALSE,
             recursive = TRUE)
  message(paste0("Plot save location:\n", outdir))
  txdb <- loadTxdb(df)
  if ("TIS" %in% region) {
    message("TIS")
    txNames <- filterTranscripts(txdb, 51, 70, 0)
    center <- loadRegion(txdb, "cds")[txNames]
    mrna <- loadRegion(txdb, "mrna")[txNames]
    heatMapL(center, mrna, df, outdir, scores = scores, upstream, downstream,
             addFracPlot = TRUE, location = "TIS", shifting = shifting,
             skip.last = FALSE, acceptedLengths = acceptedLengths, type = type)
  }
  if ("TSS" %in% region) {
    message("TSS")
    txNames <- filterTranscripts(txdb, 70, 0, 0)
    center <- loadRegion(txdb, "leaders")[txNames]
    mrna <- loadRegion(txdb, "mrna")
    len <- startSites(center, keep.names = TRUE, is.sorted = TRUE) > 51
    center <- center[len]
    mrna <- mrna[names(center)]
    if (!is.null(cage)) {
      message("Using cage file to update TSS")
      center <- reassignTSSbyCage(center, cage)
      mrna <- downstreamFromPerGroup(mrna, startSites(center, is.sorted = TRUE))
    }
    mrna <- extendLeaders(mrna, 51)
    heatMapL(center, mrna, df, outdir, scores = scores, upstream, downstream,
             addFracPlot = TRUE, location = "TSS", shifting = shifting,
             skip.last = FALSE, acceptedLengths = acceptedLengths, type = type)
  }
  if ("TTS" %in% region) {
    message("TTS")
    txNames <- filterTranscripts(txdb, 0, 51, 70)
    center <- loadRegion(txdb, "trailers")[txNames]
    mrna <- loadRegion(txdb, "mrna")[txNames]
    heatMapL(center, mrna, df, outdir, scores = scores, upstream, downstream,
             addFracPlot = TRUE, location = "TTS", shifting = shifting,
             skip.last = FALSE, acceptedLengths = acceptedLengths, type = type)
  }
  return(invisible(NULL))
}

#' Coverage heatmap of multiple libraries
#'
#' @inheritParams heatMap_single
#' @param df an ORFik \code{\link{experiment}}
#' @param type character, default: "ofst". Type of library:
#' either "default", usually bam format (the one you gave to experiment),
#' "pshifted" pshifted reads, "ofst", "bed", "bedo" optimized bed, or "wig"
#' @param outdir a character path to directory to save plot, will be named
#' from ORFik experiment columns
#' @param plot.together logical (default: FALSE), plot all in 1 plot (if TRUE)
#' @param shifting a character, default c("5prime", "3prime"), can also be
#' either or NULL (no shifting of reads)
#' @param format a character, default ".png", alternative ".pdf"
#' @param upstream 1 or 2 integers, default c(50, 30), how long upstream from 0
#' should window extend (first index is 5' end extension, second is 3' end extension).
#' If only 1 shifting, only 1 value should be given, if two are given will use first.
#' @param downstream 1 or 2 integers, default c(29, 69), how long upstream from 0
#' should window extend (first index is 5' end extension, second is 3' end extension).
#' If only 1 shifting, only 1 value should be given, if two are given will use first.
#' @param scores character vector, default c("transcriptNormalized", "sum"),
#' either of zscore, transcriptNormalized, sum, mean, median, ..
#' see ?coverageScorings for info and more alternatives.
#' @importFrom gridExtra grid.arrange
#' @return invisible(NULL), plots are saved
#' @family heatmaps
heatMapL <- function(region, tx, df, outdir, scores = "sum", upstream, downstream,
                     zeroPosition = upstream, acceptedLengths = NULL, type = "ofst",
                     legendPos = "right", colors = "default", addFracPlot = TRUE,
                     location = "TIS", shifting = NULL, skip.last = FALSE, format = ".png",
                     plot.together = TRUE, title = TRUE) {


  up <- upstream; down <- downstream
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)

  for (df in dfl) {
    heatmapList <- list()
    varNames <- bamVarName(df)
    outputLibs(df, region, type = type)
    for (i in varNames) { # For each stage
      for (score in scores) {
        for (s in seq_along(shifting)) {
          shift <- shifting[s]
          if (length(upstream) > 1) {
            up <- upstream[s]
          }
          if (length(downstream) > 1) {
            down <- downstream[s]
          }
          if (length(zeroPosition) > 1) {
            zero <- zeroPosition[s]
          }
          print(paste(i, shift, score))
          out <- paste0(outdir, df@experiment,"_hm_", location, "_",i , "_")
          out <- ifelse(!is.null(shifting),
                        paste0(out, shift, "_", score, format),
                        paste0(out, score, format))

          plot <- heatMap_single(region, tx, reads = get(i), outdir = out,
                                 shifting = shift, scores = score, upstream = up, downstream = down,
                                 zeroPosition = zero, acceptedLengths = acceptedLengths,
                                 legendPos = legendPos, colors = colors, addFracPlot = addFracPlot,
                                 location = location, skip.last = skip.last,
                                 title = ifelse(title, paste(i, shift), NULL))
          heatmapList <- c(heatmapList, list(plot))
        }
      }
    }
    # Per experiment plot together
    if (plot.together) {
      ncols <- max(1, length(shifting))
      final <- gridExtra::grid.arrange(grobs = heatmapList, ncol = ncols)
      ggsave(paste0(outdir, df@experiment, "_hm_combined_", location, format),
             plot = final, width = 5*ncols, height = ceiling(5.5*(length(heatmapList) / ncols)),
             limitsize = FALSE)
    }
  }
  return(invisible(NULL))
}

#' Coverage heatmap of single libraries
#' @inheritParams coverageHeatMap
#' @inheritParams windowPerReadLength
#' @param region #' a \code{\link{GRangesList}} object of region,
#'  usually either leaders, cds', 3' utrs or ORFs, start region, stop regions etc.
#'  This is the region that will be mapped in heatmap
#' @param upstream an integer, relative region to get upstream from.
#' @param downstream an integer, relative region to get downstream from
#' @param location a character, default "start site", will make xlabel of heatmap be
#' Position relative to "start site" or alternative given.
#' @param scores character vector, default "sum",
#' either of zscore, transcriptNormalized, sum, mean, median, ..
#' see ?coverageScorings for info and more alternatives.
#' @param skip.last skip top(highest) read length, default FALSE
#' @param shifting a character, default NULL (no shifting), can also be
#' either of c("5prime", "3prime")
#' @param returnCoverage logical, default: FALSE, return coverage, if FALSE
#' returns plot instead.
#' @param outdir a character path to save file as: not just directory,
#' but full name.
#' @return ggplot2 grob (default), data.table (if returnCoverage is TRUE)
#' @importFrom gridExtra grid.arrange
#' @family heatmaps
#' @export
heatMap_single <- function(region, tx, reads, outdir,
                           scores = "sum", upstream, downstream,  zeroPosition = upstream,
                           returnCoverage = FALSE, acceptedLengths = NULL, legendPos = "right",
                           colors = "default", addFracPlot = TRUE, location = "start site", shifting = NULL,
                           skip.last = FALSE, title = NULL) {
  if (length(scores) != 1) stop("scores must exactly only 1 score type")
  if (!is.null(shifting)) {
    reads <- convertToOneBasedRanges(reads, method = shifting, addSizeColumn = TRUE)
  }
  if (skip.last) {
    all_lengths <- sort(unique(readWidths(reads)))
    if (!is.null(acceptedLengths))
      all_lengths <- all_lengths[all_lengths %in% acceptedLengths]
    acceptedLengths <- all_lengths[-c((length(all_lengths)-0):length(all_lengths))]
  }
  dt <- windowPerReadLength(region, tx, reads, upstream = upstream, downstream = downstream,
                            zeroPosition = zeroPosition, scoring = scores,
                            acceptedLengths = acceptedLengths)

  plot <- coverageHeatMap(coverage = dt, scoring = scores, addFracPlot = addFracPlot,
                          xlab = paste0("Position relative to ", location), colors = colors,
                          legendPos = legendPos, title = title)

  ggsave(filename = outdir, plot = plot, width = 350, height = 180,
         units = "mm", dpi = 300, limitsize = FALSE)
  # return coverage or the plot
  if (returnCoverage) return(dt)
  return(plot)
}

#' Create a heatmap of coverage
#'
#' Creates a ggplot representing a heatmap of coverage:\cr
#' \itemize{
#'  \item{Rows : Position in region}
#'  \item{Columns : Read length}
#'  \item{Index intensity : (color) coverage scoring per index.}
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
#'  \item{0 reads in whole readlength : gray}
#'  \item{few reads in position : white}
#'  \item{medium reads in position : yellow}
#'  \item{many reads in position : dark blue}
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
  if (colors[1] == "default") {
    colors = c("white", "yellow2", "yellow3", "lightblue", "blue", "navy")
  } else if (colors[1] == "high") {
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
#' Pick your transcript region and plot directly \cr
#' Input CAGE file if you use TSS and want improved 5' annotation.
#'
#' @inheritParams heatMapL
#' @inheritParams filterTranscripts
#' @param region a character, default "TIS". The centering point for the heatmap
#' (what is position 0, beween -50 and 50 etc), can be any combination of the
#'  set: c("TSS", "TIS", "TTS", "TES"), which are:
#'  - Transcription start site (5' end of mrna)\cr
#'  - Translation initation site (5' end of CDS)\cr
#'  - Translation termination site (5' end of 3' UTRs)\cr
#'  - Transcription end site (3' end of 3' UTRs)\cr
#' @param outdir a character path, default: "default", saves to:
#' \code{file.path(QCfolder(df), "heatmaps/")}, a created
#' folder within the ORFik experiment data folder for plots. Change if you
#' want custom location.
#' @param cage a character path to library file or a \code{\link{GRanges}},
#' \code{\link{GAlignments}} preloaded file of CAGE data. Only used if
#' "TSS" is defined as region, to redefine 5' leaders.
#' @param longestPerGene logical, default TRUE. Use only longest transcript
#' isoform per gene. This will speed up your computation.
#' @return invisible(NULL), plots are saved
#' @family heatmaps
#' @export
#' @examples
#' # Toy example, will not give logical output, but shows how it works
#' df <- ORFik.template.experiment()[9:10,] # Subset to 2 Ribo-seq libs
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
                          type = "ofst", cage = NULL, plot.ext = ".pdf",
                          acceptedLengths = 21:75, upstream = c(50, 30),
                          downstream = c(29, 69),
                          shifting = c("5prime", "3prime"),
                          longestPerGene = TRUE, colors = "default",
                          scale_x = 5.5, scale_y = 15.5,
                          gradient.max = "default",
                          BPPARAM = BiocParallel::SerialParam()) {

  if (outdir == "default") outdir <- file.path(QCfolder(df), "heatmaps/")
  if (!(any(region %in% c("TIS", "TSS", "TTS", "TES"))))
    stop("region must be either TSS, TIS, TTS or TES")
  dir.create(outdir, showWarnings = FALSE,
             recursive = TRUE)
  message(paste0("Plot save location:\n", outdir))
  txdb <- loadTxdb(df)
  if ("TIS" %in% region) {
    message("TIS")
    minLeaderLength <- max(upstream) + 1
    minCDSLength <- max(downstream) + 1
    txNames <- filterTranscripts(txdb, minLeaderLength, minCDSLength, 0,
                                 longestPerGene = longestPerGene)
    center <- loadRegion(txdb, "cds", names.keep = txNames)
    mrna <- loadRegion(txdb, "mrna", names.keep = txNames)
    heatMapL(center, mrna, df, outdir, scores = scores, upstream, downstream,
             addFracPlot = TRUE, location = "TIS", shifting = shifting,
             skip.last = FALSE, acceptedLengths = acceptedLengths, type = type,
             plot.ext = plot.ext, colors = colors,
             scale_x = scale_x, scale_y = scale_y,
             gradient.max = gradient.max, BPPARAM = BPPARAM)
  }
  if ("TSS" %in% region) {
    message("TSS")
    minLeaderLength <- max(downstream) + 1
    minOutsideTXLength <- max(upstream) + 1
    txNames <- filterTranscripts(txdb, minLeaderLength, 0, 0,
                                 longestPerGene = longestPerGene)
    center <- loadRegion(txdb, "leaders", names.keep = txNames)
    mrna <- loadRegion(txdb, "mrna")
    len <- startSites(center, keep.names = TRUE, is.sorted = TRUE) > minOutsideTXLength
    center <- center[len]
    mrna <- mrna[names(center)]
    if (!is.null(cage)) {
      message("Using cage file to update TSS")
      center <- reassignTSSbyCage(center, cage)
      mrna <- downstreamFromPerGroup(mrna, startSites(center, is.sorted = TRUE))
    }
    mrna <- extendLeaders(mrna, minOutsideTXLength)
    heatMapL(center, mrna, df, outdir, scores = scores, upstream, downstream,
             addFracPlot = TRUE, location = "TSS", shifting = shifting,
             skip.last = FALSE, acceptedLengths = acceptedLengths, type = type,
             plot.ext = plot.ext, colors = colors,
             scale_x = scale_x, scale_y = scale_y,
             gradient.max = gradient.max, BPPARAM = BPPARAM)
  }
  if ("TTS" %in% region) {
    message("TTS")
    minCDSLength <- max(upstream) + 1
    minTrailerLength <- max(downstream) + 1
    txNames <- filterTranscripts(txdb, 0, minCDSLength, minTrailerLength,
                                 longestPerGene = longestPerGene)
    center <- loadRegion(txdb, "trailers", names.keep = txNames)
    mrna <- loadRegion(txdb, "mrna", names.keep = txNames)
    heatMapL(center, mrna, df, outdir, scores = scores, upstream, downstream,
             addFracPlot = TRUE, location = "TTS", shifting = shifting,
             skip.last = FALSE, acceptedLengths = acceptedLengths, type = type,
             plot.ext = plot.ext, colors = colors,
             scale_x = scale_x, scale_y = scale_y,
             gradient.max = gradient.max, BPPARAM = BPPARAM)
  }
  # Transcription End site
  if ("TES" %in% region) {
    message("TES")
    minTrailerLength <- max(upstream) + 1
    minOutsideTXLength <- max(downstream) + 1
    txNames <- filterTranscripts(txdb, 0, 0, minTrailerLength, longestPerGene = longestPerGene)
    center <- loadRegion(txdb, "trailer", names.keep = txNames)
    center <- stopSites(center, asGR = TRUE,
                        keep.names = TRUE, is.sorted = TRUE)
    center <- groupGRangesBy(center)
    mrna <- loadRegion(txdb, "mrna")
    len <- stopSites(center, keep.names = TRUE, is.sorted = TRUE) >
      minOutsideTXLength
    center <- center[len]
    mrna <- mrna[names(center)]
    mrna <- extendTrailers(mrna, minOutsideTXLength)
    heatMapL(center, mrna, df, outdir, scores = scores,
             upstream, downstream, zeroPosition = upstream, addFracPlot = TRUE, location = "TES",
             shifting = shifting, skip.last = FALSE, acceptedLengths = acceptedLengths,
             type = type, plot.ext = plot.ext,
             scale_x = scale_x, scale_y = scale_y,
             gradient.max = gradient.max, BPPARAM = BPPARAM)
  }
  return(invisible(NULL))
}

#' Coverage heatmap of multiple libraries
#'
#' @inheritParams heatMap_single
#' @inheritParams outputLibs
#' @param type character, default: "ofst". Type of library:
#' either "default", usually bam format (the one you gave to experiment),
#' "pshifted" pshifted reads, "ofst", "bed", "bedo" optimized bed, or "wig"
#' @param outdir a character path to directory to save plot, will be named
#' from ORFik experiment columns
#' @param plot.together logical (default: FALSE), plot all in 1 plot (if TRUE)
#' @param shifting a character, default c("5prime", "3prime"), can also be
#' NULL (no shifting of reads). If NULL, will use first index of 'upstream'
#' and 'downstream' argument.
#' @param plot.ext a character, default ".pdf", alternative ".png"
#' @param upstream 1 or 2 integers, default c(50, 30), how long upstream from 0
#' should window extend (first index is 5' end extension, second is 3' end extension).
#' If only 1 shifting, only 1 value should be given, if two are given will use first.
#' @param downstream 1 or 2 integers, default c(29, 69), how long upstream from 0
#' should window extend (first index is 5' end extension, second is 3' end extension).
#' If only 1 shifting, only 1 value should be given, if two are given will use first.
#' @param scores character vector, default \code{c("transcriptNormalized", "sum")},
#' either of zscore, transcriptNormalized, sum, mean, median, ..
#' see ?coverageScorings for info and more alternatives.
#' @param scale_x numeric, how should the width of the single plots be scaled,
#' bigger the number, the bigger the plot
#' @param scale_y numeric, how should the height of the plots be scaled,
#' bigger the number, the bigger the plot
#' @param BPPARAM a core param, default: single thread: \code{BiocParallel::SerialParam()}.
#'  Set to \code{BiocParallel::bpparam()} to use multicore. Be aware, this uses a lot of
#'  extra ram (40GB+) for larger human samples!
#' @importFrom gridExtra grid.arrange
#' @return invisible(NULL), plots are saved
#' @family heatmaps
#' @keywords internal
heatMapL <- function(region, tx, df, outdir, scores = "sum", upstream, downstream,
                     zeroPosition = upstream, acceptedLengths = NULL, type = "ofst",
                     legendPos = "right", colors = "default", addFracPlot = TRUE,
                     location = "TIS", shifting = NULL, skip.last = FALSE, plot.ext = ".pdf",
                     plot.together = TRUE, title = TRUE,
                     scale_x = 5.5, scale_y = 15.5, gradient.max = "default",
                     BPPARAM = BiocParallel::SerialParam()) {
  up <- upstream; down <- downstream
  dfl <- df
  if (!is(dfl, "list")) dfl <- list(dfl)
  if (is.null(shifting)) shifting <- "NULL"

  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, chrStyle = region, type = type)
    heatmapList <- bplapply(varNames,
             function(i, df, scores, shifting, upstream, downstream, zeroPosition, outdir,
                      location, plot.ext, acceptedLengths, legendPos, colors, addFracPlot,
                      skip.last, title) {
     heatmapListIntern <- list()
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
          if (shift == "NULL") shift <- NULL
          print(paste(i, shift, score))
          out <- file.path(outdir, paste0(name(df),"_hm_", location, "_",i , "_"))
          out <- ifelse(!is.null(shift),
                        paste0(out, shift, "_", score, plot.ext),
                        paste0(out, score, plot.ext))

          plot <- heatMap_single(region, tx, reads = get(i, envir = envExp(df), mode = "S4"), outdir = out,
                                 shifting = shift, scores = score, upstream = up, downstream = down,
                                 zeroPosition = zero, acceptedLengths = acceptedLengths,
                                 legendPos = legendPos, colors = colors, addFracPlot = addFracPlot,
                                 location = location, skip.last = skip.last,
                                 title = ifelse(title, paste(i, shift), NULL),
                                 gradient.max = gradient.max)
          heatmapListIntern <- c(heatmapListIntern, list(plot))
        }
      }
     return(heatmapListIntern)
    }, BPPARAM = BPPARAM, df = df, scores = scores, shifting = shifting, upstream = upstream,
    downstream = downstream, zeroPosition = zeroPosition, outdir = outdir,
    location = location, plot.ext = plot.ext, acceptedLengths = acceptedLengths, legendPos = legendPos,
    colors = colors, addFracPlot = addFracPlot,
    skip.last = skip.last, title = title)
    # Per experiment plot together
    if (plot.together) {
      ncols <- max(1, length(shifting))
      final <- gridExtra::grid.arrange(grobs = unlist(heatmapList, recursive = FALSE), ncol = ncols)
      ggsave(file.path(outdir, paste0(name(df), "_hm_combined_", location, plot.ext)),
             plot = final, width = scale_x * ncols,
             height = ceiling(scale_y * (length(heatmapList) / ncols)),
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
#' @param gradient.max numeric or character, default: "default", which is:
#' \code{max(coverage$score)}, the max coverage over all readlengths.
#'  If you want all plots to use same reference point
#' for max scaling, then first detect this point, look at max in plot etc,
#' and use that value, to get all plots to have same max point.
#' @return ggplot2 grob (default), data.table (if returnCoverage is TRUE)
#' @importFrom gridExtra grid.arrange
#' @family heatmaps
#' @export
heatMap_single <- function(region, tx, reads, outdir,
                           scores = "sum", upstream, downstream,  zeroPosition = upstream,
                           returnCoverage = FALSE, acceptedLengths = NULL, legendPos = "right",
                           colors = "default", addFracPlot = TRUE, location = "start site", shifting = NULL,
                           skip.last = FALSE, title = NULL, gradient.max = "default") {
  if (length(scores) != 1) stop("scores must exactly only 1 score type")
  if (!is.null(shifting)) {
    reads <- convertToOneBasedRanges(reads, method = shifting, addSizeColumn = TRUE)
  }
  if (skip.last) {
    all_lengths <- sort(unique(readWidths(reads)))
    if (!is.null(acceptedLengths))
      all_lengths <- all_lengths[all_lengths %in% acceptedLengths]
    acceptedLengths <- all_lengths[-length(all_lengths)]
  }
  drop.zero.dt <- FALSE; append.zeroes <- FALSE
  if (scores %in% c("sum", "transcriptNormalized")) {
    drop.zero.dt <- TRUE; append.zeroes <- TRUE
  }
  dt <- windowPerReadLength(region, tx, reads, upstream = upstream, downstream = downstream,
                            zeroPosition = zeroPosition, scoring = scores,
                            acceptedLengths = acceptedLengths, drop.zero.dt = drop.zero.dt,
                            append.zeroes = append.zeroes)
  if (gradient.max == "default") gradient.max <-  max(dt$score)
  plot <- coverageHeatMap(coverage = dt, scoring = scores, addFracPlot = addFracPlot,
                          xlab = paste0("Position relative to ", location), colors = colors,
                          legendPos = legendPos, title = title,
                          gradient.max = gradient.max)

  ggsave(filename = outdir, plot = plot, width = 350, height = 180,
         units = "mm", dpi = 300, limitsize = FALSE)
  # return coverage or the plot
  if (returnCoverage) return(dt)
  return(plot)
}

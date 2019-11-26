#' Make 100 bases size meta window for all libraries in experiment
#'
#' Gives you binned meta coverage plots, either saved seperatly or
#' all in one.
#' @param leaders a \code{\link{GRangesList}} of leaders (5' UTRs)
#' @param cds a \code{\link{GRangesList}} of coding sequences
#' @param trailers a \code{\link{GRangesList}} of trailers (3' UTRs)
#' @param df an ORFik experiment, to make it, see: ?experiment
#' @param outdir directory to save to
#' @param scores scoring function (default: c("sum", "zscore"))
#' @param allTogether plot all coverage plots in 1? (defualt: FALSE)
#' @param colors Which colors to use, default (skyblue4)
#' @param windowSize size of binned windows, default: 100
#' @param returnPlot return plot from function, default False
#' @return NULL, or ggplot object if returnPlot is TRUE
transcriptWindow <- function(leaders, cds, trailers, df, outdir,
                             scores = c("sum", "zscore"), allTogether = FALSE,
                             colors = rep("skyblue4", nrow(df)),
                             windowSize = min(100,
                                           min(widthPerGroup(leaders, FALSE)),
                                           min(widthPerGroup(cds, FALSE)),
                                           min(widthPerGroup(trailers, FALSE)))
                             , returnPlot = FALSE) {
  if (windowSize != 100) message(paste0("NOTE: windowSize is not 100!
                                        It is ", windowSize))

  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, leaders)
    coverage <- data.table()
    if (!allTogether) {
      stop("fix!")
      libTypes <- libraryTypes(df)
      j <- 0
      for (i in 1:nrow(df)) { # For each stage
        j <- j + 1
        print(i)
        readsList <- list()
        for (lib in libTypes) {
          # For each library of that stage (SSU, LSU, RNA-seq, RIBO-seq)
          readsList <- list(readsList, get(varNames[j]))
        }
        readsList <- readsList[-1]
        transcriptWindowPer(leaders, cds, trailers, df[i,], outdir, scores,
                            fractions, readsList)
      }
    } else { # all combined
      coverage <- data.table()
      for (i in varNames) { # For each stage
        print(i)
        coverage <- rbindlist(list(coverage,
                                   splitIn3Tx(leaders, cds, trailers,
                                              get(i), fraction = i,
                                              windowSize = windowSize)))
      }
      for(s in scores) {
        a <- windowCoveragePlot(coverage, scoring = s, colors = colors)
        ggsave(pasteDir(outdir, paste0(df@experiment,"_cp_all_", s, ".png"))
               , a, height = 10)
      }
    }
  }
  if (returnPlot) return(a)
}

#' Helper function for transcriptWindow
#'
#' Make 100 bases size meta window for one library in experiment
#'
#' Gives you binned meta coverage plots, either saved seperatly or
#' all in one.
#' @param leaders a \code{\link{GRangesList}} of leaders (5' UTRs)
#' @param cds a \code{\link{GRangesList}} of coding sequences
#' @param trailers a \code{\link{GRangesList}} of trailers (3' UTRs)
#' @param df an ORFik experiment, to make it, see: ?experiment
#' @param outdir directory to save to
#' @param scores scoring function (default: c("sum", "zscore"))
#' @param reads a GRanges / GAligment object of reads
#' @param returnCoverage return data.table with coverage (default: FALSE)
#' @param windowSize size of binned windows, default: 100
#' @return NULL, or ggplot object if returnPlot is TRUE
transcriptWindowPer <- function(leaders, cds, trailers, df,
                                outdir, scores = c("sum", "zscore"),
                                reads, returnCoverage = FALSE,
                                windowSize = 100) {

  libTypes <- libraryTypes(df)
  if (is(reads, "list") | is(reads, "GAlignmentsList") |
      is(reads, "GRangesList")) {
    if (length(libTypes) != length(reads))
      stop("not matching length of reads and lib types in df!")
  } else if(!(is(reads, "GRanges") | is(reads, "GAlignments"))) {
    stop("reads must be GRanges or GAlignments")
  }
  coverage <- data.table()

  for(i in 1:length(reads)) {
    coverage <- rbindlist(list(coverage,
                               splitIn3Tx(leaders, cds, trailers,
                                          unlist(reads[[i]]),
                                          fraction = libTypes[i],
                                          windowSize = windowSize)))
  }

  return(plotHelper(coverage, df, outdir, scores, returnCoverage))
}

#' Make 100 bases size meta window for all libraries in experiment
#'
#' Gives you binned meta coverage plots, either saved seperatly or
#' all in one.
#' @param leaders a \code{\link{GRangesList}} of leaders (5' UTRs)
#' @param cds a \code{\link{GRangesList}} of coding sequences
#' @param trailers a \code{\link{GRangesList}} of trailers (3' UTRs)
#' @param df an ORFik \code{\link{experiment}}
#' @param outdir directory to save to (default: NULL, no saving)
#' @param scores scoring function (default: c("sum", "zscore")),
#' see ?coverageScorings for possible scores.
#' @param allTogether plot all coverage plots in 1 output? (defualt: TRUE)
#' @param colors Which colors to use, default (skyblue4)
#' @param title title of ggplot
#' @param windowSize size of binned windows, default: 100
#' @param returnPlot return plot from function, default is.null(outdir),
#' so TRUE if outdir is not defined.
#' @param dfr an ORFik \code{\link{experiment}} of RNA-seq to
#' normalize against. Will add RNA normalized to plot name if this is done.
#' @param idName A character ID to add to saved name of plot,
#' if you make several plots in the same folder,
#' and same experiment, like splitting transcripts in two groups like
#' targets / nontargets etc. (default: "")
#' @param format default (".png"), do ".pdf" if you want as pdf
#' @inheritParams outputLibs
#' @export
#' @return NULL, or ggplot object if returnPlot is TRUE
#' @family experiment plots
#' @examples
#' # Make ORFik experiment
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
#' loadRegions(df) # Load leader, cds and trailers as GRangesList
#' #transcriptWindow(leaders, cds, trailers, df, outdir = "directory_to_save")
transcriptWindow <- function(leaders, cds, trailers, df, outdir = NULL,
                             scores = c("sum", "zscore"), allTogether = TRUE,
                             colors = rep("skyblue4", nrow(df)),
                             title = "Coverage metaplot",
                             windowSize = min(100,
                                           min(widthPerGroup(leaders, FALSE)),
                                           min(widthPerGroup(cds, FALSE)),
                                           min(widthPerGroup(trailers, FALSE))),
                             returnPlot = is.null(outdir),
                             dfr = NULL, idName = "", format = ".png",
                             type = "default") {
  if (windowSize != 100) message(paste0("NOTE: windowSize is not 100!
                                        It is ", windowSize))

  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, leaders, type = type)
    coverage <- data.table()
    if (!allTogether) {
      stop("fix!")
      libTypes <- libraryTypes(df)
      j <- 0
      for (i in seq(nrow(df))) { # For each stage
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
      if (!is.null(dfr)) {
        coverage <- rnaNormalize(coverage, df, dfr, cds)
        title <- paste0(title, " RNA-normalized")
      }
      for(s in scores) {
        a <- windowCoveragePlot(coverage, scoring = s, colors = colors,
                                title = title)
        if (!is.null(outdir)) {
          ggsave(pasteDir(outdir, paste0(df@experiment,"_cp_all_", s,
                                         "_", idName, format)), a,
                 height = 10)
        }
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
#' @inheritParams transcriptWindow
#' @param reads a GRanges / GAligment object of reads
#' @param returnCoverage return data.table with coverage (default: FALSE)
#' @param windowSize size of binned windows, default: 100
#' @family experiment plots
#' @return NULL, or ggplot object if returnPlot is TRUE
transcriptWindowPer <- function(leaders, cds, trailers, df,
                                outdir = NULL, scores = c("sum", "zscore"),
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

#' Meta coverage over all transcripts
#'
#' Given as single window
#' @inheritParams transcriptWindow
#' @return NULL, or ggplot object if returnPlot is TRUE
#' @family experiment plots
#'
transcriptWindow1 <- function(df, outdir = NULL,
                       scores = c("sum", "zscore"),
                       colors = rep("skyblue4", nrow(df)),
                       title = "Coverage metaplot",
                       windowSize = 100,
                       returnPlot = is.null(outdir),
                       dfr = NULL, idName = "", format = ".png",
                       type = "default") {
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, leaders, type = type)
    coverage <- data.table()
    for (f in varNames) { # For each stage
      print(f)
      temp <-  windowPerTranscript(df.rna, reads = get(f),
                                   splitIn3 = FALSE, fraction = f, )
      coverage <- rbindlist(list(coverage, temp))
    }
    if (!is.null(dfr)) {
      tx <- loadRegion(df, "tx")
      tx <- tx[widthPerGroup(tx, FALSE) >= windowSize]
      coverage <- rnaNormalize(coverage, df, dfr, )
      title <- paste0(title, " RNA-normalized")
    }
    for(s in scores) {
      a <- windowCoveragePlot(coverage, scoring = s, colors = colors,
                              title = title)
      if (!is.null(outdir)) {
        ggsave(pasteDir(outdir, paste0(df@experiment,"_cp_tx_all_", s,
                                       "_", idName, format)), a,
               height = 10)
      }
    }
  }
  if (returnPlot) return(a)
}
#' Normalize a data.table of coverage by RNA seq per position
#'
#' Normalizes per position per gene by this function:
#' (reads at position / min(librarysize, 1) * number of genes) /
#'  fpkm of that gene's RNA-seq
#'
#' Good way to compare libraries
#' @param coverage a data.table containing at least columns (count/score,
#'  position), it is possible to have additionals: (genes, fraction, feature)
#' @inheritParams transcriptWindow
#' @param tx a \code{\link{GRangesList}} of mrna transcripts
#' @param normalizeMode a character (default: "position"), how to normalize
#' library against rna library. Either on "position", normalize by
#' number of genes, sum of reads and RNA seq, on tx "region" or "feature":
#' same as position but RNA is split into the feature groups to normalize.
#' Useful if you have a list of targets and background genes.
#' @return a data.table of normalized transcripts by RNA.
#' @export
rnaNormalize <- function(coverage, df, dfr = NULL, tx, normalizeMode = "position") {
  coverage <- copy(coverage)
  if (is.null(dfr)) dfr <- df[df$libtype == "RNA",]
  if (!all(table(df$libtype) == nrow(dfr))) {
    message("RNA normalization not equal rows of experiment, returning without normalizing.")
    return(coverage)
  }

  bamVarsR <- bamVarName(df, skip.libtype = TRUE, skip.experiment = TRUE)
  bamVarsr <- bamVarName(dfr, skip.libtype = TRUE, skip.experiment = TRUE)
  matches <- bamVarsr %in% bamVarsR
  if (!any(matches)) {
    message("RNA experiment does not match any rows in df, returning without normalizing.")
    return(coverage)
  }

  outputLibs(dfr, tx)
  fullVarNameR <- bamVarName(df)
  rnaVarNames <- bamVarName(df)(dfr)[matches]
  i <- 1
  for (j in rnaVarNames) {
    fpkmsDT <- data.table(ORFik::fpkm(tx, get(j), pseudoCount = 1))

    rows <- fullVarNameR[which(bamVarsR == bamVarsr[i])]
    if (normalizeMode == "position") {
      coverage[coverage$fraction  %in% rows,
               score := (score / (min(sum(score), 1)*length(unique(genes))))
                 / fpkmsDT[genes]]
    } else if (normalizeMode == "region") {
      coverage[coverage$fraction  %in% rows,
               score := (score / fpkmsDT[genes])]
    } else if (normalizeMode == "feature") {
      print("Normalize by feature")
      f <- table(coverage$feature)
      coverage[coverage$fraction  %in% rows,
               score := (score / (min(sum(score), 1)*f[feature]))
                 / fpkmsDT[genes]]
      coverage$test <- f[coverage$feature]
    } else stop("normalizeMode must be either position, region or feature")

    i = i + 1
  }
  return(coverage)
}

#' Helper function for coverage plots
#'
#' Should only be used internally
#' @param coverage a data.table containing at least columns (count/score,
#'  position), it is possible to have additionals: (genes, fraction, feature)
#' @inheritParams transcriptWindow
#' @param returnCoverage (defualt: FALSE), return the ggplot object (TRUE)
#'  or NULL (FALSE).
#' @param title Title to give plot
#' @param plotFunction Which plot function, default: windowCoveragePlot
#' @return NULL (or ggplot object if returnCoverage is TRUE)
plotHelper <- function(coverage, df, outdir, scores, returnCoverage = FALSE,
                       title = "coverage metaplot", colors = c("skyblue4", "orange"),
                       plotFunction = "windowCoveragePlot") {
  if (!is.null(outdir)) {
    for(s in scores) {
      stage <- df$stage[1]
      type <- df$type[1]
      sample_name <- paste0(stage, "_", type, "_", s) # What happens on cds ?
      outName <- paste0(outdir, sample_name, ".png")
      if (plotFunction == "windowCoveragePlot") {
        windowCoveragePlot(coverage, output = outName, scoring = s,
                           title = title, colors = colors)
      } else if (plotFunction == "coverageHeatMap") {
        plot <- coverageHeatMap(coverage = coverage, scoring = s)
        ggsave(outName, plot = plot, width = 350, height = 180, units = "mm",
               dpi = 300, limitsize = FALSE)
      } else stop("invalid plot name")

    }
  }
  if (returnCoverage) return(coverage)
  return(NULL)
}

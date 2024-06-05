#' Make 100 bases size meta window for all libraries in experiment
#'
#' Gives you binned meta coverage plots, either saved seperatly or
#' all in one.
#' @inheritParams splitIn3Tx
#' @inheritParams outputLibs
#' @param outdir directory to save to (default: NULL, no saving)
#' @param scores scoring function (default: c("sum", "transcriptNormalized")),
#' see ?coverageScorings for possible scores.
#' @param allTogether plot all coverage plots in 1 output? (defualt: TRUE)
#' @param colors Which colors to use, default auto color from function
#' \code{\link{experiment.colors}}, new color per library type.
#' Else assign colors yourself.
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
#' @param plot.ext character, default: ".pdf". Alternatives: ".png" or ".jpg".
#' @return NULL, or ggplot object if returnPlot is TRUE
#' @export
#' @family experiment plots
#' @examples
#' df <- ORFik.template.experiment()[3,] # Only third library
#' loadRegions(df) # Load leader, cds and trailers as GRangesList
#' #transcriptWindow(leaders, cds, trailers, df, outdir = "directory/to/save")
transcriptWindow <- function(leaders, cds, trailers, df, outdir = NULL,
                             scores = c("sum", "transcriptNormalized"), allTogether = TRUE,
                             colors = experiment.colors(df),
                             title = "Coverage metaplot",
                             windowSize = min(100,
                                           min(widthPerGroup(leaders, FALSE)),
                                           min(widthPerGroup(cds, FALSE)),
                                           min(widthPerGroup(trailers, FALSE))),
                             returnPlot = is.null(outdir),
                             dfr = NULL, idName = "", plot.ext = ".pdf",
                             type = "ofst", is.sorted = FALSE, drop.zero.dt = TRUE,
                             verbose = TRUE, force = TRUE,
                             BPPARAM = bpparam()) {
  if (windowSize != 100)
    message(paste0("NOTE: windowSize is not 100! It is: ", windowSize))

  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, chrStyle = leaders, type = type, verbose = verbose, force = force, BPPARAM = BPPARAM)
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
                            fractions, readsList, drop.zero.dt = drop.zero.dt)
      }
    } else { # all combined
      coverage <- bplapply(varNames, function(x, leaders, cds, trailers,
                                              windowSize, is.sorted, drop.zero.dt) {
        message(x)
        splitIn3Tx(leaders, cds, trailers,
                   get(x, mode = "S4", envir = envExp(df)), fraction = x,
                   windowSize = windowSize, is.sorted = is.sorted,
                   drop.zero.dt = drop.zero.dt)
      }, leaders = leaders, cds = cds, trailers = trailers,
         is.sorted = is.sorted, windowSize = windowSize, drop.zero.dt = drop.zero.dt,
         BPPARAM = BPPARAM)

      coverage <- rbindlist(coverage)
      if (!is.null(dfr)) {
        coverage <- rnaNormalize(coverage, df, dfr, cds)
        title <- paste0(title, " RNA-normalized")
      }
      coverage[, feature := factor(feature, levels = c("leaders", "cds", "trailers"),
                                    ordered = TRUE)]
      a <- bplapply(scores, function(s, coverage, colors, title,
                                            idName, outdir, plot.ext, df) {
        message(s)
        a <- windowCoveragePlot(coverage, scoring = s, colors = colors,
                                title = title)
        if (!is.null(outdir)) {
          idName <- ifelse(idName == "", "", paste0("_", idName))
          ggsave(pasteDir(outdir, paste0(name(df),"_cp_all_", s,
                                         idName, plot.ext)),
                 a,
                 width = 6,
                 height = 10 + floor(nrow(df) / 5))
        }
        return(a)
      }, coverage = coverage, colors = colors, title = title,
      idName = idName, outdir = outdir, plot.ext = plot.ext, df = df,
      BPPARAM = BPPARAM)
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
#' @param reads a GRanges / GAligment object of reads, can also be a list of those.
#' @param returnCoverage return data.table with coverage (default: FALSE)
#' @param windowSize size of binned windows, default: 100
#' @family experiment plots
#' @return NULL, or ggplot object if returnPlot is TRUE
#' @keywords internal
transcriptWindowPer <- function(leaders, cds, trailers, df,
                                outdir = NULL, scores = c("sum", "zscore"),
                                reads, returnCoverage = FALSE,
                                windowSize = 100, drop.zero.dt = TRUE,
                                BPPARAM = bpparam()) {
  libTypes <- libraryTypes(df)
  if (is(reads, "list") | is(reads, "GAlignmentsList") |
      is(reads, "GRangesList")) {
    if (length(libTypes) != length(reads))
      stop("not matching length of reads and lib types in df!")
  } else if(!(is(reads, "GRanges") | is(reads, "GAlignments"))) {
    stop("reads must be GRanges or GAlignments")
  }

  coverage <- bplapply(reads, function(x, leaders, cds, trailers,
                                          fraction, windowSize, drop.zero.dt) {
    message(names(x)[1])
    splitIn3Tx(leaders, cds, trailers, x, fraction = names(x)[1],
               windowSize = windowSize, drop.zero.dt = drop.zero.dt)
  }, leaders = leaders, cds = cds, trailers = trailers,
     windowSize = windowSize, drop.zero.dt = drop.zero.dt, BPPARAM = BPPARAM)
  coverage <- rbindlist(coverage)

  return(plotHelper(coverage, df, outdir, scores, returnCoverage))
}

#' Meta coverage over all transcripts
#'
#' Given as single window
#' @inheritParams transcriptWindow
#' @return NULL, or ggplot object if returnPlot is TRUE
#' @family experiment plots
#' @keywords internal
transcriptWindow1 <- function(df, outdir = NULL,
                       scores = c("sum", "zscore"),
                       colors = experiment.colors(df),
                       title = "Coverage metaplot",
                       windowSize = 100,
                       returnPlot = is.null(outdir),
                       dfr = NULL, idName = "", plot.ext = ".pdf",
                       type = "ofst", drop.zero.dt = drop.zero.dt,
                       BPPARAM = bpparam()) {
  dfl <- df
  if(!is(dfl, "list")) dfl <- list(dfl)
  for (df in dfl) {
    varNames <- bamVarName(df)
    outputLibs(df, type = type)

    coverage <- bplapply(varNames, function(x, df, windowSize, drop.zero.dt) {
      message(x)
      windowPerTranscript(df, reads = get(x, mode = "S4", envir = envExp(df)),
                          splitIn3 = FALSE, fraction = x,
                          windowSize = windowSize, drop.zero.dt = drop.zero.dt)
    }, df = df, windowSize = windowSize, drop.zero.dt = drop.zero.dt,
    BPPARAM = BPPARAM)
    coverage <- rbindlist(coverage)

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
        idName <- ifelse(idName == "", "", paste0("_", idName))
        ggsave(pasteDir(outdir, paste0(df@experiment,"_cp_tx_all_", s,
                                       idName, plot.ext)), a,
               width = 6,
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
#' @keywords internal
plotHelper <- function(coverage, df, outdir, scores, returnCoverage = FALSE,
                       title = "coverage metaplot", plot.ext = ".pdf", colors = c("skyblue4", "orange"),
                       plotFunction = "windowCoveragePlot") {
  if (!is.null(outdir)) {
    for(s in scores) {
      stage <- df$stage[1]
      type <- df$type[1]
      sample_name <- paste0(stage, "_", type, "_", s) # What happens on cds ?
      outName <- paste0(outdir, sample_name, plot.ext)
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

#' Decide color for libraries by grouping
#'
#' Pick the grouping wanted for colors, by default only group by libtype.
#' Like RNA-seq(skyblue4) and Ribo-seq(orange).
#' @inheritParams bamVarName
#' @param color_list a character vector of colors, default "default".
#' That is the vector c("skyblue4", 'orange', "green", "red", "gray",
#'  "yellow", "blue", "red2", "orange3"). Picks number of colors needed to
#'  make groupings have unique color
#' @return a character vector of colors
#' @export
experiment.colors <- function(df, color_list = "default",
                              skip.libtype = FALSE, skip.stage = TRUE,
                              skip.replicate = TRUE,
                              skip.fraction = TRUE, skip.condition = TRUE) {
  lib_types <- bamVarName(df, skip.replicate, skip.condition,
                          skip.stage, skip.fraction, skip.libtype)
  unique_lib_types <- sort(unique(lib_types))
  n_colors <- length(unique_lib_types)

  color_list <- c("skyblue4", 'orange', "green", "red", "gray",
                  "yellow", "blue", "red2", "orange3")
  while (n_colors > length(color_list)) {
    color_list <- rep(color_list, 2)
  }

  return(color_list[match(lib_types, unique_lib_types)])
}

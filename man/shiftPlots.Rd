% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shift_footprints.R
\name{shiftPlots}
\alias{shiftPlots}
\title{Plot shifted heatmaps per library}
\usage{
shiftPlots(
  df,
  output = NULL,
  title = "Ribo-seq",
  scoring = "transcriptNormalized",
  BPPARAM = bpparam()
)
}
\arguments{
\item{df}{an ORFik \code{\link{experiment}}}

\item{output}{name to save file, full path. (Default NULL) No saving.}

\item{title}{Title for top of plot, default "Ribo-seq".
A more informative name could be "Ribo-seq zebrafish Chew et al. 2013"}

\item{scoring}{which scoring scheme to use for heatmap, default
"transcriptNormalized". Some alternatives: "sum", "zscore".}

\item{BPPARAM}{how many cores/threads to use? default: bpparam()}
}
\value{
a ggplot2 grob object
}
\description{
A good validation for you p-shifting
}
\examples{
#df <- read.experiment("ORFik_example")
#shiftFootprintsByExperiment(df, output_format = "bedo)
#shiftPlots(df, title = "Ribo-seq Human ORFik et al. 2020")
}
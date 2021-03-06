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
  pShifted = TRUE,
  upstream = if (pShifted) 5 else 20,
  downstream = if (pShifted) 20 else 5,
  type = "bar",
  addFracPlot = TRUE,
  BPPARAM = bpparam()
)
}
\arguments{
\item{df}{an ORFik \code{\link{experiment}}}

\item{output}{name to save file, full path. (Default NULL) No saving.
Sett to "auto" to save to QC_STATS folder of experiment named:
"pshifts_barplots.png" or "pshifts_heatmaps.png" depending on type argument.
Folder must exist!}

\item{title}{Title for top of plot, default "Ribo-seq".
A more informative name could be "Ribo-seq zebrafish Chew et al. 2013"}

\item{scoring}{which scoring scheme to use for heatmap, default
"transcriptNormalized". Some alternatives: "sum", "zscore".}

\item{pShifted}{a logical (TRUE), are Ribo-seq reads p-shifted to size
1 width reads? If upstream and downstream is set, this argument is
irrelevant. So set to FALSE if this is not p-shifted Ribo-seq.}

\item{upstream}{an integer (5), relative region to get upstream from.}

\item{downstream}{an integer (20), relative region to get downstream from}

\item{type}{character, default "bar". Plot as faceted bars,
gives more detailed information of read lengths,
but harder to see patterns over multiple read lengths.
Alternative: "heatmaps", better overview of patterns over
multiple read lengths.}

\item{addFracPlot}{logical, default TRUE, add positional sum plot on top
per heatmap.}

\item{BPPARAM}{how many cores/threads to use? default: bpparam()}
}
\value{
a ggplot2 grob object
}
\description{
A good validation for you p-shifting, to see shifts are corresponding
and close to the CDS TIS.
}
\examples{
df <- ORFik.template.experiment()
df <- df[3,] #lets only p-shift RFP sample at index 3
#shiftFootprintsByExperiment(df, output_format = "bedo)
#shiftPlots(df, title = "Ribo-seq Human ORFik et al. 2020")
}

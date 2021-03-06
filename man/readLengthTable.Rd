% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/report_helpers.R
\name{readLengthTable}
\alias{readLengthTable}
\title{Make table of readlengths}
\usage{
readLengthTable(df, output.dir = NULL, type = "ofst", BPPARAM = bpparam())
}
\arguments{
\item{df}{an ORFik \code{\link{experiment}}}

\item{output.dir}{NULL or character path, default: NULL, plot not saved to disc.
If defined saves plot to that directory with the name "./readLengths.csv".}

\item{type}{character, default: "ofst". Type of library:
either "default", usually bam format (the one you gave to experiment),
"pshifted" pshifted reads, "ofst", "bed", "bedo" optimized bed, or "wig"}

\item{BPPARAM}{how many cores/threads to use? default: bpparam().
To see number of threads used, do \code{bpparam()$workers}}
}
\value{
a data.table object of the the read length data with columns:
\code{c("sample", "sample_id", "read length", "counts",
 "counts_per_sample", "perc_of_counts_per_sample")}
}
\description{
Summarizing all libraries in experiment,
make a table of proportion of read lengths.
}

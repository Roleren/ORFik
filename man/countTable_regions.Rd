% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SummarizedExperiment_helpers.R
\name{countTable_regions}
\alias{countTable_regions}
\title{Make a list of count matrices from experiment}
\usage{
countTable_regions(
  df,
  out.dir = dirname(df$filepath[1]),
  longestPerGene = TRUE,
  geneOrTxNames = "tx",
  regions = c("mrna", "leaders", "cds", "trailers"),
  type = "count",
  lib.type = "ofst",
  weight = "score",
  BPPARAM = bpparam()
)
}
\arguments{
\item{df}{an ORFik \code{\link{experiment}}}

\item{out.dir}{optional output directory, default:
\code{dirname(df$filepath[1])}.
Will make a folder called "QC_STATS" with all results in this directory.}

\item{longestPerGene}{a logical (default TRUE), if FALSE all transcript
isoforms per gene.}

\item{geneOrTxNames}{a character vector (default "tx"), should row names
keep trancript names ("tx") or change to gene names ("gene")}

\item{regions}{a character vector, default:
 c("mrna", "leaders", "cds", "trailers"), make raw count matrices
of whole regions specified.}

\item{type}{default: "count" (raw counts matrix), alternative is "fpkm",
"log2fpkm" or "log10fpkm"}

\item{lib.type}{a character(default: "default"), load files in experiment
or some precomputed variant, either "ofst", "bedo", "bedoc" or "pshifted".
These are made with ORFik:::convertLibs() or shiftFootprintsByExperiment().
Can also be custom user made folders inside the experiments bam folder.}

\item{weight}{numeric or character, a column to score overlaps by. Default "score",
will check for a metacolumn called "score" in libraries. If not found,
will not use weights.}

\item{BPPARAM}{how many cores/threads to use? default: bpparam()}
}
\value{
a list of data.table, 1 data.table per region. The regions
will be the names the list elements.
}
\description{
Make a list of count matrices from experiment
}
\seealso{
Other countTable: 
\code{\link{countTable}()}
}
\concept{countTable}

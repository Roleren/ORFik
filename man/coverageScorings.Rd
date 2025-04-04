% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage_helpers.R
\name{coverageScorings}
\alias{coverageScorings}
\title{Add a coverage scoring scheme}
\usage{
coverageScorings(coverage, scoring = "zscore", copy.dt = TRUE)
}
\arguments{
\item{coverage}{a data.table containing at least columns (count, position),
it is possible to have additionals: (genes, fraction, feature)}

\item{scoring}{a character, one of (zscore, transcriptNormalized,
mean, median, sum, log2sum, log10sum, sumLength, meanPos and frameSum,
periodic, NULL). More info in details}

\item{copy.dt}{logical TRUE, copy object, to avoid overwriting original object.
Set to false to run function using reference to object,
a speed up if original object is not needed.}
}
\value{
a data.table with new scores (size dependent on score used)
}
\description{
Different scorings and groupings of a coverage representation.
}
\details{
Usually output of metaWindow or scaledWindowPositions is input in this
function.

Content of coverage data.table:
It must contain the count and position columns.

genes column: If you have multiple windows, the genes column must define
which gene/transcript grouping the different counts belong to. If there is
only a meta window or only 1 gene/transcript, then this column is
not needed.

fraction column: If you have coverage of i.e RNA-seq and Ribo-seq, or TCP
-seq of large and small subunite, divide into fractions.
Like factor(RNA, RFP)

feature column: If gene group is subdivided into parts, like gene is
transcripts, and feature column can be c(leader, cds, trailer) etc.

Given a data.table coverage of counts, add a scoring scheme.
per: the grouping given, if genes is defined,
group by per gene in default scoring.\cr
Scorings:\cr
\itemize{
\item zscore (count-mean(window))/sd(window) per) # Outlier detection
\item modzscore (count-median(window))/mad(window) per) # Outlier detection for signal with extreme values
\item transcriptNormalized (sum(count / sum of counts per)) this is scaled per group
\item mean (mean(count per))
\item median (median(count per))
\item sum (count per)
\item log2sum (count per)
\item log10sum (count per)
\item sumLength (count per) / number of windows
\item meanPos (mean per position per gene) used in scaledWindowPositions
\item sumPos (sum per position per gene) used in scaledWindowPositions
\item frameSum (sum per frame per gene) used in ORFScore
\item frameSumPerL (sum per frame per read length)
\item frameSumPerLG (sum per frame per read length per gene)
\item fracPos (fraction of counts per position per gene) non scaled version of transcriptNormalized
\item periodic (Fourier transform periodicity of meta coverage per fraction)
\item NULL (no grouping, return input directly)
}
}
\examples{
dt <- data.table::data.table(count = c(4, 1, 1, 4, 2, 3),
                             position = c(1, 2, 3, 4, 5, 6))
coverageScorings(dt, scoring = "zscore")

# with grouping gene
dt$genes <- c(rep("tx1", 3), rep("tx2", 3))
coverageScorings(dt, scoring = "zscore")

}
\seealso{
Other coverage: 
\code{\link{metaWindow}()},
\code{\link{regionPerReadLength}()},
\code{\link{scaledWindowPositions}()},
\code{\link{windowPerReadLength}()}
}
\concept{coverage}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grl_region.R
\name{startCodons}
\alias{startCodons}
\title{Get the Start codons(3 bases) from a GRangesList of orfs grouped by orfs}
\usage{
startCodons(grl, is.sorted = FALSE)
}
\arguments{
\item{grl}{a \code{\link{GRangesList}} object}

\item{is.sorted}{a boolean, a speedup if you know the ranges are sorted}
}
\value{
a GRangesList of start codons, since they might be split on exons
}
\description{
In ATGTTTTGA, get the positions ATG.
It takes care of exons boundaries, with exons < 3 length.
}
\examples{
gr_plus <- GRanges(seqnames = "chr1",
                   ranges = IRanges(c(7, 14), width = 3),
                   strand = "+")
gr_minus <- GRanges(seqnames = "chr2",
                    ranges = IRanges(c(4, 1), c(9, 3)),
                    strand = "-")
grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
startCodons(grl, is.sorted = FALSE)

}
\seealso{
Other ORFHelpers: 
\code{\link{defineTrailer}()},
\code{\link{longestORFs}()},
\code{\link{mapToGRanges}()},
\code{\link{orfID}()},
\code{\link{startSites}()},
\code{\link{stopCodons}()},
\code{\link{stopSites}()},
\code{\link{txNames}()},
\code{\link{uniqueGroups}()},
\code{\link{uniqueOrder}()}
}
\concept{ORFHelpers}

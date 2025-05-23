% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grl_region.R
\name{startSites}
\alias{startSites}
\title{Get the start sites from a GRangesList of orfs grouped by orfs}
\usage{
startSites(grl, asGR = FALSE, keep.names = FALSE, is.sorted = FALSE)
}
\arguments{
\item{grl}{a \code{\link{GRangesList}} object}

\item{asGR}{a boolean, return as GRanges object}

\item{keep.names}{a logical (FALSE), keep names of input.}

\item{is.sorted}{a speedup, if you know the ranges are sorted}
}
\value{
if asGR is False, a vector, if True a GRanges object
}
\description{
In ATGTTTTGG, get the position of the A.
}
\examples{
gr_plus <- GRanges(seqnames = c("chr1", "chr1"),
                   ranges = IRanges(c(7, 14), width = 3),
                   strand = c("+", "+"))
gr_minus <- GRanges(seqnames = c("chr2", "chr2"),
                    ranges = IRanges(c(4, 1), c(9, 3)),
                    strand = c("-", "-"))
grl <- GRangesList(tx1 = gr_plus, tx2 = gr_minus)
startSites(grl, is.sorted = FALSE)

}
\seealso{
Other ORFHelpers: 
\code{\link{defineTrailer}()},
\code{\link{longestORFs}()},
\code{\link{mapToGRanges}()},
\code{\link{orfID}()},
\code{\link{startCodons}()},
\code{\link{stopCodons}()},
\code{\link{stopSites}()},
\code{\link{txNames}()},
\code{\link{uniqueGroups}()},
\code{\link{uniqueOrder}()}
}
\concept{ORFHelpers}

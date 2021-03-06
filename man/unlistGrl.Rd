% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grl_helpers.R
\name{unlistGrl}
\alias{unlistGrl}
\title{Safe unlist}
\usage{
unlistGrl(grl)
}
\arguments{
\item{grl}{a GRangesList}
}
\value{
a GRanges object
}
\description{
Same as [AnnotationDbi::unlist2()], keeps names correctly.
Two differences is that if grl have no names, it will not
make integer names, but keep them as null. Also if the GRangesList has names
, and also the GRanges groups, then the GRanges group names will be kept.
}
\examples{
ORF <- GRanges(seqnames = "1",
               ranges = IRanges(start = c(1, 10, 20),
                                end = c(5, 15, 25)),
               strand = "+")
grl <- GRangesList(tx1_1 = ORF)
unlistGrl(grl)

}

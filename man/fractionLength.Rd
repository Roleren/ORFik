% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sequence_features.R
\name{fractionLength}
\alias{fractionLength}
\title{Fraction Length}
\usage{
fractionLength(grl, tx_len = widthPerGroup(tx, TRUE), tx = NULL)
}
\arguments{
\item{grl}{a \code{\link{GRangesList}} object
with usually either leaders,
cds', 3' utrs or ORFs. ORFs are a special case, see argument tx_len}

\item{tx_len}{the transcript lengths of the transcripts,
a named (tx names) vector of integers.
If you have the transcripts as GRangesList,
call `ORFik:::widthPerGroup(tx, TRUE)`.

If you used CageSeq to reannotate leaders, then the tss for the the leaders
have changed, therefore the tx lengths have changed. To account for that
call: `tx_len <- widthPerGroup(extendLeaders(tx, cageFiveUTRs))`
and calculate fraction length using `fractionLength(grl, tx_len)`.}

\item{tx}{default NULL, a \code{\link{GRangesList}} object of transcript to get
lengths from. Pass in for wrapping to widths inside the function.}
}
\value{
a numeric vector of ratios
}
\description{
Fraction Length is defined as
\preformatted{(widths of grl)/tx_len}
so that each group in
the grl is divided by the corresponding transcript.
}
\examples{
ORF <- GRanges(seqnames = "1",
               ranges = IRanges(start = c(1, 10, 20), end = c(5, 15, 25)),
               strand = "+")
grl <- GRangesList(tx1_1 = ORF)
# grl must have same names as cds + _1 etc, so that they can be matched.
tx <-  GRangesList(tx1 = GRanges("1", IRanges(1, 50), "+"))
fractionLength(grl, tx = tx)

}
\references{
doi: 10.1242/dev.098343
}
\seealso{
Other features: 
\code{\link{computeFeatures}()},
\code{\link{computeFeaturesCage}()},
\code{\link{countOverlapsW}()},
\code{\link{disengagementScore}()},
\code{\link{distToCds}()},
\code{\link{distToTSS}()},
\code{\link{entropy}()},
\code{\link{floss}()},
\code{\link{fpkm}()},
\code{\link{fpkm_calc}()},
\code{\link{initiationScore}()},
\code{\link{insideOutsideORF}()},
\code{\link{isInFrame}()},
\code{\link{isOverlapping}()},
\code{\link{kozakSequenceScore}()},
\code{\link{orfScore}()},
\code{\link{rankOrder}()},
\code{\link{ribosomeReleaseScore}()},
\code{\link{ribosomeStallingScore}()},
\code{\link{startRegion}()},
\code{\link{startRegionCoverage}()},
\code{\link{stopRegion}()},
\code{\link{subsetCoverage}()},
\code{\link{translationalEff}()}
}
\concept{features}

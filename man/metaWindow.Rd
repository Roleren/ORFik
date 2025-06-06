% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/coverage_helpers.R
\name{metaWindow}
\alias{metaWindow}
\title{Calculate meta-coverage of reads around input GRanges/List object.}
\usage{
metaWindow(
  x,
  windows,
  scoring = "sum",
  withFrames = FALSE,
  zeroPosition = NULL,
  scaleTo = 100,
  fraction = NULL,
  feature = NULL,
  forceUniqueEven = !is.null(scoring),
  forceRescale = TRUE,
  weight = "score",
  drop.zero.dt = FALSE,
  append.zeroes = FALSE
)
}
\arguments{
\item{x}{GRanges/GAlignment object of your reads.
Remember to resize them beforehand to width of 1 to focus on
5' ends of footprints etc, if that is wanted.}

\item{windows}{GRangesList or GRanges of your ranges}

\item{scoring}{a character, default: "sum", one of
(zscore, transcriptNormalized, mean, median, sum, sumLength, NULL),
see ?coverageScorings for info and more alternatives.}

\item{withFrames}{a logical (TRUE), return positions with the 3 frames,
relative to zeroPosition. zeroPosition is frame 0.}

\item{zeroPosition}{an integer DEFAULT (NULL), the point if all windows
are equal size, that should be set to position 0. Like leaders and
cds combination, then 0 is the TIS and -1 is last base in leader. NOTE!:
if not all windows have equal width, this will be ignored. If all have
equal width and zeroPosition is NULL, it is set to as.integer(width / 2).}

\item{scaleTo}{an integer (100), if windows have different size,
a meta window can not directly be created, since a meta window must
have equal size for all windows. Rescale (bin) all windows to scaleTo.
i.e c(1,2,3) -> size 2 -> coverage of position c(1, mean(2,3)) etc.}

\item{fraction}{a character/integer (NULL), the fraction i.e
(27) for read length 27, or ("LSU") for large sub-unit TCP-seq.}

\item{feature}{a character string, info on region. Usually either
gene name, transcript part like cds, leader, or CpG motifs etc.}

\item{forceUniqueEven, }{a logical (TRUE), if TRUE; require that all windows
are of same width and even. To avoid bugs. FALSE if score is NULL.}

\item{forceRescale}{logical, default TRUE. If TRUE, if
\code{unique(widthPerGroup(windows))} has length > 1, it will force all
windows to width of the \code{scaleTo} argument, making a binned meta
coverage.}

\item{weight}{(default: 'score'), if defined a character name
of valid meta column in subject. GRanges("chr1", 1, "+", score = 5),
would mean score column tells that this alignment region was found 5 times.
Formats which loads a score column like this:
Bigwig, wig, ORFik ofst, collapsed bam, bedoc and .bedo.
As do CAGEr CAGE files and many other package formats.
You can also assign a score column manually.}

\item{drop.zero.dt}{logical FALSE, if TRUE and as.data.table is TRUE,
remove all 0 count positions.
This greatly speeds up and most importantly, greatly reduces memory usage.
Will not change any plots, unless 0 positions are used in some sense.
(mean, median, zscore coverage will only scale differently)}

\item{append.zeroes}{logical, default FALSE. If TRUE and drop.zero.dt
is TRUE and all windows have equal length, it will add back 0 values after transformation.
Sometimes needed for correct plots, if TRUE, will call abort if not all
windows are equal length!}
}
\value{
A data.table with scored counts (score) of
reads mapped to positions (position) specified in windows along with
frame (frame) per gene (genes) per library (fraction) per transcript region
(feature). Column that does not apply is not given, but position and (score/count)
is always returned.
}
\description{
Sums up coverage over set of GRanges objects as a meta representation.
}
\examples{
library(GenomicRanges)
windows <- GRangesList(GRanges("chr1", IRanges(c(50, 100), c(80, 200)),
                               "-"))
x <- GenomicRanges::GRanges(
  seqnames = "chr1",
  ranges =  IRanges::IRanges(c(100, 180), c(200, 300)),
  strand = "-")
metaWindow(x, windows, withFrames = FALSE)

}
\seealso{
Other coverage: 
\code{\link{coverageScorings}()},
\code{\link{regionPerReadLength}()},
\code{\link{scaledWindowPositions}()},
\code{\link{windowPerReadLength}()}
}
\concept{coverage}

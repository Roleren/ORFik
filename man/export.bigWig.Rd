% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_exports.R
\name{export.bigWig}
\alias{export.bigWig}
\title{Export as bigWig format}
\usage{
export.bigWig(
  x,
  file,
  split.by.strand = TRUE,
  is_pre_collapsed = FALSE,
  seq_info = seqinfo(x)
)
}
\arguments{
\item{x}{A GRangesList, GAlignment GAlignmentPairs with score column.
Will be converted to 5' end position of original range. If score column
does not exist, will group ranges and give replicates as score column.
Since bigWig needs a score column to represent counts!}

\item{file}{a character path to valid output file name}

\item{split.by.strand}{logical, default TRUE. Split bigWig into 2 files,
one for each strand.}

\item{is_pre_collapsed}{logical, default FALSE. Have you already
collapsed reads with collapse.by.scores,
so each positions is only in 1 GRanges object with
a score column per readlength?
Set to TRUE, only if you are sure, will give a speedup.}

\item{seq_info}{a Seqinfo object, default seqinfo(x).
Must have non NA seqlengths defined!}
}
\value{
invisible(NULL) (File is saved as 2 .bigWig files)
}
\description{
Will create 2 files, 1 for + strand (*_forward.bigWig)
and 1 for - strand (*_reverse.bigWig). If all
ranges are * stranded, will output 1 file.
Can be direct input for ucsc browser or IGV
}
\examples{
x <- c(GRanges("1", c(1,3,5), "-"), GRanges("1", c(1,3,5), "+"))
seqlengths(x) <- 10
file <- file.path(tempdir(), "rna.bigWig")
# export.bigWig(x, file)
# export.bigWig(covRleFromGR(x), file)
}
\references{
https://genome.ucsc.edu/goldenPath/help/bigWig.html
}
\seealso{
Other utils: 
\code{\link{bedToGR}()},
\code{\link{convertToOneBasedRanges}()},
\code{\link{export.bed12}()},
\code{\link{export.fstwig}()},
\code{\link{export.wiggle}()},
\code{\link{fimport}()},
\code{\link{findFa}()},
\code{\link{fread.bed}()},
\code{\link{optimizeReads}()},
\code{\link{readBam}()},
\code{\link{readBigWig}()},
\code{\link{readWig}()}
}
\concept{utils}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/grl_helpers.R
\name{coverageByTranscriptFST}
\alias{coverageByTranscriptFST}
\title{Get coverage from fst large coverage format}
\usage{
coverageByTranscriptFST(grl, fst_index, columns = NULL)
}
\arguments{
\item{grl}{a GRangesList}

\item{fst_index}{a path to an existing fst index file}

\item{columns}{NULL or character, default NULL. Else must be a subset of
names in the fst files. Run ids etc.}
}
\value{
a list, each element is a data.table of coverage
}
\description{
Get coverage from fst large coverage format
}
\examples{
library(data.table)
library(ORFik)
grl <- GRangesList("1:1-5:+")
tempdir <- tempdir()
fst_index <- file.path(tempdir, "coverage_index.fst")
mock_run_names <- c("SRR1010101", "SRR1010102", "SRR1010103")
coverage_file <- file.path(tempdir, paste0("coverage_1_part1_",
 c("forward", "reverse"), ".fst"))
mock_coverage <- setnames(setDT(lapply(mock_run_names, function(x) {
 sample(seq(0, 100), 100, replace = TRUE, prob = c(0.95, rep(0.01, 100)))})),
 mock_run_names)
mock_index <- data.table(chr = "1", start = 1, end = nrow(mock_coverage),
 file_forward = coverage_file[1], file_reverse = coverage_file[2])

fst::write_fst(mock_index, fst_index)
fst::write_fst(mock_coverage, coverage_file[1])

coverageByTranscriptFST(grl, fst_index)
coverageByTranscriptFST(grl, fst_index, c("SRR1010101", "SRR1010102"))
}

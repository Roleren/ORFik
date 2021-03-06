% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{getGAlignments}
\alias{getGAlignments}
\title{Internal GAlignments loader from fst data.frame}
\usage{
getGAlignments(df)
}
\arguments{
\item{df}{a data.frame with columns minimum 4 columns:
seqnames, start ("pos" in final GA object), strand and width.\cr
Additional columns will be assigned as meta columns}
}
\value{
GAlignments object
}
\description{
Internal GAlignments loader from fst data.frame
}

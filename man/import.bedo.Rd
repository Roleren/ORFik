% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_imports.R
\name{import.bedo}
\alias{import.bedo}
\title{Load GRanges object from .bedo}
\usage{
import.bedo(path)
}
\arguments{
\item{path}{a character, location on disc (full path)}
}
\value{
GRanges object
}
\description{
.bedo is .bed ORFik, an optimized bed format for coverage reads with read lengths
.bedo is a text based format with columns (6 maximum):\cr
1. chromosome\cr 2. start\cr 3. end\cr 4. strand\cr
5. ref width (cigar # M's, match/mismatch total)\cr
6. duplicates of that read\cr
}
\details{
Positions are 1-based, not 0-based as .bed.
export with export.bedo
}

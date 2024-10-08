% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_collapse.R
\name{ofst_merge}
\alias{ofst_merge}
\title{Merge multiple ofst file}
\usage{
ofst_merge(
  file_paths,
  lib_names = sub(pattern = "\\\\.ofst$", replacement = "", basename(file_paths)),
  keep_all_scores = TRUE,
  keepCigar = TRUE,
  sort = TRUE,
  max_splits = 20
)
}
\arguments{
\item{file_paths}{Full path to .ofst files wanted to merge}

\item{lib_names}{character, the name to give the resulting score columns.
Default: \code{sub(pattern = "\\.ofst$", replacement = "", basename(file_paths))}}

\item{keep_all_scores}{logical, default TRUE, keep all library scores in the merged file. These
score columns are named the libraries full name from \code{bamVarName(df)}.}

\item{keepCigar}{logical, default TRUE. If CIGAR is defined, keep column. Setting
to FALSE compresses the file much more usually.}

\item{sort}{logical, default TRUE. Sort the ranges. Will make the file smaller and
faster to load, but some additional merging time is added.}

\item{max_splits}{integer, default 20. If number of rows to merge > 2^31,
how many times can you allow split merging to try to "rescue" the merging
process?}
}
\value{
a data.table of merged result, it is merged on all columns except "score".
The returned file will contain the scores of each file + the aggregate sum score.
}
\description{
Collapses and sums the score column of each ofst file
It is required that each file is of same ofst type.
That is if one file has cigar information, all must have it.
}

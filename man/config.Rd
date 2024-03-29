% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_directory_config.R
\name{config}
\alias{config}
\title{Read directory config for ORFik experiments}
\usage{
config(
  file = config_file(old_config_location = old_config_location),
  old_config_location = "~/Bio_data/ORFik_config.csv"
)
}
\arguments{
\item{file}{location of config csv, default:
config_file(old_config_location = old_config_location)}

\item{old_config_location}{path, old config location before BiocFileCache
implementation.
Will copy this to cache directory and delete old version.
This is done to follow bioc rules on not writing to user home directory.}
}
\value{
a named character vector of length 3
}
\description{
Defines a folder for:
1. fastq files (raw data)\cr
2. bam files (processed data)\cr
3. references (organism annotation and STAR index)\cr
4. experiments (Location to store and load all \code{\link{experiment}} .csv files)
Update or use another config using \code{config.save()} function.
}
\examples{
## Make with default config path
#config()
}

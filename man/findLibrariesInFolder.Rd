% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experiment_files.R
\name{findLibrariesInFolder}
\alias{findLibrariesInFolder}
\title{Get all library files in folder/folders of given types}
\usage{
findLibrariesInFolder(dir, types, pairedEndBam = FALSE)
}
\arguments{
\item{dir}{Which directory / directories to create experiment from,
must be a directory with NGS data from your experiment. Will include
all files of file type specified by "types" argument. So do not mix
files from other experiments in the same folder!}

\item{types}{Default \code{c("bam", "bed", "wig", "bigWig","ofst")},
which types of libraries to allow as NGS data.}

\item{pairedEndBam}{logical FALSE, else TRUE, or a logical list of
TRUE/FALSE per library you see will be included (run first without and check
what order the files will come in) 1 paired end file, then two single will
be c(T, F, F). If you have a SRA metadata csv file, you can set this argument to
study$LibraryLayout == "PAIRED", where study is the SRA metadata for all files
that was aligned.}
}
\value{
(data.table) All files found from types parameter.
With 2 extra column (logical), is it wig pairs, and paired bam files.
}
\description{
Will try to guess paired / unpaired wig, bed, bam files.
}
\details{
Set pairedEndBam if you have paired end reads as a single bam file.
}
\keyword{internal}

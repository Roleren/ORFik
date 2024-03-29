% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RiboORFs_detect_and_classify.R
\name{riboORFsFolder}
\alias{riboORFsFolder}
\title{Define folder for prediction output}
\usage{
riboORFsFolder(df, parrent_dir = resFolder(df))
}
\arguments{
\item{df}{ORFik experiment}

\item{parrent_dir}{Parrent directory of computed study results, default:
resFolder(df)}
}
\value{
a file path (full path)
}
\description{
Define folder for prediction output
}
\examples{
df <- ORFik.template.experiment()
df <- df[df$libtype == "RFP",][c(1,2),]
riboORFsFolder(df)
riboORFsFolder(df, tempdir())
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils_format_conversion.R
\name{convert_to_covRleList}
\alias{convert_to_covRleList}
\title{Convert libraries to covRleList objects}
\usage{
convert_to_covRleList(
  df,
  in_files = filepath(df, "pshifted"),
  out_dir = file.path(libFolder(df), "cov_RLE_List"),
  out_dir_merged = file.path(libFolder(df), "cov_RLE"),
  split.by.strand = TRUE,
  seq_info = seqinfo(df),
  weight = "score",
  format = "qs",
  verbose = TRUE
)
}
\arguments{
\item{df}{an ORFik \code{\link{experiment}}, or NULL is allowed if
both in_files and out_dir is specified manually.}

\item{in_files}{paths to input files, default pshifted files:
\code{filepath(df, "pshifted")} in ofst format}

\item{out_dir}{paths to output files, default
\code{file.path(libFolder(df), "cov_RLE_List")}.}

\item{out_dir_merged}{character vector of paths, default:
\code{file.path(libFolder(df), "cov_RLE")}.
Paths to merged output files, Set to NULL to skip making merged covRle.}

\item{split.by.strand}{logical, default TRUE, split into forward and reverse
strand RleList inside covRle object.}

\item{seq_info}{SeqInfo object, default \code{seqinfo(findFa(df))}}

\item{weight}{integer, numeric or single length character. Default "score".
Use score column in loaded in_files.}

\item{format}{chatacter, default "qs", alternative "rds". File format to
save R object.}

\item{verbose}{logical, default TRUE, message about library output status.}
}
\value{
invisible(NULL), files saved to disc
}
\description{
Useful to store reads separated by readlength, for much faster
coverage calculation.
Saved by default in folder "cov_RLE_List" relative to default
libraries of experiment
}
\examples{
df <- ORFik.template.experiment()[10,]
## Usually do default folder, here we use tmpdir
folder_to_save <- file.path(tempdir(), "cov_RLE_List")
folder_to_save_merged <- file.path(tempdir(), "cov_RLE")
ORFik:::convert_to_covRleList(df, out_dir = folder_to_save,
out_dir_merged = folder_to_save_merged)
fimport(file.path(folder_to_save, "RFP_Mutant_rep2.covqs"))
}
\seealso{
Other lib_converters: 
\code{\link{convertLibs}()},
\code{\link{convert_bam_to_ofst}()},
\code{\link{convert_to_bigWig}()},
\code{\link{convert_to_covRle}()}
}
\concept{lib_converters}

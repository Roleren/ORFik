% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/experiment.R
\name{simpleLibs}
\alias{simpleLibs}
\title{Converted format of NGS libraries}
\usage{
simpleLibs(
  df,
  out.dir = libFolder(df),
  addScoreColumn = TRUE,
  addSizeColumn = TRUE,
  must.overlap = NULL,
  method = "None",
  type = "ofst",
  input.type = "ofst",
  reassign.when.saving = FALSE,
  envir = envExp(df),
  force = TRUE,
  library.names = bamVarName(df),
  libs = outputLibs(df, type = input.type, chrStyle = must.overlap, library.names =
    library.names, output.mode = "list", force = force, BPPARAM = BPPARAM),
  BPPARAM = bpparam()
)
}
\arguments{
\item{df}{an ORFik \code{\link{experiment}}}

\item{out.dir}{optional output directory, default: libFolder(df),
if it is NULL, it will just reassign R objects to simplified libraries.
Will then create a final folder specfied as: paste0(out.dir, "/", type, "/").
Here the files will be saved in format given by the type argument.}

\item{addScoreColumn}{logical, default TRUE, if FALSE will not add
replicate numbers as score column, see ORFik::convertToOneBasedRanges.}

\item{addSizeColumn}{logical, default TRUE, if FALSE will not add
size (width) as size column, see ORFik::convertToOneBasedRanges.
Does not apply for (GAlignment version of.ofst) or .bedoc. Since they
contain the original cigar.}

\item{must.overlap}{default (NULL), else a GRanges / GRangesList object, so
only reads that overlap (must.overlap) are kept. This is useful when you
only need the reads over transcript annotation or subset etc.}

\item{method}{character, default "None", the method to reduce ranges,
for more info see \code{\link{convertToOneBasedRanges}}}

\item{type}{character, output format, default "ofst".
Alternatives: "ofst", "bigWig", "wig","bedo" or "bedoc". Which format you want.
Will make a folder within out.dir with this name containing the files.}

\item{input.type}{character, input type "ofst". Remember this function
uses the loaded libraries if existing, so this argument is usually ignored.
Only used if files do not already exist.}

\item{reassign.when.saving}{logical, default FALSE. If TRUE, will reassign
library to converted form after saving. Ignored when out.dir = NULL.}

\item{envir}{environment to save to, default
\code{envExp(df)}, which defaults to .GlobalEnv, but can be set with
\code{envExp(df) <- new.env()} etc.}

\item{force}{logical, default TRUE If TRUE, reload library files even if
matching named variables are found in environment used by experiment
 (see \code{\link{envExp}}) A simple way to make
sure correct libraries are always loaded. FALSE is faster if data
is loaded correctly already.}

\item{library.names}{character vector, names of libraries, default:
name_decider(df, naming)}

\item{libs}{list, output of outputLibs as list of
GRanges/GAlignments/GAlignmentPairs objects. Set input.type and force arguments to define parameters.}

\item{BPPARAM}{how many cores/threads to use? default: bpparam().
To see number of threads used, do \code{bpparam()$workers}.
You can also add a time remaining bar, for a more detailed pipeline.}
}
\value{
invisible NULL (saves files to disc or R .GlobalEnv)
}
\description{
Export as either .ofst, .wig, .bigWig,.bedo (legacy format) or .bedoc (legacy format) files:\cr
Export files as .ofst for fastest load speed into R.\cr
Export files as .wig / bigWig for use in IGV or other genome browsers.\cr
The input files are checked if they exist from: \code{envExp(df)}.\cr
}
\details{
We advice you to not use this directly, as other function are more safe
for library type conversions. See family description below. This is
mostly used internally in ORFik. It is only adviced to use if large bam files
are already loaded in R and conversions are wanted from those.

See \code{\link{export.ofst}}, \code{\link{export.wiggle}},
\code{\link{export.bedo}} and \code{\link{export.bedoc}}
for information on file formats.\cr
If libraries of the experiment are
already loaded into environment (default: .globalEnv) is will export
using those files as templates. If they are not in environment the
.ofst files from the bam files are loaded (unless you are converting
to .ofst then the .bam files are loaded).
}
\examples{
df <- ORFik.template.experiment()
#convertLibs(df, out.dir = NULL)
# Keep only 5' ends of reads
#convertLibs(df, out.dir = NULL, method = "5prime")
}
\seealso{
Other lib_converters: 
\code{\link{convert_bam_to_ofst}()},
\code{\link{convert_to_bigWig}()},
\code{\link{convert_to_covRle}()},
\code{\link{convert_to_covRleList}()}
}

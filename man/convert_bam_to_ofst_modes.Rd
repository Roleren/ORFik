\name{convert_bam_to_ofst_modes}
\alias{convert_bam_to_ofst_modes}
\title{Convert one BAM sample to requested OFST mapper modes with one import}
\usage{
convert_bam_to_ofst_modes(df, modes = c("all", "unique"), out_dirs = NULL,
  read_strand = "auto", after_mode = NULL)
}
\arguments{
\item{df}{One-row ORFik experiment.}
\item{modes}{Character subset of \code{c("all", "unique")}.}
\item{out_dirs}{Named output directories for requested modes; defaults to experiment library folders.}
\item{read_strand}{RNA orientation selection; default \code{"auto"} uses validated source provenance.}
\item{after_mode}{Optional function(mode, path, distribution), called after publishing each output. The distribution is score-weighted; complex paired inputs provide NULL so the caller may import the output.}
}
\value{Invisible named paths to the completed output files.}
\description{Imports alignments once and filters NH tags before unique export. Atomic outputs, RNA strand normalization, selected-mate behaviour and weights follow the ordinary converter. Complex paired inputs use that converter as a fallback. Successful-mode callbacks allow independent resume flags.}
\seealso{\code{\link{convert_bam_to_ofst}}}

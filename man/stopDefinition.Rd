% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_ORFs.R
\name{stopDefinition}
\alias{stopDefinition}
\title{Returns stop codon definitions}
\usage{
stopDefinition(transl_table)
}
\arguments{
\item{transl_table}{numeric.  NCBI genetic code number for translation.}
}
\value{
A string of STOP sites separatd with "|".
}
\description{
According to:
<http://www.ncbi.nlm.nih.gov/Taxonomy/taxonomyhome.html/
index.cgi?chapter=tgencodes#SG1>
ncbi genetic code number for translation. This version is a cleaned up
version, unknown indices removed.
}
\examples{
stopDefinition
stopDefinition(1)

}
\seealso{
Other findORFs: 
\code{\link{findMapORFs}()},
\code{\link{findORFs}()},
\code{\link{findORFsFasta}()},
\code{\link{findUORFs}()},
\code{\link{startDefinition}()}
}
\concept{findORFs}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/txdbHelpers.R
\name{add_pseudo_5utrs_txdb_if_needed}
\alias{add_pseudo_5utrs_txdb_if_needed}
\title{add_pseudo_5utrs_txdb_if_needed}
\usage{
add_pseudo_5utrs_txdb_if_needed(
  txdb,
  pseudo_5UTRS_if_needed = NULL,
  minimum_5UTR_percentage = 30
)
}
\arguments{
\item{txdb}{a TxDb object}

\item{pseudo_5UTRS_if_needed}{integer, default NULL. If defined > 0,
will add pseudo 5' UTRs of maximum this length if 'minimum_5UTR_percentage" (default 30%) of
mRNAs (coding transcripts) do not have a leader. (NULL and 0 are both the ignore command)}

\item{minimum_5UTR_percentage}{numeric, default 30. What minimum percentage
of mRNAs most have a 5' UTRs (leaders), to not do the pseudo_UTR addition.
If percentage is higher, addition is ignored, set to 101 to always do it.}
}
\value{
txdb (new txdb if it was done, old if not)
}
\description{
add_pseudo_5utrs_txdb_if_needed
}

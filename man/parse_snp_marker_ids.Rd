% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{parse_snp_marker_ids}
\alias{parse_snp_marker_ids}
\title{Parse SNP marker IDs.}
\usage{
parse_snp_marker_ids(ids)
}
\arguments{
\item{ids}{Vector of SNP marker IDs.}
}
\value{
A \code{data.frame} with loci
corresponding to the SNP marker IDs.
}
\description{
This function parses an input vector of SNP marker IDs, and returns a \code{data.frame} with
the locus in each row derived from the corresponding marker ID in the input vector. An error
is raised if any of the input values cannot be parsed as a SNP marker ID.
}
\keyword{internal}

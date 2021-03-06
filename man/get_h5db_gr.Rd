% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/snakemake_functions.R
\name{get_h5db_gr}
\alias{get_h5db_gr}
\title{Get HDF5 database of DNAm fractions, from h5se rg data}
\usage{
get_h5db_gr(
  files.dpath = "recount-methylation-files",
  comp.dname = "compilations",
  ngsm.block = 50
)
}
\arguments{
\item{files.dpath}{Path to instance files directory.}

\item{comp.dname}{Name of compilations directory.}

\item{ngsm.block}{Number of GSM IDs per processed data block.}
}
\value{
NULL, generates an h5 gr database file as side effect.
}
\description{
Get HDF5 database of DNAm fractions, from h5se rg data
}

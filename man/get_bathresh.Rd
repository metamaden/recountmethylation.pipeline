% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bacontrol_format.R
\name{get_bathresh}
\alias{get_bathresh}
\title{Get BeadArray control outcomes from a matrix of metric signals}
\usage{
get_bathresh(rm, pass.lab = "PASS", fail.lab = "FAIL")
}
\arguments{
\item{rm}{BeadArray signals returned by `bactrl()` (matrix, rows = samples, 
cols = metrics).}

\item{pass.lab}{Label for passing samples ("PASS", chracter).}

\item{fail.lab}{Label for failing samples ("FAIL", chracter).}
}
\value{
Matrix of threshold assessments, either 'FAIL' or 'PASS', 
of same dimensions as input matrix.
}
\description{
Apply BeadArray thresholds to determine passing and failing samples. 
This function takes a matrix of signals and performs controls using 
minimum signal thresholds from the BeadArray Controls Reporter 
Software Guide ewastools.
}
\seealso{
bactrl
}

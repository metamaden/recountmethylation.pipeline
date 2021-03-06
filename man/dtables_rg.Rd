% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rmpipeline.R
\name{dtables_rg}
\alias{dtables_rg}
\title{Get red and grn signal data tables}
\usage{
dtables_rg(
  platform = c("hm450k", "epic"),
  version,
  ts,
  verbose = TRUE,
  gsmint = 60,
  overwrite = TRUE,
  fnstem = "mdat.compilation",
  sepval = " ",
  idatspath = file.path("recount-methylation-files", "idats"),
  destpath = file.path("recount-methylation-files", "compilations")
)
}
\arguments{
\item{platform}{The DNAm array platform of files to read/write (either 
"hm450k" or "epic").}

\item{version}{Version of the run for data table filenames.}

\item{ts}{NTP timestamp integer, for filenames (see get_metadata function).}

\item{verbose}{Whether to return verbose notifications.}

\item{gsmint}{Number of GSMs to process at a time, typically runs best near 50 samples.}

\item{overwrite}{Whether to overwrite existing data.table files with same destpath (default TRUE).}

\item{fnstem}{Filename stem for data tables.}

\item{sepval}{Separator symbol for data being written (default " ").}

\item{idatspath}{Path to IDAT files to read.}

\item{getnb}{Whether to get noob-normalized Beta-vlaues (default: FALSE).}
}
\value{

}
\description{
Generates 2-channel signal data tables. Validates IDATs, handles new paths 
and data.table options, and provides verbose input from main for loop to 
write data.
}
\examples{
#version = "0.0.1"
#timestamp = get_metadata("title", version)[["timestamp"]]
#dtables_rg(version = version, timestamp = timestamp)
}

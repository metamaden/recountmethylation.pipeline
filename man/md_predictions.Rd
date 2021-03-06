% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/md_dnam.R
\name{md_predictions}
\alias{md_predictions}
\title{Make model-based metadata predictions}
\usage{
md_predictions(
  ts,
  rgset.fname,
  grset.fname,
  mdmod.fname = "mdmod_dnam-predictions",
  files.dname = "recount-methylation-files",
  md.dname = "metadata",
  nsamp.block = 50,
  comp.dname = "compilations",
  verbose = TRUE
)
}
\arguments{
\item{ts}{Timestamp for the preprocessed metadata table to output 
(integer or character).}

\item{rgset.fname}{Name of the compilation file containing red signals 
extracted from IDATs (character).}

\item{grset.fname}{Name of the compilation file containing green signals 
extracted from IDATs (character).}

\item{mdmod.fname}{Name of the table of model-based predictions produced 
("mdmod_dnam-predictions").}

\item{files.dname}{Main recountmethylation instance files directory 
("recount-methylation-files").}

\item{md.dname}{Name of directory, in files.dname, containing the instance
metadata files ("metadata).}

\item{nsamp.block}{Number of samples per data block processed (integer, 50).}

\item{comp.dname}{Name of directory, in files.dname, containing the
compilation files with red and green signals for prediction calculations 
("compilations")}

\item{verbose}{Whether to show status messages (TRUE).}

\item{rgset.path}{Path}
}
\value{
NULL, produces table of model-based metadata predictions.
}
\description{
Make model-based metadata predictions using DNAm assays. Predictions for 
age, sex, and blood cell fractions are produced. Sex and blood cell 
predictions use raw/unnormalized red/green signal from h5se rg data 
(rgset.fname argument), while age predictions use noob-normalized 
Beta-values from h5se gr data (grset.fname argument).
}
\seealso{
get_qcmetrics(); md_preprocess(); md_postprocess()
}

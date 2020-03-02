#!/usr/bin/env R

# Run the Recount Methylation Pipeline to generate data objects for the R package.

require(rmpipeline)

# make data tables from IDAT files
# rmpipeline::h5db.fromsignal(version = commandArgs(T)[1])

vers <- "00.00.01"
rmd <- get_metadata("newrun", vers)
ts <- rmd[["timestamp"]]

dtables_fromsignal(version = "0.0.1", timestamp = , idats.path = "idats", dest.path = "")

# make an HDF5 databse from data tables
rmpipeline::make.h5db(dbn = "remethdb")

# make HDF5-SummarizedExperiment objects
rmpipeline::make.h5se()

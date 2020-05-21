#!/usr/bin/env R

# Run the Recount Methylation Pipeline to generate data objects for the R package.

require(rmpipeline)

# from main recount-methylation (base directory)
dtables_rg(version, timestamp)

# from recount-methylation/recount-methylation-analysis/files/mdata/compilations
makeh5db_rg(dbfnstem = "remethdb", version = "0.0.1", ts = 1589820348, 
            mdpath = "mdpost_all-gsm-md.rda", fnpath = ".",
            fnl = c("redsignal_1589820348_0-0-1.mdat.compilation",
                    "greensignal_1589820348_0-0-1.mdat.compilation"))

# make data tables from IDAT files
rmpipeline::h5db.fromsignal(version = commandArgs(T)[1])

# make an HDF5 databse from data tables
rmpipeline::make.h5db(dbn = "remethdb")

# make HDF5-SummarizedExperiment objects
rmpipeline::make.h5se()
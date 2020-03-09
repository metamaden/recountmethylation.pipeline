#!/usr/bin/env R

# Run the Recount Methylation Pipeline to generate data objects for the R package.

require(rmpipeline)

# make data tables from IDAT files
# rmpipeline::h5db.fromsignal(version = commandArgs(T)[1])

vers <- "00.00.01"
rmd <- get_metadata("newrun", vers)
ts <- rmd[["timestamp"]]

dtables_fromsignal(version = "0.0.1", timestamp = ts,
                   idatspath = "idats", destpath = "compilations")

# make an HDF5 databse from data tables
make_h5db(dbfnstem = "remethdb",
          version = "0.0.1", ts = ts,
          fnl = list.files("compilations"),
          addmd = TRUE, mdpath = "mdpost_all-gsm-md.rda",
          fnpath = "compilations", rmax = 2)

# make HDF5-SummarizedExperiment objects
make_h5se("remethdb-seh5", "0.0.1", "1123", se = "rg",
          dbn = "remethdb_1123_0-0-1.h5",
          dsn.data1 = "redsignal", dsn.data2 = "greensignal",
          dsn.rn = "redsignal.rownames", addpheno = TRUE, dsn.md = "mdpost",
          dsn.cn = "redsignal.colnames")




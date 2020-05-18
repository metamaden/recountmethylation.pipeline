#!/usr/bin/env R

# Run the Recount Methylation Pipeline to generate data objects for the R package.

require(rmpipeline)

# make data tables from IDAT files
# rmpipeline::h5db.fromsignal(version = commandArgs(T)[1])

vers <- "00.00.01"
rmd <- get_metadata("newrun", vers)
ts <- rmd[["timestamp"]]

# run this from recount-methylation main dir (e.g. `cd recount-methylation`)
dtables_fromsignal(version = "0.0.1", timestamp = ts, getnb = FALSE)

# make an HDF5 databse from data tables
make_h5db(dbfnstem = "remethdb",
          version = "0.0.1", ts = ts,
          fnl = list.files("compilations"),
          addmd = TRUE, mdpath = "mdpost_all-gsm-md.rda",
          fnpath = "compilations", rmax = 2)

# make HDF5-SummarizedExperiment objects

# RGset
make_h5se("remethdb-seh5", "0.0.1", "1123", se = "rg",
          dbn = "remethdb_1123_0-0-1.h5",
          dsn.data1 = "redsignal", dsn.data2 = "greensignal",
          dsn.rn = "redsignal.rownames", addpheno = TRUE, dsn.md = "mdpost",
          dsn.cn = "redsignal.colnames")

# GRset
make_h5se("remethdb-seh5", "0.0.1", "1123", se = "gr",
          dbn = "remethdb_1123_0-0-1.h5",
          dsn.data1 = "redsignal", dsn.data2 = "greensignal",
          dsn.rn = "redsignal.rownames", addpheno = TRUE, dsn.md = "mdpost",
          dsn.cn = "redsignal.colnames")

# GMset
make_h5se("remethdb-seh5", "0.0.1", "1123", se = "gm",
          dbn = "remethdb_1123_0-0-1.h5",
          dsn.data1 = "redsignal", dsn.data2 = "greensignal",
          dsn.rn = "redsignal.rownames", addpheno = TRUE, dsn.md = "mdpost",
          dsn.cn = "redsignal.colnames")

#-------------------
# debugging at scale
#-------------------
library(rmpipeline)

vers <- "0.0.1"
rmd <- get_metadata("newrun", vers)
ts <- rmd[["timestamp"]]

idatspath <- "recount-methylation-files/idats"
fnpath = "recount-methylation-analysis/files/mdata/compilations"

dtables_fromsignal(version = "0.0.1", timestamp = ts,
                   idatspath = idatspath, destpath = fnpath)

fnl = c("redsignal_1583780004_0-0-1.mdat.compilation", 
        "greensignal_1583780004_0-0-1.mdat.compilation")

make_h5db(dbfnstem = "remethdb",
          version = "0.0.1", ts = ts,
          fnl = fnl, addmd = TRUE, mdpath = "mdpost_all-gsm-md.rda",
          fnpath = fnpath, rmax = 1000)





#!/usr/bin/env R

# Run the Recount Methylation Pipeline to generate data objects for the R package.

require(rmpipeline)

# make data tables from IDAT files
# rmpipeline::h5db.fromsignal(version = commandArgs(T)[1])

vers <- "00.00.01"
rmd <- get_metadata("newrun", vers)
ts <- rmd[["timestamp"]] # 1583780004
semd <- list(version = vers, timestamp = ts)

dtables_fromsignal(version = "0.0.1", timestamp = ts,
                   idatspath = "idats", destpath = "compilations")

# make an HDF5 databse from data tables
make_h5db(dbfnstem = "remethdb",
          version = "0.0.1", ts = ts,
          fnl = list.files("compilations"),
          addmd = TRUE, mdpath = "mdpost_all-gsm-md.rda",
          fnpath = "compilations", rmax = 2)

# make HDF5-SummarizedExperiment objects
dbn <- "remethdb2.h5"
# make rg h5se
make_h5se(dbn = dbn, newfnstem = "remethdb_h5se_rg", version = vers,
          ts = ts, se = "rg", dsn.data1 = "redsignal", dsn.data2 = "greensignal",
          addpheno = TRUE, dsn.md = "mdpost", dsn.rn = "redsignal.rownames",
          dsn.cn = "redsignal.colnames")

# make gr h5se
make_h5se(dbn = dbn, newfnstem = "remethdb_h5se_gr", version = vers,
          ts = ts, se = "gr", dsn.data1 = "noobbeta", addpheno = TRUE,
          phenopath = "mdpost_all-gsm-md.rda",
          dsn.md = "mdpost", dsn.rn = "redsignal.rownames", dsn.cn = "redsignal.colnames")

# make gm h5se
make_h5se(dbn = dbn, newfnstem = "remethdb_h5se_gm",
          version = vers, ts = ts, se = "gr",
          dsn.data1 = "methylated_signal", dsn.data2 = "unmethylated_signal",
          addpheno = TRUE, phenopath = "mdpost_all-gsm-md.rda",
          dsn.md = "mdpost", dsn.rn = "redsignal.rownames", dsn.cn = "redsignal.colnames")


#---------

library(rmpipeline)

vers <- "0.0.1"
rmd <- get_metadata("newrun", vers)
ts <- rmd[["timestamp"]]

dtables_fromsignal(version = "0.0.1", timestamp = ts,
                   idatspath = "recount-methylation-files/idats",
                   destpath = "recount-methylation-analysis/files/mdata/compilations")


dtables_fromsignal2(version = "0.0.1", timestamp = ts)

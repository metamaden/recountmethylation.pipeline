#!/usr/bin/env R

library(rmpipeline)
library(rhdf5)

fn <- "remethdb_h5se-test_gr_00-00-01_1583780004"


fnl = c("newred.comp", "1553728068.rawgrn.compilation.mdat", "1553728595.noobbeta.compilation.mdat")
dsnl = c("redsignal", "greensignal", "noobbeta")
dbn = "remethdb_test.h5"
h5createFile(dbn)
rmax = 10; cmax = 10 # 622399

for(d in dsnl){
  h5delete(dbn, d)
  h5delete(dbn, paste0(d, ".colnames"))
  h5delete(dbn, paste0(d, ".rownames"))
}









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
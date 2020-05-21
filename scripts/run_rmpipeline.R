#!/usr/bin/env R

require(rmpipeline)
require(recountmethylation)

# Run the Recount Methylation Pipeline to generate data objects for the R package.
#
# Purpose: version dataset files and coerce data into usable data objects.
# 
# Deliverables: 
# * 2 signal tables (red and green channels)  
# * 1 h5 database file containing signal tables and samples metadata  
# * 3 h5se objects of different HDF5-SummarizedExperiment object classes:
# ** 1 unnormalized RGChannelSet
# ** 1 unnormalized GenomicMethylSet
# ** 1 noob-normalized GenomicRangesSet

#--------------------------
# get new datasets metadata
#--------------------------
# filenames example
# md[["timestamp"]] # [1] "1590090412"
# version <- "0.0.1"
# ts <- 1590090412
# dbn <- "remethdb_1590090412_0-0-1"

md <- get_metadata(title = "newrun", version = "0.0.1")

#--------------------------
# make the rg signal tables
#--------------------------
# navigate to main recount-methylation dir/base dir.
# e.g. 
# > cd recount-methylation

dtables_rg(version, timestamp)

#-----------------
# make the h5 file
#-----------------
# navigate to compilations dir
# e.g. 
# > cd recount-methylation/recount-methylation-analysis/files/mdata/compilations

makeh5db_rg(dbfnstem = "remethdb", version = "0.0.1", ts = 1589820348, 
            mdpath = "mdpost_all-gsm-md.rda", fnpath = ".",
            fnl = c("redsignal_1589820348_0-0-1.mdat.compilation",
                    "greensignal_1589820348_0-0-1.mdat.compilation"))

#------------------
# make h5se objects
#------------------
# examples arguments:
# dbn <- "remethdb_1590090412_0-0-1.h5"
# fnstem <- "remethdb_h5se-rg"

rmpipeline::make_h5se(dbn = dbn, newfnstem = fnstem, 
                      dsn.data1 = "redsignal", dsn.data2 = "greensignal",
                      version = version, ts = ts, se = "rg")

# make gm and gr objects from rg h5se






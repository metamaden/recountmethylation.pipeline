#!/usr/bin/env R

# Author: Sean Maden
# 
# Title:
# Run the Recount Methylation Pipeline to generate data objects for the R package.
#
# Purpose:
# This script outputs versioned, timestamped database files from a set of GEO IDATs.
# It is assumed the user has run recount-methylation-server and has a structured 
# directory tree containing the sample IDAT files and the sample postprocessed metadata.
# From these files, several database files are generated. 
# 
# First, all database files for raw red and green signals are created. This approach 
# expedites database creation so full data objects can be viewed as quickly as possible.
# First, large flat tables of raw/unnormalized red and green signals are created. 
# Second, an h5 HDF5 database containing the red/grn signals and sample metadata is created
# with blocking. Third, an h5se RGChannelSet is created. Note the h5 and h5se files are
# made available at recount.bio for the recountmethylation R package.
# 
# Next, h5 and h5se objects are made for raw meth/unmeth signals and noob-normalized 
# Beta-values. These are also made available at recount.bio/data for recountmethylation. 
# First, h5 files are generated from the h5 red/grn data in blocks. Second, the h5se 
# objects are created from the corresponding h5 objects with process delay (e.g. data 
# chunk management is performed silently).
# 
# Output: 
# 1. Signal tables (raw, red/grn signal)
# 2. 3 h5 files (raw red/grn signal, raw meth/unmeth signal, noob-norm beta-values)
# 3. 3 h5se objects (raw RGChannelSet, raw MethylSet, raw RatioSet)

library(rmpipeline)

#------------------
# datasets metadata
#------------------
# set new version
newversion <- "0.0.2"
# get new datasets metadata
md <- get_metadata(title = "newrun", version = newversion)
versionfn <- md[["version"]]; timestamp <- md[["timestamp"]]


idats.path <- file.path("recount-methylation-files", "idats")
idatsv <- dt_checkidat(idatspath=idats.path, verbose = TRUE)

#--------------------------
# filter idats on file size
#--------------------------
dirpath <- file.path("recount-methylation-files", "idats")
fnv <- list.files(dir.path)
fnv <- fnv[grepl(".*idat$", fnv) & grepl(".*hlink.*", fnv)]
dat <- file.info(file.path(dirpath, fnv[1]))
fnvf <- fnv[2:length(fnv)]
for(fi in seq(fnvf)){
  dat <- rbind(dat, file.info(file.path(dirpath, fnvf[fi])))
  message(fi)
}
datf <- dat[dat$size<1.2e7,]; dim(datf)
for(r in rownames(datf)){file.remove(r)}
#file.exists(file.path(dirpath, rownames(datf)[1]))
#rownames(dat)

#----------------------
# red/grn signal data
#----------------------
# navigate to main recount-methylation dir/base dir.
# e.g. 
# > cd recount-methylation
dtables_rg_epic(versionfn, timestamp, destpath = "compilations")

# make the h5 file
# navigate to compilations dir
# e.g. 
# > cd recount-methylation/recount-methylation-analysis/files/mdata/compilations
make_h5db_rg(dbfnstem = "remethdb", version = versionfn, ts = timestamp, 
            mdpath = "mdpost_all-gsm-md.rda", fnpath = ".",
            fnl = c("redsignal_1589820348_0-0-1.mdat.compilation",
                    "greensignal_1589820348_0-0-1.mdat.compilation"))



fnl <- c("greensignal_1606324405_0-0-2.mdat.compilation", 
    "redsignal_1606324405_0-0-2.mdat.compilation")

makeh5db_rg(dbfnstem = "remethdb", version = "0.0.2", ts = "1589820348", 
            mdpath = "mdpost_all-gsm-md.rda", fnpath = ".", fnl = fnl)

# make the h5se file
make_h5se(dbn = dbn, newfnstem = fnstem, version = versionfn, ts = timestamp)

#------------------------------
# meth/unmeth and betavals data
#------------------------------
library(rmpipeline)

ts <- "1589820348"; version <- "0.0.2"
dbn <- paste0("remethdb_", ts, "_", gsub("\\.", "-", version), ".h5")
# dbn <- "remethdb_1589820348_0-0-2.h5"


make_h5se_rg(dbn = dbn, ts = ts, version = version, cmax = 1052641,
    rmax = 14000)



# make new h5 files
h5name.gm <- make_h5_gm(dbn = dbn, ts = ts, version = version)
h5name.gr <- make_h5_gr()
# make new h5se objects
make_h5se_gm(h5name.gm)
make_h5se_gr(h5name.gr)






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

library(rmpipeline)

read.path <- file.path("home","metamaden","recount-methylation-epic","compilations")
write.path <- file.path("eternity","recount-methylation","recount-methylation-epic")

fnl <- c("greensignal_1606324405_0-0-2.mdat.compilation",
    "redsignal_1606324405_0-0-2.mdat.compilation")

make_h5db_rg(dbfnstem = "remethdb", dbpath = write.path, version = "0.0.2", 
    ts = "1589820348", mdpath = "mdpost_all-gsm-md.rda", fnpath = read.path, fnl = fnl,
    ngsm.block = 50, cmax = 1052641, rmax = 14000)

fnl <- "greensignal_1606324405_0-0-2.mdat.compilation"
make_h5db_rg(dbfnstem = "remethdb", dbpath = write.path, version = "0.0.2", 
    ts = "1589820348", mdpath = "mdpost_all-gsm-md.rda", fnpath=read.path, fnl=fnl,
    ngsm.block = 50, cmax = 1052641, rmax = 14000, dsnl = "greensignal")



platform = "epic" # c("hm450k", "epic")
version = "0.0.2"
ts = 1589820348
dbn = "remethdb_1589820348_0-0-2.h5"
newfnstem = "remethdb_h5se-rg"
dsnv = c("redsignal", "greensignal")
add.metadata=FALSE
mdpath=NULL
dsn.md="mdpost"
dsn.md.cn=paste0(dsn.md,".colnames")
verbose = TRUE
replace.opt = TRUE
dsn.rnv = c(paste0(dsnv[1], ".rownames"), paste0(dsnv[2], ".rownames"))
dsn.cnv = c(paste0(dsnv[1], ".colnames"), paste0(dsnv[2], ".colnames"))
semd=list("title"="RGChannelSet HDF5-SummarizedExperiment object",
    "preprocessing"="raw")

#--------------
# rgchannel set
#--------------
# hdf5 db
fnl <- "redsignal_1606324405_0-0-2.mdat.compilation"
make_h5db_rg(dbfnstem = "remethdb", dbpath = write.path, version = "0.0.2", 
    ts = "1589820348", mdpath = "mdpost_all-gsm-md.rda", fnpath = read.path, 
    fnl = fnl, ngsm.block = 50, cmax = 1052641, rmax = 14000, dsnl="redsignal")
fnl <- "greensignal_1606324405_0-0-2.mdat.compilation"
make_h5db_rg(dbfnstem = "remethdb", dbpath = write.path, version = "0.0.2", 
    ts = "1589820348", mdpath = "mdpost_all-gsm-md.rda", fnpath=read.path, 
    fnl=fnl, ngsm.block = 50, cmax = 1052641, rmax = 14000, dsnl = "greensignal")

# h5se file
make_h5se_rg(max.sample = 12650, platform = "epic", version = "0.0.2", 
    ts = 1589820348, dbn = "remethdb_1589820348_0-0-2.h5", 
    newfnstem = "remethdb_h5se-rg")

#-------------------
# genomic methyl set
#-------------------
# hdf5 db
make_h5_gm(dbn = "remethdb_h5se-rg_epic_0-0-2_1589820348", version = "0.0.2", 
    ts = 1589820348, num.samp = 12650, blocksize = 65, platform = "epic", 
    newfnstem = "remethdb_h5-gm", verbose = TRUE, replace.opt = TRUE)

# h5se file
make_h5se_gm(dbn = "remethdb_h5-gm_epic_0-0-2_1589820348.h5", version = "0.0.2", 
    ts = 1589820348, platform = "epic", replaceopt = TRUE, verbose = TRUE, 
    add.metadata = FALSE, pdata = NULL, newdnstem = "remethdb_h5se-gm")

#------------------
# genomic methylset
#------------------
# hdf5 db
make_h5_gr(dbn = "remethdb_h5se-rg_epic_0-0-2_1589820348", version = "0.0.2", 
    ts = 1589820348, num.samp = 12650, blocksize = 10, platform = "epic", 
    newfnstem = "remethdb_h5-gr", verbose = TRUE, replace.opt = TRUE)

# h5se file
make_h5se_gr(dbn = "remethdb_h5-gr_epic_0-0-2_1589820348.h5", version = "0.0.2", 
    ts = 1589820348, platform = "epic", replaceopt = TRUE, 
    verbose = TRUE, add.metadata = FALSE, pdata = NULL, 
    newdnstem = "remethdb_h5se-gr", 
  semd=list("title"="GenomicMethylSet HDF5-SummarizedExperiment object",
    "preprocessing"="Normalization with out-of-band signal (noob)"))



#--------------
# HM450K arrays
#--------------
library(rmpipeline)

# datasets metadata
platform <- "hm450k"
newversion <- "0.0.2"
md <- get_metadata(title = "newrun", version = newversion)
versionfn <- md[["version"]]; timestamp <- md[["timestamp"]]

timestamp <- ts <- 1607018051
version <- "0.0.2"

# idats paths
idats.path <- file.path("home", "metamaden", "recount-methylation-hm450k", 
    "recount-methylation-files", "idats")
idatsv <- dt_checkidat(idatspath=idats.path, verbose = TRUE)

# filter idats on file size
#dirpath <- file.path("recount-methylation-files", "idats")
#fnv <- list.files(dir.path)
#fnv <- fnv[grepl(".*idat$", fnv) & grepl(".*hlink.*", fnv)]
#dat <- file.info(file.path(dirpath, fnv[1]))
#fnvf <- fnv[2:length(fnv)]
#for(fi in seq(fnvf)){
#  dat <- rbind(dat, file.info(file.path(dirpath, fnvf[fi])))
#  message(fi)
#}
#datf <- dat[dat$size<1.2e7,]; dim(datf)
#for(r in rownames(datf)){file.remove(r)}
#file.exists(file.path(dirpath, rownames(datf)[1]))
#rownames(dat)

# red/grn signal data
writepath <- file.path("eternity", "recount-methylation", 
    "recount-methylation-hm450k")
readpath <- file.path("home", "metamaden", "recount-methylation-hm450k", 
    "recount-methylation-files", "idats")

dtables_rg(platform = platform, version = version, timestamp = timestamp, 
    idatspath = readpath, destpath = writepath)








dtables_fromsignal(version = "0.0.1", timestamp = ts, getnb = FALSE)

#----------------------
# dtables function, new
#----------------------

dtables_rg()

dtables_rg <- function(version, timestamp, verbose = TRUE, gsmint = 60,
                               overwrite = TRUE, fnstem = "mdat.compilation", sepval = " ",
                               idatspath = paste("recount-methylation-files", "idats", sep = "/"),
                               destpath = paste("recount-methylation-analysis",
                                                "files", "mdata", "compilations", sep = "/")){
  # get valid gsms idats dir
  idatinfo <- dt_checkidat(idatspath = idatspath, verbose = verbose)
  gpath <- idatinfo[["gpath"]]
  gsmu <- idatinfo[["gsmu"]]
  if(verbose){message("Found ", length(gsmu), " GSM IDs with valid idats.")}
  gsmii <- getblocks(length(gsmu), gsmint) # get list of id vectors, e.g. "blocks"
  if(verbose){message("Making new data tables")}
  dtinfo <- dt_makefiles(gpath = gpath, idatspath = idatspath, 
                         destpath = destpath, version = version, 
                         nts = timestamp, overwrite = overwrite,
                         verbose = verbose)
  dtcond <- dtinfo[["dtcond"]]
  if(dtcond){
    reds.path <- dtinfo[["reds.path"]]
    grns.path <- dtinfo[["grns.path"]]
    if(verbose){message("Appending new table data...")}
    tt <- Sys.time()
    for(i in 1:length(gsmii)){
      dt_write_rg(gi = gsmii[i], idatspath = idatspath,
                  gpath = gpath, reds.path = reds.path,
                  grns.path = grns.path, tt = tt, verbose = verbose)}
  } else{stop("Problem encountered handling data tables.")}
  if(verbose){message("Successfully completed data tables.")}
  return(NULL)
}

# test the new function
version <- "0.0.1"
timestamp <- 1589762825
verbose = TRUE
gsmint = 60
overwrite = TRUE
fnstem = "mdat.compilation"
sepval = " "
idatspath = paste("recount-methylation-files", "idats", sep = "/")
destpath = paste("recount-methylation-analysis",
                 "files", "mdata", "compilations", sep = "/")


# get valid gsms idats dir
idatinfo <- dt_checkidat(idatspath = idatspath, verbose = verbose)
gpath <- idatinfo[["gpath"]]
gsmu <- idatinfo[["gsmu"]]
if(verbose){message("Found ", length(gsmu), " GSM IDs with valid idats.")}
gsmii <- getblocks(length(gsmu), gsmint) # get list of id vectors, e.g. "blocks"
if(verbose){message("Making new data tables")}
dtinfo <- dt_makefiles(gpath = gpath, idatspath = idatspath, 
                       destpath = destpath, version = version, 
                       nts = timestamp, overwrite = overwrite,
                       verbose = verbose)
dtcond <- dtinfo[["dtcond"]]
if(dtcond){
  reds.path <- dtinfo[["reds.path"]]
  grns.path <- dtinfo[["grns.path"]]
  if(verbose){message("Appending new table data...")}
  tt <- Sys.time()
  for(i in 1:length(gsmii)){
    dt_write_rg(gi = gsmii[i], idatspath = idatspath, gpath = gpath, tt = tt,
                reds.path = reds.path, grns.path = grns.path, verbose = verbose)
    }
} else{stop("Problem encountered handling data tables.")}
if(verbose){message("Successfully completed data tables.")}


#------------------------------------------------
# function to handle making new data.table files
#------------------------------------------------
dt_makefiles <- function(gpath, idatspath, destpath, version, nts, 
                         overwrite = TRUE, fnstem = "mdat.compilation",
                         verbose = TRUE){
  version.fn <- gsub("\\.", "-", version)
  reds.fn <- paste("redsignal", nts, version.fn, sep = "_")
  reds.fn <- paste(reds.fn, fnstem, sep = ".")
  grns.fn <- paste("greensignal", nts, version.fn, sep = "_")
  grns.fn <- paste(grns.fn, fnstem, sep = ".")
  reds.path = paste(destpath, reds.fn, sep = "/")
  grns.path = paste(destpath, grns.fn, sep = "/")
  cn = c("gsmi")
  rgi = minfi::read.metharray(c(paste(idatspath, gpath[1:2], sep = "/")))
  rgcni = colnames(t(minfi::getRed(rgi))) # cgids as colnames
  rgcn = matrix(c(cn, rgcni), nrow = 1)
  if(overwrite){
    if(verbose){message("Making/verifying data tables...")}
    dt1 <- try(data.table::fwrite(rgcn, reds.path, sep = sepval, 
                                  append = FALSE, col.names = F))
    dt2 <- try(data.table::fwrite(rgcn, grns.path, sep = sepval, 
                                  append = FALSE, col.names = F))
    dtcond <- is.null(dt1) & is.null(dt2)
  } else{
    dt1 <- try(file.exists(reds.path))
    dt2 <- try(file.exists(grns.path))
    dtcond <- dt1 == TRUE & dt1 == TRUE
  }
  lr <- list("dtcond" = dtcond,
             "reds.path" = reds.path,
             "grns.path" = grns.path)
  return(lr)
}

#----------------------------------
# making a new idats check function
#----------------------------------

dt_checkidat <- function(idatspath, verbose = TRUE){
  # checks for valid idat file pairs
  # checks for latest timestamps, files with matching timestamps,
  # and valid red/grn IDAT file pairs
  # get valid gsms idats dir
  if(verbose){message("Getting valid GSM IDs from IDAT filenames...")}
  idats.lf = list.files(idatspath)
  which.valid1 = grepl("\\.idat$", substr(idats.lf,
                                          nchar(idats.lf) - 4,
                                          nchar(idats.lf))) # idat ext
  which.valid2 = grepl(".*hlink.*", idats.lf) # hlink
  # timestamp
  which.valid3 = length(gsub("\\..*", "", gsub(".*\\.", "", idats.lf))) > 0 
  idats.valid = idats.lf[which.valid1 & which.valid2 & which.valid3]
  gsmu = gsub("\\..*", "", idats.valid)
  gsmu = gsmu[grepl("^GSM.*", gsmu)]
  gsmu = unique(gsmu)
  # red/grn channel files
  if(verbose){message("Finding paired red and grn channel files...")}
  gstr <- gsub(".*hlink\\.", "", gsub("(_Red.idat$|_Grn.idat$)", "", idats.valid))
  rsub <- idats.valid[grepl(".*_Red.idat$", idats.valid)]
  gsub <- idats.valid[grepl(".*_Grn.idat$", idats.valid)]
  if(verbose){message("Intersecting .hlink and _Red.idat or _Grn.idat")}
  gstr.rg <- intersect(gsub(".*hlink\\.", "", gsub("_Red.idat$", "", rsub)),
                       gsub(".*hlink\\.", "", gsub("_Grn.idat$", "", gsub)))
  idats.validf <- idats.valid[gstr %in% gstr.rg]
  gpath <- unique(gsub("(_Red.idat|_Grn.idat)", "", idats.validf))
  # filter on timestamps
  if(verbose){message("Filtering files on latest timestamps...")}
  gpath.ts <- c()
  gpath.gid <- c()
  gidv <- gsub("\\..*", "", gpath)
  tidv <- gsub(".*\\.", "", gsub("\\.hlink.*", "", gpath))
  # iterate on timestamps
  tsu <- unique(tidv)
  tsu <- tsu[rev(order(tsu))] # sort so that latest comes first
  for(t in tsu){
    gpathf <- gpath[tidv == t & !gidv %in% gpath.gid]
    if(verbose){message("Removing redundant gsm ids in this iter")}
    gpathf.gid <- gsub("\\..*", "", gpathf)
    gpathf <- gpathf[!duplicated(gpathf.gid)]
    if(verbose){message("Appending new fn")}
    gpath.ts <- c(gpath.ts, gpathf)
    if(verbose){message("remaking gid list")}
    gpath.gid <- gsub("\\..*", "", gpath.ts)
  }
  gpath <- gpath.ts
  gsmu <- gsub("\\..*", "", gpath)
  lr <- list("gsmu" = gsmu, "gpath" = gpath)
  return(lr)
}

#----------------------------
# making a new write function
#----------------------------

dt_write_rg <- function(gi, idatspath, gpath, reds.path, grns.path, tt, 
                        sepval = " ", verbose = TRUE){
  # gi a vector of gsm ids
  # idatspath a path of idat files
  # Notes: currently just skips chunks that throw read.metharray error
  gi = gsmii[[i]]
  if(verbose){message("Reading in new data")}
  pathl = paste(idatspath, gpath[gi], sep = "/")
  rgi = try(minfi::read.metharray(c(pathl)))
  if(!class(rgi) == "RGChannelSet"){
    message("There was a problem reading data chunk num. ", i)
  } else{
    if(verbose){message("getting data matrices")}
    redi = matrix(c(colnames(rgi), 
                    t(minfi::getRed(rgi))), ncol = nrow(rgi) + 1)
    grni = matrix(c(colnames(rgi), 
                    t(minfi::getGreen(rgi))), ncol = nrow(rgi) + 1)
    if(verbose){message("appending new data")}
    data.table::fwrite(redi, reds.path, sep = sepval, append = TRUE)
    data.table::fwrite(grni, grns.path, sep = sepval, append = TRUE)
  }
  if(verbose){
    te <- Sys.time() - tt
    message("Finished chunk ", i ," / ", 
            length(gsmii), ", time elapsed: ", te)
  }
  return(NULL)
}

#-------------------------------------------
# debugging readBin error for read_metharray
#-------------------------------------------
# notes about error
# (prints: 20281 to 20340, time elapsed: 8.5225774894158, then err)
# “Error in readBin(con, what = "integer", n = n, size = 4, endian = "little",  :
#                     invalid 'n' argument”

version = "0.0.1"
timestamp = ts
getnb = FALSE

verbose = TRUE
gsmint = 60
fnstem = "mdat.compilation"
sepval = " "
idatspath = paste("recount-methylation-files", "idats", sep = "/")
destpath = paste("recount-methylation-analysis",
                 "files", "mdata", "compilations", sep = "/")

# get run metadata
#if(is.null(timestamp)){
#  runmd <- get_metadata(version, "notitle")
#  nts = runmd[["timestamp"]]
#} else{
  nts <- timestamp
#}

# get valid gsms idats dir
#if(verbose){message("Getting valid GSM IDs from IDAT filenames...")}
idats.lf = list.files(idatspath)
which.valid1 = grepl("\\.idat$", substr(idats.lf,
                                        nchar(idats.lf) - 4,
                                        nchar(idats.lf))) # idat ext
which.valid2 = grepl(".*hlink.*", idats.lf) # hlink
which.valid3 = length(gsub("\\..*", "",
                           gsub(".*\\.", "", idats.lf))) > 0 # timestamp
idats.valid = idats.lf[which.valid1 &
                         which.valid2 &
                         which.valid3]
gsmu = gsub("\\..*", "", idats.valid)
gsmu = gsmu[grepl("^GSM.*", gsmu)]
gsmu = unique(gsmu)
# red/grn channel files
#if(verbose){message("Finding paired red and grn channel files...")}
gstr <- gsub(".*hlink\\.", "", gsub("(_Red.idat$|_Grn.idat$)", "", idats.valid))
rsub <- idats.valid[grepl(".*_Red.idat$", idats.valid)]
gsub <- idats.valid[grepl(".*_Grn.idat$", idats.valid)]
# intersect strings between ".*hlink\\." and "(_Red.idat|_Grn.idat)"
gstr.rg <- intersect(gsub(".*hlink\\.", "", gsub("_Red.idat$", "", rsub)),
                     gsub(".*hlink\\.", "", gsub("_Grn.idat$", "", gsub)))
idats.validf <- idats.valid[gstr %in% gstr.rg]
gpath <- unique(gsub("(_Red.idat|_Grn.idat)", "", idats.validf))
# filter on timestamps
#if(verbose){message("Filtering files on latest timestamps...")}
gpath.ts <- c()
gpath.gid <- c()
gidv <- gsub("\\..*", "", gpath)
tidv <- gsub(".*\\.", "", gsub("\\.hlink.*", "", gpath))
# iterate on timestamps
tsu <- unique(tidv)
tsu <- tsu[rev(order(tsu))] # sort so that latest comes first
for(t in tsu){
  gpathf <- gpath[tidv == t & !gidv %in% gpath.gid]
  # rm redundant gsm ids in this iter
  gpathf.gid <- gsub("\\..*", "", gpathf)
  gpathf <- gpathf[!duplicated(gpathf.gid)]
  # append new fn
  gpath.ts <- c(gpath.ts, gpathf)
  # remake gid list
  gpath.gid <- gsub("\\..*", "", gpath.ts)
}
gpath <- gpath.ts
gsmu <- gsub("\\..*", "", gpath)
#if(verbose){message("Finished filtering filenames. ",
#                    "Found ", length(gsmu),
#                    " valid file pairs.")
#}
# get GSM ID indices by interval, as list
gsmii = getblocks(length(gsmu), gsmint)
# make new data paths
#if(verbose){message("Making and instantiating new data table files...")}
reds.fn <- paste(paste("redsignal", nts,
                       gsub("\\.", "-", version),
                       sep = "_"),
                 fnstem, sep = ".") # "redsignal_1589762825_0-0-1.mdat.compilation"
grns.fn <- paste(paste("greensignal", nts,
                       gsub("\\.", "-", version),
                       sep = "_"),
                 fnstem, sep = ".")
reds.path = paste(destpath, reds.fn, sep = "/")
grns.path = paste(destpath, grns.fn, sep = "/")
#if(getnb){
#  nb.fn <- paste(paste("noobbeta", nts,
#                       gsub("\\.", "-", version),
#                       sep = "_"),
#                 fnstem, sep = ".")
#  nb.path = paste(destpath, nb.fn, sep = "/")
#}
# instantiate new empty data tables with probes as colnames
cn = c("gsmi")
rgi = minfi::read.metharray(c(paste(idatspath, gpath[1:2], sep = "/")))
rgcni = colnames(t(minfi::getRed(rgi))) # grcni = colnames(t(getGreen(rgi)))
# rgcni == grcni
rgcn = matrix(c(cn, rgcni), nrow = 1)
data.table::fwrite(rgcn, reds.path, sep = sepval, append = FALSE, col.names = F)
data.table::fwrite(rgcn, grns.path, sep = sepval, append = FALSE, col.names = F)
#if(getnb){
#  if(verbose){message("Instantiating noob beta table...")}
#  nbi = minfi::preprocessNoob(rgi)
#  nbcn = matrix(c(cn, colnames(t(getBeta(nbi)))), nrow = 1)
#  data.table::fwrite(nbcn, nb.path, sep = sepval, append = FALSE, col.names = F)
#}
# append new methdata
tt = Sys.time()
for(i in 1:length(gsmii)){
  gi = gsmii[[i]]
  # read in new data
  pathl = paste(idatspath, gpath[gi], sep = "/")
  # rgi = minfi::read.metharray(c(pathl))
  rgi = try(minfi::read.metharray(c(pathl)))
  if(!class(rgi) == "RGChannelSet"){
      message("There was a problem reading data chunk num. ", i)
  } else{
    if(verbose){message("getting data matrices")}
    redi = matrix(c(colnames(rgi), 
                    t(minfi::getRed(rgi))), ncol = nrow(rgi) + 1)
    grni = matrix(c(colnames(rgi), 
                    t(minfi::getGreen(rgi))), ncol = nrow(rgi) + 1)
    
    if(verbose){message("appending new data")}
    data.table::fwrite(redi, reds.path, sep = sepval, append = TRUE)
    data.table::fwrite(grni, grns.path, sep = sepval, append = TRUE)
    
    # parse noob-normalized data option
    if(getnb){
      gsi = minfi::preprocessNoob(rgi)
      nbi = matrix(c(colnames(gsi), t(minfi::getBeta(rgi))), ncol = nrow(gsi) + 1)
      data.table::fwrite(nbi, nb.path, sep = sepval, append = TRUE)
    }
    }
  if(verbose){
    message("Finished gsmi ", gi[1],
            " to ", gi[length(gi)],
            ", time elapsed: ", Sys.time() - tt)
  }
}


# get GSM IDs again




dtables_fromsignal <- function(version, timestamp = NULL, verbose = TRUE, gsmint = 60,
                               fnstem = "mdat.compilation", sepval = " ", getnb = FALSE,
                               idatspath = paste("recount-methylation-files", "idats", sep = "/"),
                               destpath = paste("recount-methylation-analysis",
                                                "files", "mdata", "compilations", sep = "/")){
  # get run metadata
  if(is.null(timestamp)){
    runmd <- get_metadata(version, "notitle")
    nts = runmd[["timestamp"]]
  } else{
    nts <- timestamp
  }
  # get valid gsms idats dir
  if(verbose){message("Getting valid GSM IDs from IDAT filenames...")}
  idats.lf = list.files(idatspath)
  which.valid1 = grepl("\\.idat$", substr(idats.lf,
                                          nchar(idats.lf) - 4,
                                          nchar(idats.lf))) # idat ext
  which.valid2 = grepl(".*hlink.*", idats.lf) # hlink
  which.valid3 = length(gsub("\\..*", "",
                             gsub(".*\\.", "", idats.lf))) > 0 # timestamp
  idats.valid = idats.lf[which.valid1 &
                           which.valid2 &
                           which.valid3]
  gsmu = gsub("\\..*", "", idats.valid)
  gsmu = gsmu[grepl("^GSM.*", gsmu)]
  gsmu = unique(gsmu)
  # red/grn channel files
  if(verbose){message("Finding paired red and grn channel files...")}
  gstr <- gsub(".*hlink\\.", "", gsub("(_Red.idat$|_Grn.idat$)", "", idats.valid))
  rsub <- idats.valid[grepl(".*_Red.idat$", idats.valid)]
  gsub <- idats.valid[grepl(".*_Grn.idat$", idats.valid)]
  # intersect strings between ".*hlink\\." and "(_Red.idat|_Grn.idat)"
  gstr.rg <- intersect(gsub(".*hlink\\.", "", gsub("_Red.idat$", "", rsub)),
                       gsub(".*hlink\\.", "", gsub("_Grn.idat$", "", gsub)))
  idats.validf <- idats.valid[gstr %in% gstr.rg]
  gpath <- unique(gsub("(_Red.idat|_Grn.idat)", "", idats.validf))
  # filter on timestamps
  if(verbose){message("Filtering files on latest timestamps...")}
  gpath.ts <- c()
  gpath.gid <- c()
  gidv <- gsub("\\..*", "", gpath)
  tidv <- gsub(".*\\.", "", gsub("\\.hlink.*", "", gpath))
  # iterate on timestamps
  tsu <- unique(tidv)
  tsu <- tsu[rev(order(tsu))] # sort so that latest comes first
  for(t in tsu){
    gpathf <- gpath[tidv == t & !gidv %in% gpath.gid]
    # rm redundant gsm ids in this iter
    gpathf.gid <- gsub("\\..*", "", gpathf)
    gpathf <- gpathf[!duplicated(gpathf.gid)]
    # append new fn
    gpath.ts <- c(gpath.ts, gpathf)
    # remake gid list
    gpath.gid <- gsub("\\..*", "", gpath.ts)
  }
  gpath <- gpath.ts
  gsmu <- gsub("\\..*", "", gpath)
  if(verbose){message("Finished filtering filenames. ",
                      "Found ", length(gsmu),
                      " valid file pairs.")
  }
  # get GSM ID indices by interval, as list
  gsmii = getblocks(length(gsmu), gsmint)
  # make new data paths
  if(verbose){message("Making and instantiating new data table files...")}
  reds.fn <- paste(paste("redsignal", nts,
                         gsub("\\.", "-", version),
                         sep = "_"),
                   fnstem, sep = ".")
  grns.fn <- paste(paste("greensignal", nts,
                         gsub("\\.", "-", version),
                         sep = "_"),
                   fnstem, sep = ".")
  reds.path = paste(destpath, reds.fn, sep = "/")
  grns.path = paste(destpath, grns.fn, sep = "/")
  if(getnb){
    nb.fn <- paste(paste("noobbeta", nts,
                         gsub("\\.", "-", version),
                         sep = "_"),
                   fnstem, sep = ".")
    nb.path = paste(destpath, nb.fn, sep = "/")
  }
  # instantiate new empty data tables with probes as colnames
  cn = c("gsmi")
  rgi = minfi::read.metharray(c(paste(idatspath, gpath[1:2], sep = "/")))
  rgcni = colnames(t(minfi::getRed(rgi))) # grcni = colnames(t(getGreen(rgi)))
  # rgcni == grcni
  rgcn = matrix(c(cn, rgcni), nrow = 1)
  data.table::fwrite(rgcn, reds.path, sep = sepval, append = FALSE, col.names = F)
  data.table::fwrite(rgcn, grns.path, sep = sepval, append = FALSE, col.names = F)
  if(getnb){
    if(verbose){message("Instantiating noob beta table...")}
    nbi = minfi::preprocessNoob(rgi)
    nbcn = matrix(c(cn, colnames(t(getBeta(nbi)))), nrow = 1)
    data.table::fwrite(nbcn, nb.path, sep = sepval, append = FALSE, col.names = F)
  }
  # append new methdata
  tt = Sys.time()
  for(i in 1:length(gsmii)){
    gi = gsmii[[i]]
    # read in new data
    pathl = paste(idatspath, gpath[gi], sep = "/")
    rgi = minfi::read.metharray(c(pathl))
    
    # get data matrices
    redi = matrix(c(colnames(rgi), t(minfi::getRed(rgi))), ncol = nrow(rgi) + 1)
    grni = matrix(c(colnames(rgi), t(minfi::getGreen(rgi))), ncol = nrow(rgi) + 1)
    
    # append new data
    data.table::fwrite(redi, reds.path, sep = sepval, append = TRUE)
    data.table::fwrite(grni, grns.path, sep = sepval, append = TRUE)
    
    # parse noob-normalized data option
    if(getnb){
      gsi = minfi::preprocessNoob(rgi)
      nbi = matrix(c(colnames(gsi), t(minfi::getBeta(rgi))), ncol = nrow(gsi) + 1)
      data.table::fwrite(nbi, nb.path, sep = sepval, append = TRUE)
    }
    
    if(verbose){
      message("Finished gsmi ", gi[1],
              " to ", gi[length(gi)],
              ", time elapsed: ", Sys.time() - tt)
    }
  }
  return(NULL)
}
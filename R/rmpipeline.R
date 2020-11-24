#!/usr/bin/env R

# Author: Sean Maden
# Functions for `rmpipeline` and `recountmethylation` packages.

#----------
# Utilities
#----------

#' Get list of index blocks
#'
#' Get list of index blocks, allowing for remainders.
#'
#' @param slength Total length of index vector
#' @param bsize Size of index blocks along length
#' @return List of index blocks of min length `slength`/`bsize`
#' @export
getblocks <- function(slength, bsize){
  iv <- list()
  if(slength < bsize){
    iv[[1]] <- seq(1, slength, 1)
  } else{
    sc <- 1; ec <- sc + bsize - 1
    nblocks <- slength %/% bsize
    for(b in 1:nblocks){
      iv[[b]] <- seq(sc, ec, 1)
      sc <- ec + 1; ec <- ec + bsize
    }
    # add final indices
    if(nblocks < (slength/bsize)){
      iv[[length(iv) + 1]] <- seq(sc, slength, 1)
    }
  }
  return(iv)
}

#' Append metadata to HDF5 or SE object
#'
#' Passes version info and timestamp from Python to object metadata
#' @param title Object title
#' @param version Numeric version to be passed, should conform to ##.##.## 
#' nomenclature
#' @param pname Name of pipeline package (default: "rmpipeline")
#' @param sname Name of Python script (default: "get_timestamp.y")
#' @return Metadata content for the object
#' @export
get_metadata <- function(title, version, pname = "rmpipeline",
                         sname = "get_timestamp.py"){
  mdl <- list(title = title, version = version)
  # get timestamp from package python script
  path <- paste(system.file(package = pname), sname, sep="/")
  ts <- system(paste("python", path, sep = " "),
               intern = TRUE, wait = FALSE)
  mdl[["timestamp"]] <- ts
  return(mdl)
}

#----------------------------
# Make data tables from IDATs
#----------------------------

#' Get red and grn signal data tables
#' 
#' Generates 2-channel signal data tables. Validates IDATs, handles new paths 
#' and data.table options, and provides verbose input from main for loop to 
#' write data. 
#'
#' @param version Version of the run for data table filenames.
#' @param timestamp NTP timestamp integer, for filenames (see get_metadata function).
#' @param verbose Whether to return verbose notifications.
#' @param gsmint Number of GSMs to process at a time, typically runs best near 50 samples.
#' @param overwrite Whether to overwrite existing data.table files with same destpath (default TRUE).
#' @param fnstem Filename stem for data tables.
#' @param sepval Separator symbol for data being written (default " ").
#' @param idatspath Path to idat files to read
#' @param getnb Whether to get noob-normalized Beta-vlaues (default: FALSE).
#' @return 
#' @examples
#' #version = "0.0.1"
#' #timestamp = get_metadata("title", version)[["timestamp"]]
#' #dtables_rg(version = version, timestamp = timestamp)
#' @export
dtables_rg <- function(version, timestamp, verbose = TRUE, gsmint = 60,
  overwrite = TRUE, fnstem = "mdat.compilation", sepval = " ",
  idatspath = file.path("recount-methylation-files", "idats"),
  destpath = file.path("recount-methylation-analysis",
    "files", "mdata", "compilations")){
  # get valid gsms idats dir
  idatinfo <- dt_checkidat(idatspath = idatspath, verbose = verbose)
  gpath <- idatinfo[["gpath"]]; gsmu <- idatinfo[["gsmu"]]
  if(verbose){message("Found ", length(gsmu), " GSM IDs with valid idats.")}
  gsmii <- getblocks(length(gsmu), gsmint) # get list of id vectors, e.g. "blocks"
  if(verbose){message("Making new data tables...")}
  dtinfo <- dt_makefiles(gpath = gpath, idatspath = idatspath, 
    destpath = destpath, version = version, nts = timestamp, 
    overwrite = overwrite, sepval = sepval, verbose = verbose)
  if(verbose){message("Wrote data with ", dtinfo[["num.assays"]], " assays.")}
  dtcond <- dtinfo[["dtcond"]]; num.assays = dtinfo[["num.assays"]]
  if(dtcond){rpath <- dtinfo[["reds.path"]]; gpath <- dtinfo[["grns.path"]]
    if(verbose){message("Appending new data for ", length(gsmii)," chunks...")}
    tt <- Sys.time()
    for(i in 1:length(gsmii)){
      dt_write_rg(gi = gsmii[[i]], idatspath = idatspath, gpath = gpath, 
        reds.path = reds.path, grns.path = grns.path, verbose = verbose,
        num.assays = num.assays); te <- Sys.time() - tt
      if(verbose){message("Finished chunk ", i , " time elapsed: ", te)}
    }
  } else{stop("Problem encountered handling data tables.")}
  if(verbose){message("Successfully completed data tables.")}; return(NULL)
}

#' Make and manage data.table files for red and grn signal
#'
#' Handles options for signal tables, returns paths used to write data chunks.
#' @param hlinkv Vector of GSM IDAT hlinks, or the basenames used to read data.
#' @param idatspath Path to idat files to read
#' @param destpath Destination path to new signal data tables.
#' @param version File version information for file names.
#' @param nts NTP timestamp integer, for filenames (see get_metadata function).
#' @param overwrite Whether to overwrite existing data.table files with same 
#' destpath (default TRUE).
#' @param fnstem Filename stem for data tables.
#' @param sepval Separator symbol for tables (" ").
#' @param verbose Whether to return verbose notifications.
#' @return list containing dtcond (try conditions results), and new dt paths 
#' (reds.path and grns.path)
#' @export
dt_makefiles <- function(hlinkv, idatspath, destpath, version, nts, 
  overwrite = TRUE, fnstem = "mdat.compilation", sepval = " ", verbose = TRUE){
  version.fn <- gsub("\\.", "-", version)
  reds.fn <- paste("redsignal", nts, version.fn, sep = "_")
  reds.fn <- paste(reds.fn, fnstem, sep = ".")
  grns.fn <- paste("greensignal", nts, version.fn, sep = "_")
  grns.fn <- paste(grns.fn, fnstem, sep = ".")
  reds.path = file.path(destpath, reds.fn)
  grns.path = file.path(destpath, grns.fn); cn = c("gsmi")
  rgi = minfi::read.metharray(c(file.path(idatspath, hlinkv[1:2])), force = TRUE)
  rgcni = colnames(t(minfi::getRed(rgi))); rgcn = matrix(c(cn, rgcni), nrow = 1)
  if(overwrite){if(verbose){message("Making/verifying data tables...")}
    dt1 <- try(data.table::fwrite(rgcn, reds.path, sep = sepval, 
                                  append = FALSE, col.names = F))
    dt2 <- try(data.table::fwrite(rgcn, grns.path, sep = sepval, 
                                  append = FALSE, col.names = F))
    dtcond <- is.null(dt1) & is.null(dt2)
  } else{dt1 <- try(file.exists(reds.path)); dt2 <- try(file.exists(grns.path))
    dtcond <- dt1 == TRUE & dt1 == TRUE
  }
  lr <- list("dtcond" = dtcond, "reds.path" = reds.path, 
    "grns.path" = grns.path, "num.assays" = nrow(rgi))
  return(lr)
}

#' Validate idats at path, returns gsmu and gpath
#' 
#' Checks for valid idat file pairs, checks for latest timestamps, 
#' files with matching timestamps, and valid red/grn IDAT file pairs, 
#' gets valid gsms idats dir.
#' @param idatspath Path to idat files to read
#' @return list containing gsmu and gpath
#' @export
dt_checkidat <- function(idatspath, verbose = TRUE){
  if(verbose){message("Getting valid GSM IDs from IDAT filenames...")}
  idats.lf = list.files(idatspath)
  which.valid1 = grepl("\\.idat$", substr(idats.lf, nchar(idats.lf) - 4,
                                          nchar(idats.lf))) # idat ext
  which.valid2 = grepl(".*hlink.*", idats.lf) # hlink
  which.valid3 = length(gsub("\\..*", "", gsub(".*\\.", "", idats.lf))) > 0 
  idats.valid = idats.lf[which.valid1 & which.valid2 & which.valid3]
  gsmu = gsub("\\..*", "", idats.valid); gsmu = gsmu[grepl("^GSM.*", gsmu)]
  gsmu = unique(gsmu)
  if(verbose){message("Finding paired red and grn channel files...")}
  gstr <- gsub(".*hlink\\.", "", gsub("(_Red.idat$|_Grn.idat$)", "", idats.valid))
  rsub <- idats.valid[grepl(".*_Red.idat$", idats.valid)]
  gsub <- idats.valid[grepl(".*_Grn.idat$", idats.valid)]
  if(verbose){message("Intersecting .hlink and _Red.idat or _Grn.idat")}
  gstr.rg <- intersect(gsub(".*hlink\\.", "", gsub("_Red.idat$", "", rsub)),
                       gsub(".*hlink\\.", "", gsub("_Grn.idat$", "", gsub)))
  idats.validf <- idats.valid[gstr %in% gstr.rg]
  gpath <- unique(gsub("(_Red.idat|_Grn.idat)", "", idats.validf))
  if(verbose){message("Filtering files on latest timestamps...")}
  gpath.ts <- c()
  gpath.gid <- c()
  gidv <- gsub("\\..*", "", gpath)
  tidv <- gsub(".*\\.", "", gsub("\\.hlink.*", "", gpath))
  tsu <- unique(tidv); tsu <- tsu[rev(order(tsu))] # iterate on timestamps
  for(t in tsu){
    gpathf <- gpath[tidv == t & !gidv %in% gpath.gid]
    if(verbose){message("Removing redundant gsm ids in this iter")}
    gpathf.gid <- gsub("\\..*", "", gpathf)
    gpathf <- gpathf[!duplicated(gpathf.gid)]
    if(verbose){message("Appending new fn")}; gpath.ts <- c(gpath.ts, gpathf)
    if(verbose){message("remaking gid list")}
    gpath.gid <- gsub("\\..*", "", gpath.ts)
  }; gpath <- gpath.ts; gsmu <- gsub("\\..*", "", gpath)
  lr <- list("gsmu" = gsmu, "gpath" = gpath); return(lr)
}

#' Writes chunk of red and grn signal data
#'
#' @param gi Vector of IDAT indices to read from hlinkv.
#' @param hlinkv Vector of GSM IDAT hlink file names, or basenames used to read
#' in data.
#' @param idatspath Path to idat files directory.
#' @param reds.path Path to new red signal data table.
#' @param grns.path Path to new grn signal data table.
#' @param min.cols Minimum data columns to write DNAm data (integer, )
#' @param sepval Separator symbol for data tables.
#' @param num.assays Number of assays to expect before writing (1052641 for 
#' EPIC).
#' @param verbose Whether to include verbose messages.
#' @return NULL, writes data chunks as side effect
#' @export
dt_write_rg <- function(gi, hlinkv, idatspath, gpath, reds.path, grns.path, 
  num.assays = 1052641, sepval = " ", verbose = TRUE){
  if(verbose){message("Reading data...")};pathl=file.path(idatspath, hlinkv[gi])
  upathl <- unique(pathl); rgi = try(minfi::read.metharray(upathl,force = TRUE))
  cond <- nrow(rgi) == num.assays & class(rgi) == "RGChannelSet"
  if(cond){if(verbose){message("getting data matrices")}
    redi = matrix(c(colnames(rgi), 
                    t(minfi::getRed(rgi))), ncol = nrow(rgi) + 1)
    grni = matrix(c(colnames(rgi), 
                    t(minfi::getGreen(rgi))), ncol = nrow(rgi) + 1)
    if(verbose){message("appending new data")}
    data.table::fwrite(redi, reds.path, sep = sepval, append = TRUE)
    data.table::fwrite(grni, grns.path, sep = sepval, append = TRUE)
  } else{if(verbose){message("Error reading data. Skipping...")}}; return(NULL)
}

#------------------------------------
# Make and populate the HDF5 database
#------------------------------------

#' Append sample metadata to HDF5 or SE object
#'
#' Add sample metadata table to the HDF5 database.
#' @param dbn Name of H5 dataset, passed from make_h5db().
#' @param mdpath Path to metadata file.
#' @param dsn Name of new entity in HDF5 database.
#' @param verbose Whether to return verbose status updates.
#' @return Adds metadata table and colnames object to HDF5 database
#' @export
h5_addmd = function(dbn, mdpath, dsn = "mdpost", verbose = TRUE){
  cnn = paste(dsn, "colnames", sep = ".")
  if(verbose){message("Reading in sample metadata from file: ", mdpath)}
  mdt <- get(load(mdpath)); if(verbose){message("Formatting metadata...")}
  mmf = as.matrix(mdt); class(mmf) = "character"; mmf.colnames = colnames(mmf)
  if(verbose){message("Making new entities for HDF5 db...")}
  rhdf5::h5createDataset(dbn, dsn, dims = c(nrow(mmf), ncol(mmf)),
                  maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                  storage.mode = "character",
                  level = 5, chunk = c(10, 16), size = 256)
  rhdf5::h5createDataset(dbn, cnn, dims = length(mmf.colnames),
                  maxdims = c(rhdf5::H5Sunlimited()),storage.mode = "character",
                  level = 5, chunk = c(5), size = 256)
  if(verbose){message("Populating new HDF5 entities...")}
  rhdf5::h5write(mmf, file = dbn, name = dsn,
          index = list(1:nrow(mmf), 1:ncol(mmf)))
  rhdf5::h5write(mmf.colnames, file = dbn, name = cnn,
          index = list(1:length(mmf.colnames)))
  if(verbose){message("Finished adding metadata.")}; rhdf5::h5closeAll()
  return(NULL)
}

#' Append data tables to HDF5 or SE object
#'
#' Add signal data (red and green channel) to the HDF5 database.
#' @param dbn Name of H5 dataset, passed from make_h5db().
#' @param fnl List of signal data table filenames, corresponding (1:1) to names 
#' for new h5 datasets in dsnl.
#' @param fnpath Path to signal data tables.
#' @param dsnl List of data set names in HDF5 database to be populated.
#' @param rmax Total rows to append to data sets, reflecting total samples.
#' @param cmax Total columns to append to data sets, reflecting total assays or 
#' probes.
#' @param verbose Whether to print verbose progress messages.
#' @param nr.inc Number of samples to append at a time (default: 10).
#' @return Populates the HDF5 database
#' @export
h5_addtables = function(dbn, fnl, fnpath, dsnl, rmax, cmax,
                    verbose = TRUE, nr.inc = 10){
  for(di in 1:length(dsnl)){
    fnread = fnl[di]; dsn = dsnl[di]; tt <- Sys.time()
    rhdf5::h5createDataset(dbn, dsn, dims = c(rmax, cmax),
                           maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                           storage.mode = "double", level = 5, chunk = c(1, 5))
    rn = cn = c(); con <- file(paste(fnpath, fnread, sep = "/"), "r")
    cn = unlist(strsplit(readLines(con, n = 1), " ")); cn = cn[2:length(cn)]
    cn = gsub("\n", "",gsub('\"', '', cn[1:cmax]));nri = getblocks(rmax, nr.inc)
    for(ni in nri){wgsm = c()
      i = ni[1]; dati = unlist(strsplit(readLines(con, n = length(ni)), " "))
      wdi = which(grepl(".*GSM.*", dati)); dff = matrix(nrow = 0, ncol = cmax)
      ngsm = gsub("\n", "", gsub('\"', '', gsub("\\..*", "", dati[wdi])))
      for(wi in 1:length(wdi)){
        if(!ngsm[wi] %in% rn){wadd = wdi[wi] + 1
          dff = rbind(dff, matrix(dati[wadd:(wadd + cmax - 1)], nrow = 1))
          wgsm = c(wgsm, ngsm[wi])
        }
      }; rn = c(rn, wgsm); class(dff) = "numeric"
      rhdf5::h5write(dff, file = dbn, name = dsn,
                     index = list(ni[1]:ni[length(ni)], 1:cmax))
      rhdf5::h5closeAll()
      message("For ds ", dsn,", finished reading index ", i,
              " to ", ni[length(ni)], ", time; ", Sys.time() - tt)
    }
    message("Adding row and column names for ds ", dsn)
    cnn = paste0(dsn, ".colnames"); rnn = paste0(dsn, ".rownames");
    rhdf5::h5createDataset(dbn, cnn, dims = length(cn), maxdims = c(rhdf5::H5Sunlimited()),
                           storage.mode = "character", level = 5, # compression level, 1-9
                           chunk = c(20), size = 256) # chunk dims
    message("Added colnames.")
    rhdf5::h5createDataset(dbn, rnn, dims = length(rn), maxdims = c(rhdf5::H5Sunlimited()),
                           storage.mode = "character", level = 5, # compression level, 1-9
                           chunk = c(20), size = 256) # chunk dims
    message("Added rownames.")
    rhdf5::h5write(cn, file = dbn, name = cnn, index = list(1:length(cn)))
    rhdf5::h5write(rn, file = dbn, name = rnn, index = list(1:length(rn)))
    rhdf5::h5closeAll(); message("Completed writing data for ds, ", dsn)
  }
  if(verbose){message("Finished adding red and green channel data sets.")}
  rhdf5::h5closeAll(); return(NULL)
}

#' Make and populate a new HDF5 database with red and green signal, and metadata
#'
#' Add signal data (red and green channel) to the HDF5 database.
#' @param dbfnstem Stem of filename for h5 database (to which ts and version 
#' appended)
#' @param version Version of new database file.
#' @param ts Timestamp of new database file.
#' @param fnl Vector of signal tables containing data to be added.
#' @param fnpath Path to dir containing files in fnl.
#' @param dsnl Vector of data set names in HDF5 database to be populated.
#' @param rmax Total rows to append to data sets, reflecting total samples.
#' @param cmax Total columns to append to data sets, reflecting total assays or 
#' probes. Should be 622399 for raw red/grn signal.
#' @param newtables Whether to also add new data tables (noob-norm. Beta-values, 
#' meth. and unmeth. signal).
#' @param mdpath If addmd, the path to the metadata file to load.
#' @param nr.inc Number of samples to append at a time (default: 10, passed to 
#' `h5_addtables()`).
#' @param ngsm.block Number of GSMs (samples) per process block (default 50, 
#' passed to `h5_newtables()`).
#' @param verbose Whether to show verbose status updates.
#' @param dsnl Vector of new h5 data set objects in the h5db object, 
#' corresponding (1:1) to the files declared in fnl.
#' @return Populates the HDF5 database
#' @export
makeh5db_rg <- function(dbfnstem, version, ts, fnl, fnpath,
                      rmax = 35300, cmax = 622399,
                      rmoldh5 = TRUE, newtables = FALSE, mdpath = NULL,
                      nr.inc = 10, ngsm.block = 50, verbose = TRUE,
                      dsnl = c("redsignal", "greensignal", "mdpost")){
  require(rhdf5);fn <- paste(dbfnstem, ts, gsub("\\.", "-", version), sep = "_")
  dbn <- paste(fn, "h5", sep = ".")
  if(verbose){message("Making new h5 db file: ", dbn)}
  suppressMessages(try(rhdf5::h5createFile(dbn), silent = TRUE))
  if(rmoldh5){if(verbose){message("Removing old data...")}
    for(d in dsnl){suppressMessages(try(rhdf5::h5delete(dbn, d), silent = TRUE))
      suppressMessages(try(rhdf5::h5delete(dbn, paste0(d, ".colnames")), 
        silent = TRUE))
      suppressMessages(try(rhdf5::h5delete(dbn, paste0(d, ".rownames")), 
        silent = TRUE))
    }
  }
  # add red and grn signal tables
  if(verbose){message("Adding and populating data tables to HDF5 database")}
  h5_addtables(dbn = dbn, fnl = fnl, dsnl = c("redsignal", "greensignal"),
               rmax = rmax, cmax = cmax, nr.inc = ngsm.block,
               fnpath = fnpath, verbose = verbose)
  if(verbose){message("Finished adding red and green channel data.")}
  if(!is.null(mdpath)){if(verbose){message("Adding sample metadata to db...")}
    h5_addmd(dbn, mdpath, verbose = verbose)
  } else{if(verbose){message("No mdpath specified!")}}
  rhdf5::h5closeAll();if(verbose){message("Finished all processes. Returning.")}
  return(NULL)
}





#-----------------------
# Make the SE-H5 objects
#-----------------------

#' Retrieve samples metadata from HDF5 db.
#'
#' Retrieves sample metadata from an HDF5 database.
#'
#' @param dbn Path to HDF5 database file.
#' @param dsn Name or group path to HDF5 dataset containing postprocessed 
#' metadata.
#' @return Postprocessed metadata as a `data.frame`.
#' @export
data_mdpost = function(dbn, dsn){
  mdp = as.data.frame(rhdf5::h5read(file = dbn, name = dsn), stringsAsFactors = F)
  colnames(mdp) = rhdf5::h5read(file = dbn, name = paste(dsn, "colnames", sep = "."))
  return(mdp)
}

#' Append phenotype data to a SummarizedExperiment object
#'
#' Append phenotype data to a SummarizedExperiment object.
#' Note, this data should be a data.frame with certain data specified by 
#' colnames "basename" and "gsm".
#'
#' @param version The data file version
#' @param ts Data file NTP timestamp (integer)
#' @param dbn Name of h5 file to load data matrices from.
#' @param fnstem File name stem for new H5SE file
#' @param dsnv Vector of the data set file names to read from h5 file.
#' @param dsn.md Name of metadata file in h5 file.
#' @param verbose Whether to show verbose status messages (default TRUE).
#' @param addmd Whether to add sample metadata.
#' @param mdpath External path to sample metadata. If FALSE, looks for dsn.md 
#' in the dbn h5 file.
#' @param dsn.rnv Vector for rownames datasets for dsnv sets.
#' @param dsn.cnv Vector of colnames datasets for dsnv sets.
#' @param semd H5SE object metadata information.
#' @returns New H5SE object.
#' @export
se_addpheno <- function(phenopath, se, pdat = NULL, verbose = TRUE){
  if(is.null(pdat)){
    if(verbose){message("Loading pheno data from path...")}
    mdp <- get(load(phenopath))
  } else{
    mdp <- pdat
  }
  # Adds pheno data to a SummarizedExperiment objects
  mdp <- mdp[mdp$basename %in% colnames(se),]
  bnv <- colnames(se)
  gsmv <- gsub("\\..*", "", bnv)
  gsmfilt <- !gsmv %in% mdp$gsm
  gsm.ov <- gsmv[gsmfilt]
  bnv.ov <- bnv[gsmfilt]
  # add na values
  nacol <- rep("NA", length(gsm.ov))
  mdp.ov <- matrix(c(gsm.ov,
                     rep(nacol, 6),
                     bnv.ov,
                     rep(nacol, 11)),
                   nrow = length(gsm.ov))
  colnames(mdp.ov) <- colnames(mdp)
  mdp <- rbind(mdp, mdp.ov)
  mdp <- mdp[order(match(mdp$gsm, gsmv)),]
  if(identical(mdp$gsm, gsmv)){
    # add valid basenames
    mdp$basename <- colnames(se)
    rownames(mdp) <- colnames(se)
    identical(rownames(mdp), colnames(se))
    # append pheno to se object
    if(verbose){message("Adding metadata to se.")}
    minfi::pData(se) <- S4Vectors::DataFrame(mdp)
    return(se)
  } else{
    stop("There was a problem matching GSM IDs between mdp, se.")
  }
  return(NULL)
}

#' Blah
#'
#' Blah
#' @param version The data file version.
#' @param ts Data file NTP timestamp (integer).
#' @returns New h5se object.
#' @export
make_h5se_rg <- function(version = "0.0.1", ts = 1590090412, dbn = "remethdb_1590090412_0-0-1.h5",
                         newfnstem = "remethdb_h5se-rg", dsnv = c("redsignal", "greensignal"), 
                         dsn.md = "mdpost", dsn.md.cn = paste0(dsn.md, ".colnames"), 
                         verbose = TRUE, addmd = TRUE, mdpath = NULL, replace.opt = TRUE,
                         dsn.rnv = c(paste0(dsnv[1], ".rownames"), paste0(dsnv[2], ".rownames")),
                         dsn.cnv = c(paste0(dsnv[1], ".colnames"), paste0(dsnv[2], ".colnames")),
                         semd = list("title" = "Recount methylation H5-SE object", "version" = version, 
                                     "timestamp" = ts, "preprocessing" = "raw")){
  newfn <- paste(newfnstem, gsub("\\.", "-", version), ts, sep = "_")
  if(verbose){message("Setting annotation info...")}
  anno = c("IlluminaHumanMethylation450k", "ilmn12.hg19")
  names(anno) = c("array", "annotation")
  if(verbose){message("Getting h5 datasets...")}
  ldat <- list()
  for(i in 1:length(dsnv)){
    nb <- HDF5Array::HDF5Array(dbn, dsnv[i])
    rn <- rhdf5::h5read(dbn, dsn.rnv[i])
    cn <- rhdf5::h5read(dbn, dsn.cnv[i])
    nb <- nb[c(1:length(rn)),c(1:length(cn))]
    rownames(nb) <- as.character(rn)
    colnames(nb) <- as.character(cn)
    nb <- t(nb)
    ldat[[dsnv[i]]] <- nb
  }
  # Make the new H5-SE set(s)
  if(verbose){message("Making the H5SE set...")}
  gri <- minfi::RGChannelSet(Red = ldat[[1]], Green = ldat[[2]], anno = anno)
  S4Vectors::metadata(gri) <- semd
  # Parse md file options
  if(addmd){
    if(verbose){message("Adding samples metadata...")}
    if(is.null(mdpath)){
      if(dsn.md %in% h5ls(dbn)$name & dsn.md.cn %in% h5ls(dbn)$name){
        if(verbose){message("Adding metadata from dbn...")}
        mdp <- recountmethylation::data_mdpost(dbn, dsn.md)
        gri <- se_addpheno(pdat = mdp, se = gri)
      }
    } else{
      if(verbose){message("Adding metadata from phenopath...")}
      gri <- se_addpheno(mdpath, se = gri)
    }
  }
  # start to instantiate the H5SE object
  t1 <- Sys.time()
  if(verbose){message("Adding data to new file ", newfn, "...")}
  HDF5Array::saveHDF5SummarizedExperiment(gri, dir = newfn, replace = replace.opt)
  if(verbose){message("Data additions complete, time elapsed:", Sys.time() - t1)}
  return(NULL)
}

#' Make an HDF5 h5 file with MethylSet data, from an RGChannelSet 
#' h5se file
#'
#' Makes an HDF5 h5 database file with raw/unnormalized methylated and 
#' unmethylated signal, including rows and columns as datasets.
#' @param dbn Path to the RGChannelSet h5se object.
#' @param version Version for new h5 file.
#' @param ts Timestamp for new h5 file.
#' @param newfnstem New filename stem for h5 file.
#' @param verbose Whether to return verbose messages (default TRUE).
#' @param replace.opt Whether to replace h5 files with duplicate filename 
#' (default TRUE).
#' @param blocksize Size of blocks/number of samples to process at a time 
#' (default 65).
#' @returns Path to new h5 object.
#' @export
make_h5_gm <- function(dbn, version, ts, newfnstem = "remethdb_h5", 
                        verbose = TRUE, replace.opt = TRUE, blocksize = 65){
  rg <- HDF5Array::loadHDF5SummarizedExperiment(dbn)
  message("Making new h5 file and datasets.")
  h5dbn <- paste(newfnstem, se, version, ts, sep = "_")
  h5dbn <- paste0(c(h5dbn, "h5"), collapse = ".") 
  rhdf5::h5createFile(h5dbn)
  dsv <- c("meth", "unmeth")
  for(d in dsv){
    rhdf5::h5createDataset(h5dbn, d, dims = list(485512, 36240), 
      maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()), 
      storage.mode = "double", level = 5, chunk = c(1000, 50))
  }
  rhdf5::h5createDataset(h5dbn, "rownames", dims = list(485512), 
      maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character", size = 256, 
      level = 5, chunk = c(1000))
  rhdf5::h5createDataset(h5dbn, "colnames", dims = list(36240), 
    maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character", size = 256, 
    level = 5, chunk = c(1000))
  if(verbose){message("Getting blocks of sample indices.")}
  blocks <- getblocks(ncol(rg), bsize = blocksize)
  if(verbose){message("Writing first block of data.")}
  b <- 1
  cindices <- blocks[[b]]
  rgi <- rg[, cindices]
  gm <- minfi::preprocessRaw(rgi)
  meth <- as.matrix(minfi::getMeth(gm))
  unmeth <- as.matrix(minfi::getUnmeth(gm))
  rhdf5::h5write(meth, file = h5dbn, name = "meth", 
          index = list(1:485512, cindices[1]:cindices[length(cindices)]))
  rhdf5::h5write(unmeth, file = h5dbn, name = "unmeth", 
          index = list(1:485512, cindices[1]:cindices[length(cindices)]))
  if(verbose){message("Writing rows and cols.")}
  rhdf5::h5write(colnames(rg), file = h5dbn, name = "colnames", index = list(1:ncol(rg)))
  rhdf5::h5write(rownames(gm), file = h5dbn, name = "rownames", index = list(1:485512))
  message("Writing data blocks.")
  t1 <- Sys.time()
  for(b in 2:length(blocks)){
    cindices <- blocks[[b]]
    gm <- minfi::preprocessRaw(rg[,cindices])
    meth <- as.matrix(minfi::getMeth(gm))
    unmeth <- as.matrix(minfi::getUnmeth(gm))
    write.indices <- cindices[1]:cindices[length(cindices)]
    rhdf5::h5write(meth, index = list(1:485512, write.indices), file = h5dbn, 
        name = "meth")
    rhdf5::h5write(unmeth, index = list(1:485512, write.indices), file = h5dbn, 
        name = "unmeth")
    if(verbose){
      message("Finished block ", b, ", time elapsed = ", Sys.time() - t1)
    }
  }
  message("Finished writing data blocks. Returning h5 name.")
  return(h5dbn)
}


#' Make an HDF5 h5 file with GenomicRatioSet data, from an RGChannelSet 
#' h5se file
#'
#' Makes an HDF5 h5 file containing noob-normalized Beta-values, including
#' rows and columns as separate datasets.
#' @param dbn Path to the RGChannelSet h5se object.
#' @param version Version for new h5 file.
#' @param ts Timestamp for new h5 file.
#' @param newfnstem New filename stem for h5 file.
#' @param replace.opt Whether to replace h5 files with duplicate filename 
#' (default TRUE).
#' @param blocksize Size of blocks/number of samples to process at a time 
#' (default 65).
#' @param verbose Whether to return verbose messages (default TRUE).
#' @returns Path to new h5 object.
#' @export
make_h5_gr <- function(dbn, version, ts, newfnstem = "remethdb_h5",
                        replace.opt = TRUE, blocksize = 65, verbose = TRUE){
  rg <- HDF5Array::loadHDF5SummarizedExperiment(dbn)
  message("Making new h5 file and datasets.")
  h5dbn <- paste(newfnstem, se, version, ts, sep = "_")
  h5dbn <- paste0(c(h5dbn, "h5"), collapse = ".") 
  rhdf5::h5createFile(h5dbn)
  rhdf5::h5createDataset(h5dbn, "noobbeta", dims = list(485512, 36240), 
      maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()), 
      storage.mode = "double", level = 5, chunk = c(1000, 50))
  rhdf5::h5createDataset(h5dbn, "rownames", dims = list(485512), 
      maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character", size = 256, 
      level = 5, chunk = c(1000))
  rhdf5::h5createDataset(h5dbn, "colnames", dims = list(36240), 
    maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character", size = 256, 
    level = 5, chunk = c(1000))
  if(verbose){message("Getting blocks of sample indices.")}
  blocks <- getblocks(ncol(rg), bsize = blocksize)
  if(verbose){message("Writing first block of data.")}
  b <- 1
  cindices <- blocks[[b]]
  rgi <- rg[, cindices]
  gr <- minfi::preprocessNoob(rgi)
  nbeta <- as.matrix(minfi::getBeta(gr))
  rhdf5::h5write(nbeta, file = h5dbn, name = "noobbeta", 
          index = list(1:485512, cindices[1]:cindices[length(cindices)]))
  if(verbose){message("Writing rows and cols.")}
  rhdf5::h5write(colnames(rg), file = h5dbn, name = "colnames", index = list(1:ncol(rg)))
  rhdf5::h5write(rownames(gr), file = h5dbn, name = "rownames", index = list(1:485512))
  message("Writing data blocks.")
  t1 <- Sys.time()
  for(b in 2:length(blocks)){
    cindices <- blocks[[b]]
    gr <- minfi::preprocessNoob(rg[,cindices])
    nbeta <- as.matrix(minfi::getBeta(gr))
    write.indices <- cindices[1]:cindices[length(cindices)]
    rhdf5::h5write(nbeta, file = h5dbn, name = "noobbeta", 
          index = list(1:485512, cindices[1]:cindices[length(cindices)]))
    if(verbose){
      message("Finished block ", b, ", time elapsed = ", Sys.time() - t1)
    }
  }
  message("Finished writing data blocks. Returning h5 name.")
  return(h5dbn)
}

#' Make an h5se object from an h5se gm file
#' 
#' @param version Version for new filenames.
#' @param ts NTP timestamp.
#' @param h5name Name of the h5 file to read data from.
#' @param pdata Metadata object obtained from pData(se). 
#' @param replaceopt Whether to overwrite existing file of same name as new file
#' (default = TRUE).
#' @param newdnstem Stem for new h5se object file.
#' @return NULL, generates a new h5se file.
#' @export
make_h5se_gm <- function(version, ts, h5name, pdata, replaceopt = TRUE,
    newdnstem = "remethdb-h5se"){
  require(minfiData)
  message("Getting new name.")
  newfn <- paste0(newdnstem, "gm", version, ts, sep = "_")
  message("Getting granges.")
  ms <- minfi::mapToGenome(get(data("MsetEx")))
  anno <- minfi::annotation(ms)
  gr <- GenomicRanges::granges(ms)
  message("Reading datasets.")
  meth <- HDF5Array::HDF5Array(h5name, "meth")
  unmeth <- HDF5Array::HDF5Array(h5name, "unmeth")
  cn <- rhdf5::h5read(h5name, "colnames")
  rn <- rhdf5::h5read(h5name, "rownames")
  rownames(meth) <- rownames(unmeth) <- as.character(rn)
  colnames(meth) <- colnames(unmeth) <- as.character(cn)
  gr <- gr[order(match(names(gr), rownames(meth)))]
  icond <- identical(rownames(meth), names(gr))
  if(!icond){stop("Problem matching gr names to h5 meth rows")}
  icond <- identical(rownames(unmeth), names(gr))
  if(!icond){stop("Problem matching gr names to h5 unmeth rows")}
  icond <- identical(colnames(meth), rownames(pdata))
  if(!icond){stop("Problem matching pdata rows to h5 meth cols")}
  icond <- identical(colnames(unmeth), rownames(pdata))
  if(!icond){stop("Problem matching pdata rows to h5 unmeth cols")}
  message("Writing new data to h5se object. This may take awhile.")
  meth <- HDF5Array::HDF5Array(h5name, "meth")
  unmeth <- HDF5Array::HDF5Array(h5name, "unmeth")
  cn <- rhdf5::h5read(h5name, "colnames")
  rn <- rhdf5::h5read(h5name, "rownames")
  rownames(meth) <- rownames(unmeth) <- as.character(rn)
  colnames(meth) <- colnames(unmeth) <- as.character(cn)
  gm <- minfi::GenomicMethylSet(gr = gr, Meth = meth, Unmeth = unmeth, 
                          annotation = anno)
  pData(gm) <- pdata
  HDF5Array::saveHDF5SummarizedExperiment(gm, dir = newfn, replace = replaceopt)
  message("Finished writing h5se data.")
  return(NULL)
}

#' Make an h5se object from an h5se gr file
#' 
#' @param version Version for new filenames.
#' @param ts NTP timestamp.
#' @param h5name Name of the h5 file to read data from.
#' @param pdata Metadata object obtained from pData(se). 
#' @param replaceopt Whether to overwrite existing file of same name as new file
#' (default = TRUE).
#' @param newdnstem Stem for new h5se object file.
#' @return NULL, generates a new h5se file.
#' @export
make_h5se_gr <- function(version, ts, h5name, pdata, replaceopt = TRUE,
    newdnstem = "remethdb-h5se"){
  require(minfiData)
  message("Getting new name.")
  newfn <- paste0(newdnstem, "gr", version, ts, sep = "_")
  message("Getting granges and annotation.")
  ms <- minfi::mapToGenome(get(data("MsetEx")))
  anno <- minfi::annotation(ms)
  gr <- GenomicRanges::granges(ms)
  message("Reading datasets.")
  nbeta <- HDF5Array::HDF5Array(h5name, "noobbeta")
  cn <- rhdf5::h5read(h5name, "colnames")
  rn <- rhdf5::h5read(h5name, "rownames")
  rownames(nbeta) <- as.character(rn)
  colnames(nbeta) <- as.character(cn)
  gr <- gr[order(match(names(gr), rownames(nbeta)))]
  icond <- identical(rownames(nbeta), names(gr))
  if(!icond){stop("Problem matching gr names to h5 nbeta rows")}
  icond <- identical(colnames(nbeta), rownames(pdata))
  if(!icond){stop("Problem matching pdata rows to h5 nbeta cols")}
  message("Writing new data to h5se object. This may take awhile.")
  nbeta <- HDF5Array::HDF5Array(h5name, "noobbeta")
  cn <- rhdf5::h5read(h5name, "colnames")
  rn <- rhdf5::h5read(h5name, "rownames")
  rownames(nbeta) <- as.character(rn)
  colnames(nbeta) <- as.character(cn)
  gmi <- minfi::GenomicRatioSet(gr = gr, Beta = nbeta, annotation = anno)
  minfi::pData(gmi) <- pdata
  HDF5Array::saveHDF5SummarizedExperiment(gmi, dir = newfn, replace = replaceopt)
  message("Finished writing h5se data.")
  return(NULL)
}
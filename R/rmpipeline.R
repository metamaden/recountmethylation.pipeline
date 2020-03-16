#!/usr/bin/env R

# Functions for the Recount Methylation Pipeline and `recountmethylation` R package.

#----------
# Utilities
#----------

#' Get list of index blocks
#'
#' Get list of index blocks allowing for remainders.
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
#' @param version Numeric version to be passed, should conform to ##.##.## nomenclature
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

#' Make data tables from IDATs
#'
#' Make data tables from IDATs, including red and green channel signals and noob-normalized Beta-values.
#' This function should be run from the "recount-methylation" base directory.
#'
#' @param version Version of the run for data table filenames.
#' @param verbose Whether to return verbose notifications.
#' @param gsmint Number of GSMs to process at a time, typically runs best near 50 samples.
#' @param fnstem Filename stem for data tables.
#' @param sepval Separator symbol for data being written.
#' @param getnb Whether to get noob-normalized Beta-vlaues (default: FALSE).
#' @return Creates new HDF5 database from DNAm signal data.
#' @export

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
  # name of new colnames entity
  cnn = paste(dsn, "colnames", sep = ".")
  # read in data
  if(verbose){message("Reading in sample metadata from file: ", mdpath)}
  mdt <- get(load(mdpath))
  # format data
  if(verbose){message("Formatting sample metadata...")}
  mmf = as.matrix(mdt)
  class(mmf) = "character"
  mmf.colnames = colnames(mmf)
  # instantiate new entities in db
  if(verbose){message("Making new entities for HDF5 db...")}
  rhdf5::h5createDataset(dbn, dsn, dims = c(nrow(mmf), ncol(mmf)),
                  maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                  storage.mode = "character",
                  level = 5, chunk = c(10, 16), size = 256)
  rhdf5::h5createDataset(dbn, cnn, dims = length(mmf.colnames),
                  maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character",
                  level = 5, chunk = c(5), size = 256)
  # write new data
  if(verbose){message("Populating new HDF5 entities...")}
  rhdf5::h5write(mmf, file = dbn, name = dsn,
          index = list(1:nrow(mmf), 1:ncol(mmf)))
  rhdf5::h5write(mmf.colnames, file = dbn, name = cnn,
          index = list(1:length(mmf.colnames)))
  if(verbose){message("Finished adding metadata to HDF5 db.")}
  rhdf5::h5closeAll() # close open connections to db
  return(NULL)
}

#' Append data tables to HDF5 or SE object
#'
#' Add signal data (red and green channel) to the HDF5 database.
#' @param dbn Name of H5 dataset, passed from make_h5db().
#' @param fnl List of signal data table filenames, corresponding (1:1) to names for new h5 datasets in dsnl.
#' @param fnpath Path to signal data tables.
#' @param dsnl List of data set names in HDF5 database to be populated.
#' @param rmax Total rows to append to data sets, reflecting total samples.
#' @param cmax Total columns to append to data sets, reflecting total assays or probes.
#' @param verbose Whether to print verbose progress messages.
#' @param nr.inc Number of samples to append at a time (default: 10).
#' @return Populates the HDF5 database
#' @export
h5_addtables = function(dbn, fnl, fnpath, dsnl, rmax, cmax,
                    verbose = TRUE, nr.inc = 10){
  for(di in 1:length(dsnl)){
    fnread = fnl[di]; dsn = dsnl[di]
    rhdf5::h5createDataset(dbn, dsn, dims = c(rmax, cmax),
                           maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                           storage.mode = "double", level = 5, chunk = c(1, 5))
    rn = cn = c()
    con <- file(paste(fnpath, fnread, sep = "/"), "r") # make new connection object
    cn = unlist(strsplit(readLines(con, n = 1), " ")) # read first line (colnames)
    cn = cn[2:length(cn)] # filt first value
    cn = gsub("\n", "",gsub('\"', '', cn[1:cmax])) # grab the max indexed value
    nri = getblocks(rmax, nr.inc)
    tt <- Sys.time()
    for(ni in nri){
      i = ni[1]
      # read new lines; note this is contextual/resumes at next new line each iter
      dati = unlist(strsplit(readLines(con, n = length(ni)), " "))
      wdi = which(grepl(".*GSM.*", dati))
      dff = matrix(nrow = 0, ncol = cmax)
      # get stripped GSM IDs from rownames (first col each line)
      ngsm = gsub("\n", "",
                  gsub('\"', '',
                       gsub("\\..*", "", dati[wdi])))
      wgsm = c() # new gsm ids to write
      for(wi in 1:length(wdi)){
        # filter redundant gsm ids
        if(!ngsm[wi] %in% rn){
          wadd = wdi[wi] + 1
          dff = rbind(dff,
                      matrix(dati[wadd:(wadd + cmax - 1)],
                             nrow = 1))
          wgsm = c(wgsm, ngsm[wi])
        }
      }
      rn = c(rn, wgsm) # gsm ids for new data
      class(dff) = "numeric"
      rhdf5::h5write(dff, file = dbn, name = dsn,
                     index = list(ni[1]:ni[length(ni)], 1:cmax))
      rhdf5::h5closeAll()
      message("For ds ", dsn,", finished reading index ", i,
              " to ", ni[length(ni)],
              ", time; ", Sys.time() - tt)
    }
    message("Adding row and column names for ds ", dsn)
    cnn = paste0(dsn, ".colnames"); rnn = paste0(dsn, ".rownames");
    rhdf5::h5createDataset(dbn, cnn, dims = length(cn), maxdims = c(rhdf5::H5Sunlimited()),
                           storage.mode = "character", level = 5, # compression level, 1-9
                           chunk = c(20), size = 256) # chunk dims
    message("Added colnames...")
    rhdf5::h5createDataset(dbn, rnn, dims = length(rn), maxdims = c(rhdf5::H5Sunlimited()),
                           storage.mode = "character", level = 5, # compression level, 1-9
                           chunk = c(20), size = 256) # chunk dims
    message("Added rownames...")
    rhdf5::h5write(cn, file = dbn, name = cnn, index = list(1:length(cn)))
    rhdf5::h5write(rn, file = dbn, name = rnn, index = list(1:length(rn)))
    rhdf5::h5closeAll()
    message("Completed writing data, rownames, colnames for ds, ", dsn)
  }
  if(verbose){message("Finished adding red and green channel data sets.")}
  rhdf5::h5closeAll() # close open connections to db
  return(NULL)
}

#' Make new minfi data tables from base h5 db
#'
#' Uses minfi to populate new data tables (noobbeta, methylated and unmethylated signal) from h5 db with red and green signal data.
#'
#' @param dbn Name or path to h5 db.
#' @param dsn.nb Name of the new noobbeta h5 dataset.
#' @param dsn.meth Name of the new methylated h5 dataset.
#' @param dsn.unmeth Name of the new unmethylated h5 dataset
#' @param dsn.red Name of existing red signal h5 dataset.
#' @param dsn.grn Name of existing green signal h5 dataset.
#' @param verbose Whether to show verbose status messages.
#' @param ngsm.block Number of GSMs (samples) per process block (default 50).
#' @param ncol.chunk Number of columns (GSMs/samples) per chunk in saved h5 datasets.
#' @return Adds new minfi h5 datasets to specified h5 db (dbn).
#' @export
h5_newtables <- function(dbn, dsn.nb = "noobbeta",
                         dsn.meth = "methylated_signal",
                         dsn.unmeth = "unmethylated_signal",
                         dsn.red = "redsignal", dsn.grn = "greensignal",
                         verbose = TRUE, ngsm.block = 50,
                         ncol.chunk = 5000){
  # Generate noobbeta and meth/unmeth signal tables
  # get dimensions from red and grn signal data
  if(verbose){message("Getting red and green signal h5 object info...")}
  rs.rn <- rhdf5::h5read(dbn, paste0(dsn.red, ".rownames"))
  rs.cn <- rhdf5::h5read(dbn, paste0(dsn.red, ".colnames"))
  gs.rn <- rhdf5::h5read(dbn, paste0(dsn.grn, ".rownames"))
  gs.cn <- rhdf5::h5read(dbn, paste0(dsn.grn, ".colnames"))
  if(verbose){message("Getting blocked indices of row data to process...")}
  sbv <- getblocks(length(rs.rn), ngsm.block)
  # get new cg  dimensions
  if(verbose){message("Getting manifest and setting new colname dim...")}
  anno.name = "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  man = eval(parse(text = paste(anno.name, "Manifest", sep = "::")))
  ncg = nrow(man)
  # new h5 data params
  newdims <- c(length(rs.rn), ncg)
  chunkvars <- c(5, ncol.chunk)
  # make new tables
  rhdf5::h5createDataset(dbn, "unmethylated_signal", dims = newdims,
                         maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                         storage.mode = "double", level = 5, chunk = chunkvars)
  rhdf5::h5createDataset(dbn, "methylated_signal", dims = newdims,
                         maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                         storage.mode = "double", level = 5, chunk = chunkvars)
  rhdf5::h5createDataset(dbn, "noobbeta", dims = newdims,
                         maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                         storage.mode = "double", level = 5, chunk = chunkvars)
  # append new data
  if(verbose){message("Appending new data to h5 datasets...")}
  t1 <- Sys.time()
  for(i in 1:length(sbv)){
    b = sbv[[i]]
    gsmvi <- gsub("\\..*", "", rs.rn[b])
    se.rgi = recountmethylation::getrg(gsmv = gsmvi,
                                       cgv = "all", dbn = dbn,
                                       data.type = "se",
                                       metadata = FALSE,
                                       verbose = FALSE)
    # new se and data objects
    se.nb <- minfi::preprocessNoob(se.rgi)
    methb <- t(minfi::getMeth(se.nb))
    unmethb <- t(minfi::getUnmeth(se.nb))
    nb <- t(minfi::getBeta(se.nb))
    # append new data to h5 data
    writei <- list(b[1]:b[length(b)], 1:ncg)
    rhdf5::h5write(unmethb, file = dbn,
                   name = dsn.unmeth, index = writei)
    rhdf5::h5write(methb, file = dbn,
                   name = dsn.meth, index = writei)
    rhdf5::h5write(nb, file = dbn,
                   name = dsn.nb, index = writei)
    if(verbose){
      message("finished block ", i," of ",
              length(sbv),", time elapsed: ",
              Sys.time() - t1)
    }
  }
  if(verbose){message("Completed addition of new tables!")}
  rhdf5::h5closeAll() # close open connections to db
  return(NULL)
}

#' Make and populate a new HDF5 database
#'
#' Add signal data (red and green channel) to the HDF5 database.
#' @param dbfnstem Stem of filename for h5 database (to which ts and version appended)
#' @param version Version of new database file.
#' @param ts Timestamp of new database file.
#' @param fnl Vector of signal tables containing data to be added.
#' @param fnpath Path to dir containing files in fnl.
#' @param dsnl Vector of data set names in HDF5 database to be populated.
#' @param rmax Total rows to append to data sets, reflecting total samples.
#' @param cmax Total columns to append to data sets, reflecting total assays or probes. Should be 622399 for raw red/grn signal.
#' @param newtables Whether to also add new data tables (noob-norm. Beta-values, meth. and unmeth. signal).
#' @param addmd Whether to add metadata to the h5db object.
#' @param mdpath If addmd, the path to the metadata file to load.
#' @param nr.inc Number of samples to append at a time (default: 10, passed to `h5_addtables()`).
#' @param ngsm.block Number of GSMs (samples) per process block (default 50, passed to `h5_newtables()`).
#' @param verbose Whether to show verbose status updates.
#' @param dsnl Vector of new h5 data set objects in the h5db object, corresponding (1:1) to the files declared in fnl.
#' @return Populates the HDF5 database
#' @export
make_h5db <- function(dbfnstem, version, ts, fnl, fnpath,
                      rmax = 35300, cmax = 622399,
                      rmoldh5 = TRUE, newtables = TRUE,
                      addmd = TRUE, mdpath = NULL,
                      nr.inc = 10, ngsm.block = 50, verbose = TRUE,
                      dsnl = c("redsignal", "greensignal", "mdpost",
                               "unmethylated_signal", "methylated_signal",
                               "noobbeta")){

  # make new h5 db filename
  dbn <- paste(paste(dbfnstem, ts,
                     gsub("\\.", "-", version),
                     sep = "_"), "h5",
               sep = ".")
  if(verbose){message("Making new h5 db file: ", dbn)}
  suppressMessages(try(rhdf5::h5createFile(dbn), silent = TRUE))

  # remove old data, if present
  if(rmoldh5){
    if(verbose){message("Removing old data...")}
    for(d in dsnl){
      suppressMessages(try(rhdf5::h5delete(dbn, d), silent = TRUE))
      suppressMessages(try(rhdf5::h5delete(dbn, paste0(d, ".colnames")), silent = TRUE))
      suppressMessages(try(rhdf5::h5delete(dbn, paste0(d, ".rownames")), silent = TRUE))
    }
  }

  # add red and grn signal tables
  if(verbose){message("Adding and populating data tables to HDF5 database")}
  h5_addtables(dbn = dbn, fnl = fnl, dsnl = c("redsignal", "greensignal"),
               rmax = rmax, cmax = cmax, nr.inc = ngsm.block,
               fnpath = fnpath, verbose = verbose)
  if(verbose){message("Finished adding red and green channel data.")}

  # optionally add meth, unmeth, and noobbeta tables
  if(newtables){
    if(verbose){message("Adding new data sets from h5 signal tables...")}
    h5_newtables(dbn, verbose = verbose)
  }

  # optionally add sample metadata
  # note: FINISH THIS PART!
  if(addmd & !is.null(mdpath)){
    if(verbose){message("Adding sample metadata to db...")}
    h5_addmd(dbn, mdpath, verbose = verbose)
  } else{
    if(verbose){message("Failed adding metadata.",
                        "Option addmd is TRUE, but no mdpath specified!",
                        "Continuing..")}
  }

  # finally, close open connections
  rhdf5::h5closeAll()

  if(verbose){message("Finished all processes. Returning.")}
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
#' @param dsn Name or group path to HDF5 dataset containing postprocessed metadata.
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
#' Note, this data should be a data.frame with certain data specified by colnames "basename" and "gsm".
#'
#' @param phenopath Path to properly formatted pheno data file (see details).
#' @param pdat Object containing sample metadata/phenodata (rows = samples, NULL by default).
#' @param se SummarizedExperiment object.
#' @param verbose Whether to display status messages.
#' @return SummarizedExperiment object with appended phenotype data.
#' @export
se_addpheno <- function(phenopath, pdat = NULL, se, verbose = TRUE){
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
  identical(mdp$gsm, gsmv)
  # add valid basenames
  mdp$basename <- colnames(se)
  rownames(mdp) <- colnames(se)
  identical(rownames(mdp), colnames(se))
  # append pheno to se object
  minfi::pData(se) <- DataFrame(mdp)
  return(se)
}

#' Use DelayedArray function to store H5 SummarizedExperiment directory
#'
#' Use DelayedArray function to store H5 SummarizedExperiment directory from h5 db. Handles 3 classes of Summarized Experiment objects.
#'
#' @param newfnstem Filename stem of H5-SE directory.
#' @param version Version of object, to be appended to fnstem.
#' @param ts Timestamp of object, to be appended to fnstem.
#' @param make.from.rg Whether to make new data from red and grn signal (TRUE/FALSE).
#' @param se Type of SummarizedExperiment to parse (either `rg` for RGChannelSet, `gr` for GenomicRatioSet, or `gm` for GenomicMethylSet).
#' @param dbn Name or path to h5 db.
#' @param dsn.data1 Name of first dataset in h5 db (required).
#' @param dsn.data2 Name of second dataset in h5 db (required if se either rg or gm).
#' @param phenopath Sample metadata file (columns are fields, rows are samples).
#' @param dsn.rn Name of object in h5 db with the rownames (GSM or sample basenames).
#' @param dsn.cn Optional name of h5db object containing the column names (cg ids/addresses, default: NULL).
#' @param semd Metadata for new H5-SE object. Should be a list including timestamp, version, description, etc.
#' @param verbose Whether to show verbose status messages.
#' @param replace.opt Whether to replace/overwrite any existing H5SE directories of the same name as `newfn`.
#' @return List of index blocks of min length `slength`/`bsize`
#' @export
make_h5se <- function(dbn, newfnstem, version, ts, make.from.rg = FALSE,
                      dsn.data1, dsn.md = "mdpost", mdpath = NULL,
                      se = c("rg", "gr", "gm"), dsn.rn = "redsignal.rownames",
                      dsn.cn = "redsignal.colnames",
                      semd = list("title" = "Recount Methylation H5-SE Object",
                                  "version" = version,
                                  "timestamp" = ts),
                      dsn.data2 = NULL, addpheno = FALSE, phenopath = NULL,
                      verbose = TRUE, replace.opt = TRUE){
  # Creates an SE-H5 object from an HDF5 db
  # make the new filename
  newfn <- paste(newfnstem, gsub("\\.", "-", version), ts, sep = "_")
  # check the specified se
  if(length(se) > 1){stop("Specify a single se set to process per run.")}
  # anno for se sets
  if(verbose){message("Setting annotation info...")}
  anno = c("IlluminaHumanMethylation450k", "ilmn12.hg19")
  names(anno) = c("array", "annotation")
  # get granges object
  if("gr" %in% se | "gm" %in% se){
    if(verbose){message("Getting the GRanges object...")}
    mset <- get(data(MsetEx, package = "minfiData"))
    mrset <- minfi::mapToGenome(mset)
    grcg <- GenomicRanges::granges(mrset)
  }
  # load data table
  if(verbose){message("Getting dsn.data1...")}
  ldat <- list()
  nb <- HDF5Array::HDF5Array(dbn, dsn.data1)
  rn <- rhdf5::h5read(dbn, dsn.rn)
  cn <- rhdf5::h5read(dbn, dsn.cn)
  rn <- rn[1:nrow(nb)] # subset rows
  cn <- cn[1:ncol(nb)] # subset cols
  nb <- nb[1:length(rn), 1:length(cn)] # trim data
  rownames(nb) <- as.character(rn)
  colnames(nb) <- as.character(cn)
  nb <- t(nb) # transpose so cols are samples, rows are assays
  ldat[[dsn.data1]] <- nb
  # sanity checks and parse data2
  if(!is.null(dsn.data2)){
    if(verbose){message("Getting dsn.data2...")}
    nb <- HDF5Array::HDF5Array(dbn, dsn.data2)
    nb <- nb[1:length(rn), 1:length(cn)] # trim data
    rownames(nb) <- as.character(rn)
    colnames(nb) <- as.character(cn)
    nb <- t(nb) # transpose
    if(!identical(nrow(nb), nrow(ldat[[1]])) |
       !identical(ncol(nb), ncol(ldat[[1]]))){
      stop("Matrix dsn.data2 not similar dim to dsn.data1!")
    }
    ldat[[dsn.data2]] <- nb
  } else if (se %in% c("rg", "gm")){
    stop("Must provide dsn.data2 for se as rg or gm!")
  }
  # get probe ids or addresses
  if(verbose){message("Getting probe ids/addresses...")}
  if(se %in% c("gm", "gr")){
    man.package <- "IlluminaHumanMethylation450kanno.ilmn12.hg19"
    man <- get(data(Manifest, package = man.package))
    cgrn <- rownames(man)
    ldat <- lapply(ldat, function(x){
      rownames(x) <- cgrn
      return(x)
    })
  }
  # make the new H5-SE set(s)
  if(verbose){message("Making the new se object...")}
  # verify se md is list
  if(!class(semd) == "list"){
    semd <- list(semd)
  }
  # make the se object
  if("rg" %in% se){
    if(verbose){message("Making RGChannelSet...")}
    gri <- minfi::RGChannelSet(Red = ldat[[1]],
                               Green = ldat[[2]],
                               anno = anno)
    metadata(gri) <- semd
  } else if ("gr" %in% se){
    if(verbose){message("Making GenomicRatioSet...")}
    gri <- minfi::GenomicRatioSet(gr = grcg,
                                  Beta = ldat[[1]],
                                  anno = anno)
    metadata(gri) <- semd
  } else{
    if(verbose){message("Making GenomicMethylSet...")}
    gri <- minfi::GenomicMethylSet(gr = grcg,
                                   Meth = ldat[[1]],
                                   Unmeth = ldat[[2]],
                                   anno = anno)
    metadata(gri) <- semd
  }
  # append pheno data
  if(addpheno){
    if(is.null(phenopath)){
      message("No phenopath provided, checking ",
              "HDF5 db for sample metadata...")
      if(dsn.md %in% h5ls(dbn) &
         paste0(dsn.md, ".colnames") %in% h5ls(dbn)){
        if(verbose){message("Sample metadata detected in dbn. ",
                            "Adding sample metadata as pheno data...")}
        mdp <- data_mdpost(dbn, dsn.md)
        gri <- se_addpheno(pdat = mdp, se = gri)
      }
      message("Couldn't add pheno data!",
              " Specify phenopath or ensure",
              " the dsn.md entity exists in dbn.",
              " Continuing...")
    } else{
      if(verbose){message("Adding sample metadata from phenopath...")}
      gri <- se_addpheno(phenopath, se = gri)
    }
  }
  # start the run and save the new H5SE set
  t1 <- Sys.time()
  if(verbose){message("Starting process to make new file ", newfn, "...")}
  HDF5Array::saveHDF5SummarizedExperiment(gri,
                                          dir = newfn,
                                          replace = replace.opt)
  if(verbose){message("Save complete, time elapsed:", Sys.time() - t1)}
  return(NULL)
}

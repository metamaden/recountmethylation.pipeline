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
#' @param platform The DNAm array platform of files to read/write (either 
#' "hm450k" or "epic").
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
dtables_rg <- function(platform = c("hm450k", "epic"), version, timestamp, 
  verbose = TRUE, gsmint = 60, overwrite = TRUE, fnstem = "mdat.compilation", 
  sepval = " ", idatspath = file.path("recount-methylation-files", "idats"),
  destpath = file.path("recount-methylation-analysis",
    "files", "mdata", "compilations")){
  idatinfo <- dt_checkidat(idatspath = idatspath, verbose = verbose)
  hlinkv <- idatinfo[["hlinkv"]]; gsmu <- idatinfo[["gsmu"]]
  if(verbose){message("Found ", length(gsmu), " GSM IDs with valid idats.")}
  gsmii <- getblocks(length(gsmu), gsmint) # get id vectors list
  if(verbose){message("Making new data tables...")}
  dtinfo <- dt_makefiles(platform = platform, hlinkv = hlinkv, 
    idatspath = idatspath, destpath = destpath, version = version, nts=timestamp, 
    overwrite = overwrite, sepval = sepval, verbose = verbose)
  if(verbose){message("Wrote data with ", dtinfo[["num.assays"]], " assays.")}
  dtcond <- dtinfo[["dtcond"]]; num.assays = dtinfo[["num.assays"]]
  if(dtcond){red.path <- dtinfo[["reds.path"]];grn.path <- dtinfo[["grns.path"]]
    if(verbose){message("Appending new data for ", length(gsmii)," chunks...")}
    if(platform == "hm450k"){
      require(minfiData); data(RGsetEx); rgi <- RGsetEx
      } else if(platform == "epic"){
        require(minfiDataEPIC); data(RGsetEPIC); rgi <- RGsetEPIC
        } else(stop("Error, invalid platform provided."))
    tt <- Sys.time(); probesv <- rownames(rgi)
    for(i in 1:length(gsmii)){hlinkvi <- hlinkv[gsmii[[i]]]
      dt_write_rg(probesv = probesv, hlinkv = hlinkvi, 
        idatspath = idatspath, reds.path = red.path, grns.path = grn.path, 
        verbose = verbose, num.assays = num.assays)
      if(verbose){message("Finished chunk ", i , " time: ", Sys.time() - tt)}
    }
  } else{stop("Problem encountered handling data tables.")}
    if(verbose){message("Successfully completed data tables.")}; return(NULL)
  }

#' Make and manage data.table files for red and grn signal
#'
#' Handles options for signal tables, returns paths used to write data chunks.
#' @param platform The DNAm array platform of files to read/write (either 
#' "hm450k" or "epic").
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
dt_makefiles <- function(platform = c("hm450k", "epic"), hlinkv, idatspath, 
  destpath, version, nts, overwrite = TRUE, fnstem = "mdat.compilation", 
  sepval = " ", verbose = TRUE){
  version.fn <- gsub("\\.", "-", version)
  reds.fn <- paste("redsignal", nts, version.fn, sep = "_")
  reds.fn <- paste(reds.fn, fnstem, sep = ".")
  grns.fn <- paste("greensignal", nts, version.fn, sep = "_")
  grns.fn <- paste(grns.fn, fnstem, sep = ".")
  reds.path = file.path(destpath, reds.fn)
  grns.path = file.path(destpath, grns.fn); cn = c("gsmi")
  rgi = minfi::read.metharray(c(file.path(idatspath, hlinkv[1:2])),force=TRUE)
  if(verbose){message("Getting platform data...")}
  if(platform == "hm450k"){
    require(minfiData); data(RGsetEx); rgi <- RGsetEx
    } else if(platform == "epic"){
      require(minfiDataEPIC); data(RGsetEPIC); rgi <- RGsetEPIC
      } else{stop("Error, invalid platform provided.")}
  rgcni = colnames(t(minfi::getRed(rgi)));rgcn = matrix(c(cn, rgcni),nrow=1)
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

#' Validate idats at path, returns gsmu and hlinkv
#' 
#' Checks for valid idat file pairs, checks for latest timestamps, 
#' files with matching timestamps, and valid red/grn IDAT file pairs, 
#' gets valid gsms idats dir.
#' @param idatspath File path to directory containing IDATs.
#' @param verbose Whether to return status messages (boolean, TRUE).
#' @return list containing gsmu and hlinkv
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
  gpath.ts <- c(); gpath.gid <- c(); gidv <- gsub("\\..*", "", gpath)
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
  }; hlinkv <- gpath.ts; gsmu <- gsub("\\..*", "", hlinkv)
  lr <- list("gsmu" = gsmu, "hlinkv" = hlinkv); return(lr)
}

#' Writes chunk of red and grn signal data
#'
#' @param probesv Vector of probe labels.
#' @param hlinkv Vector of GSM IDAT hlink file names, or file basenames.
#' @param idatspath Path to idat files directory.
#' @param reds.path Path to new red signal data table.
#' @param grns.path Path to new grn signal data table.
#' @param min.cols Minimum data columns to write DNAm data (integer, )
#' @param sepval Separator symbol for data tables.
#' @param num.assays Number of assays to expect before writing (1052641 for 
#' epic, 622399 for hm450k).
#' @param verbose Whether to include verbose messages.
#' @return NULL, writes data chunks as side effect
#' @export
dt_write_rg <- function(probesv, hlinkv, idatspath, reds.path, grns.path, 
  num.assays = c(622399, 1052641), sepval = " ", verbose = TRUE){
  if(verbose){message("Reading data...")};
  pathl=unique(file.path(idatspath, hlinkv))
  rgi=try(minfi::read.metharray(pathl,force=TRUE))
  if(class(rgi) == "RGChannelSet"){if(verbose){message("Getting signal data..")}
    red.dat <- minfi::getRed(rgi); grn.dat <- minfi::getGreen(rgi)
    probesv.out <- probesv[!probesv %in% rownames(rgi)]
    if(length(probesv.out) > 0){
      outdat <- rep(rep("NA", length(probesv.out)), ncol(red.dat))
      outm <- matrix(outdat, ncol = ncol(red.dat));rownames(outm) <- probesv.out
      red.dat <- rbind(red.dat, outm); grn.dat <- rbind(grn.dat, outm)
    }
    reorder.red <- order(match(rownames(red.dat), probesv))
    reorder.grn <- order(match(rownames(grn.dat), probesv))
    rdat.order <- red.dat[reorder.red,, drop = FALSE]
    gdat.order <- grn.dat[reorder.grn,, drop = FALSE]
    cond <- identical(rownames(rdat.order), probesv) & 
      identical(rownames(gdat.order), probesv)
    if(cond){if(verbose){message("Appending new data...")}
      redi = matrix(c(colnames(rgi),t(rdat.order)), ncol = nrow(rdat.order)+1)
      grni = matrix(c(colnames(rgi),t(gdat.order)), ncol = nrow(gdat.order)+1)
      data.table::fwrite(redi, reds.path, sep = sepval, append = TRUE)
      data.table::fwrite(grni, grns.path, sep = sepval, append = TRUE)
    } else{if(verbose){message("Error matching probe vectors. Skipping...")}}
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

#' @name make_h5db_rg
#' @rdname make_h5db_rg
#' Make and populate a new HDF5 database with red and green signal, and metadata
#'
#' Add signal data (red and green channel) to the HDF5 database.
#' @param dbfnstem Stem of filename for h5 database (to which ts and version 
#' appended)
#' @param dbpath Database file path to write to.
#' @param version Version of new database file.
#' @param ts Timestamp of new database file.
#' @param fnl Vector of signal tables containing data to be added.
#' @param fnpath Path to dir containing files in fnl.
#' @param dsnl Vector of data set names in HDF5 database to be populated.
#' @param rmax Total rows to append to data sets, reflecting total samples.
#' @param cmax Total columns to append to data sets, reflecting total assays or 
#' probes. For raw red/grn signal, can be either 622399 (HM450K platform) or 
#' 1052641 (EPIC platform). 
#' @param newtables Whether to also add new data tables (noob-norm. Beta-values, 
#' meth. and unmeth. signal).
#' @param mdpath If addmd, the path to the metadata file to load.
#' @param ngsm.block Number of GSMs (samples) per process block (default 50, 
#' passed to `h5_newtables()`).
#' @param verbose Whether to show verbose status updates.
#' @param dsnl Vector of new h5 data set objects in the h5db object, 
#' corresponding (1:1) to the files declared in fnl.
#' @return Populates the HDF5 database
#' @export
make_h5db_rg <- function(dbfnstem, dbpath, version, ts, fnpath,
  fnl, dsnl = c("redsignal", "greensignal"), rmax = 35300, 
  cmax = c(622399, 1052641), rmoldh5 = TRUE, newtables = FALSE, mdpath = NULL,
  ngsm.block = 50, verbose = TRUE){
  require(rhdf5);fn <- paste(dbfnstem, ts, gsub("\\.", "-", version), sep = "_")
  dbn <- file.path(dbpath, paste(fn, "h5", sep = "."))
  if(verbose){message("Making new h5 db file: ", dbn)}
  suppressMessages(try(rhdf5::h5createFile(dbn), silent = TRUE))
  if(rmoldh5){if(verbose){message("Removing old data...")}
    for(d in dsnl){suppressMessages(try(rhdf5::h5delete(dbn, d), silent = TRUE))
      suppressMessages(try(rhdf5::h5delete(dbn, paste0(d, ".colnames")),
        silent = TRUE))
      suppressMessages(try(rhdf5::h5delete(dbn, paste0(d, ".rownames")), 
        silent = TRUE))}}
  # add red and grn signal tables
  if(verbose){message("Adding and populating data tables to HDF5 database")}
  for(di in seq(length(dsnl))){
    h5_add_tables(dsn = dsnl[di], fnl = fnl, fnpath = fnpath, dbn = dbn,
      rmax = rmax, cmax = cmax, ngsm.block = ngsm.block, verbose = verbose)}
  if(verbose){message("Finished adding red and green channel data.")}
  if(!is.null(mdpath)){if(verbose){message("Adding sample metadata to db...")}
    h5_addmd(dbn, mdpath, verbose = verbose)
  } else{if(verbose){message("No mdpath specified!")}}
  rhdf5::h5closeAll();if(verbose){message("Finished all processes. Returning.")}
  return(NULL)
}
#' @rdname getexpfractions
#' @examples
#' #
#' @export
h5_add_tables = function(dsn,fnl,fnpath,dbn,rmax,cmax,ngsm.block,verbose){
  if(dsn == "redsignal"){fnread = fnl[grepl("red", fnl)]}else{
    fnread <- fnl[grepl("green", fnl)]}
  if(verbose){message("for dsn ", dsn, " using fnread ", fnread, "...")}
  tt<-Sys.time();cnn=paste0(dsn,".colnames");rnn=paste0(dsn, ".rownames")
  rn=cn=c();con<-file(paste(fnpath, fnread, sep = "/"), "r")
  if(verbose){"Making dataset...."}
  rhdf5::h5createDataset(dbn,dsn,dims=c(rmax, cmax),
                         maxdims=c(rhdf5::H5Sunlimited(),rhdf5::H5Sunlimited()),
                         storage.mode="double",level=5,chunk=c(100, 5000))
  if(verbose){message("Adding colnames...")}
  cn = unlist(strsplit(readLines(con, n = 1), " ")); cn = cn[2:length(cn)]
  cn = gsub("\n", "",gsub('\"', '', cn[1:cmax]))
  rhdf5::h5createDataset(dbn,cnn,dims=length(cn),storage.mode="character", 
    maxdims = c(rhdf5::H5Sunlimited()), level = 5, chunk = c(20), size = 256)
  rhdf5::h5write(cn, file = dbn, name = cnn, index = list(1:length(cn)))
  if(verbose){message("Making rownames dataset....")}
  rhdf5::h5createDataset(dbn, rnn, dims = rmax, 
        maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character", 
        level = 5, chunk = c(20), size = 256)
  if(verbose){message("Working on sample indices...")}
  # nri = getblocks(rmax, ngsm.block)
  for(ni in seq(1, rmax, ngsm.block)){if(verbose){message("Reading lines...")}
    start.index <- ni; end.index <- ni + ngsm.block -1
    dati = matrix(unlist(strsplit(readLines(con, n = ngsm.block), " ")),
      nrow = ngsm.block, byrow = TRUE)
    dati.filt <- dati[,c(2:ncol(dati))]; class(dati.filt) <- "numeric"
    wt <- try(rhdf5::h5write(dati.filt, file = dbn, name = dsn,
        index = list(start.index:end.index, 1:cmax)))
    if(!class(wt) == "type-error"){rn = c(rn, dati[,1])
      if(verbose){message("Successful table write. Writing new row entries...")}
      rhdf5::h5write(dati[,1],file=dbn,name=rnn,index=list(start.index:end.index))
    } else{if(verbose){message("Error writing new table data. Skipping...")}}
    message("For ds ", dsn,", finished reading index ", start.index, 
      " to ", end.index, ", time; ", Sys.time() - tt)
  };if(verbose){message("Completed writing data for ds, ", dsn)}
  rhdf5::h5closeAll();if(verbose){message("Finished adding signal dataset")}
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

#' make_h5se_rg
#'
#' Make an RGChannelSet HDF5-SummarizedExperiment object from an HDF5 database.
#'
#' @param max.sample Max index of samples to write to new data (integer). Setting
#' this will trim the samples on reading the HDF5 db. This should reflect the 
#' final successfully written sample index, and/or the final unique sample id 
#' in the rownames vector that is read.
#' @param platform Array platform (either "hm450k" or "epic").
#' @param version The data file version.
#' @param ts Data file NTP timestamp (character).
#' @param dbn Name/path of HDF5 database file.
#' @param newfnstem File name stem for new database file.
#' @param dsnv Table names vector.
#' @param add.metadata Whether to add samples metadata (boolean, FALSE).
#' @param mdpath Path to metadata pile (if NULL, search the HDF5 database).
#' @param dsn.md Metadata table name, used if add.metadata is TRUE.
#' @param dsn.md.cn Metadata columns object, used if add.metadata is TRUE.
#' @param verbose Whether to show status messages (boolean, TRUE).
#' @param replace.opt Whether to replace/overwrite old H5SE object (TRUE).
#' @param dsn.rnv Names vector of rownames objects in HDF5 database.
#' @param dsn.cnv Names vector of colnames objects in HDF5 database.
#' @param semd Metadata vector for new H5SE object to be written.
#' @returns New h5se object.
#' @export
make_h5se_rg <- function(max.sample, platform = c("hm450k", "epic"),
  version = "0.0.1", ts = "1590090412", dbn = "remethdb_1590090412_0-0-1.h5", 
  newfnstem = "remethdb_h5se-rg", dsnv = c("redsignal", "greensignal"),
  add.metadata=FALSE,mdpath=NULL,dsn.md="mdpost",
  dsn.md.cn=paste0(dsn.md,".colnames"), verbose = TRUE, replace.opt = TRUE,
  dsn.rnv = c(paste0(dsnv[1], ".rownames"), paste0(dsnv[2], ".rownames")),
  dsn.cnv = c(paste0(dsnv[1], ".colnames"), paste0(dsnv[2], ".colnames")),
  semd=list("title"="RGChannelSet HDF5-SummarizedExperiment object",
    "preprocessing"="raw")){
  newfn <- paste(newfnstem, platform, gsub("\\.", "-", version), ts, sep = "_")
  if(verbose){message("Setting annotation info...")}
  if(platform == "hm450k"){
    anno = c("IlluminaHumanMethylation450k", "ilmn12.hg19")
    } else if(platform == "epic"){
      anno = c("IlluminaHumanMethylationEPIC", "ilm10b2.hg19")
      } else{stop("Error, didn't recognize provided platform.")}
  names(anno) = c("array", "annotation")
  if(verbose){message("Getting h5 datasets...")};ldat <- list()
  for(i in 1:length(dsnv)){nb <- HDF5Array::HDF5Array(dbn, dsnv[i])
    rn <- rhdf5::h5read(dbn, dsn.rnv[i]);if(max.sample){rn <- rn[1:max.sample]}
    cn <- rhdf5::h5read(dbn, dsn.cnv[i])
    nb <- nb[c(1:length(rn)),c(1:length(cn))]
    rownames(nb) <- as.character(rn); colnames(nb) <- as.character(cn)
    nb <- t(nb);ldat[[dsnv[i]]] <- nb
  };if(verbose){message("Making the H5SE set...")}
  gri <- minfi::RGChannelSet(Red = ldat[[1]], Green = ldat[[2]], anno = anno)
  S4Vectors::metadata(gri) <- semd # H5SE object metadata
  if(add.metadata){if(verbose){message("Adding metadata...")}
    if(is.null(mdpath)){
      if(dsn.md %in% h5ls(dbn)$name & dsn.md.cn %in% h5ls(dbn)$name){
        if(verbose){message("Adding metadata from dbn...")}
        mdp <- recountmethylation::data_mdpost(dbn, dsn.md) 
        gri <- se_addpheno(pdat = mdp, se = gri)
      } else{message("Couldn't find md objects in h5 db, skipping..")}
    } else{if(verbose){message("Adding metadata from phenopath...")}
      gri <- try(se_addpheno(mdpath, se = gri))
    }
  }else{if(verbose){message("Skipping metadata addition step...")}}
  if(verbose){message("Adding signals data to file ", newfn)};t1 <- Sys.time()
  HDF5Array::saveHDF5SummarizedExperiment(gri,dir=newfn,replace=replace.opt)
  if(verbose){message("Data additions complete, time elapsed:",Sys.time()-t1)}
  return(NULL)
}

#' Make an HDF5 h5 file with MethylSet data, from an RGChannelSet 
#' h5se file
#'
#' Makes an HDF5 h5 database file with raw/unnormalized methylated and 
#' unmethylated signal, including rows and columns as datasets.
#'
#' @param dbn Path to the RGChannelSet H5SE data directory.
#' @param version Version for new h5 file.
#' @param ts Timestamp for new h5 file.
#' @param num.samp Number of samples in HDF5 dataset.
#' @param blocksize Size of blocks/number of samples to process at a time 
#' (default 65).
#' @param platform Array platform (either "hm450k" or "epic").
#' @param newfnstem New filename stem for h5 file.
#' @param verbose Whether to return verbose messages (default TRUE).
#' @param dsv Vector of new datasets to add to h5 HDF5 database.
#' @param replace.opt Whether to replace h5 files with duplicate filename 
#' (default TRUE).
#' @returns Path to new h5 object.
#' @export
make_h5_gm <- function(dbn, version, ts, num.samp, blocksize = 65,
  platform = c("hm450k", "epic"), newfnstem = "remethdb_h5-gm", verbose = TRUE, 
  dsv = c("meth", "unmeth"), replace.opt = TRUE){
  if(platform == "hm450k"){num.assays = 485512} else if(platform == "epic"){
    num.assays = 866836} else{stop("Error, didn't recognize platform.")}
  rg <- HDF5Array::loadHDF5SummarizedExperiment(dbn)
  h5dbn <- paste(newfnstem, platform, gsub("\\.", "-", version), ts, sep = "_")
  h5dbn <- paste0(c(h5dbn, "h5"), collapse = ".");message("Making h5 db ",h5dbn)
  if(!file.exists(h5dbn)){rhdf5::h5createFile(h5dbn)}
  for(d in dsv){rhdf5::h5createDataset(h5dbn, d,dims=list(num.assays,num.samp), 
      maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()), 
      storage.mode = "double", level = 5, chunk = c(1000, 50))}
  rhdf5::h5createDataset(h5dbn, "rownames", dims = list(num.assays), 
      maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character", size = 256, 
      level = 5, chunk = c(1000))
  rhdf5::h5createDataset(h5dbn, "colnames", dims = list(num.samp), 
    maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character", size = 256, 
    level = 5, chunk = c(1000))
  if(verbose){message("Getting blocks of sample indices.")}
  blocks <- getblocks(ncol(rg), bsize = blocksize)
  b<-1;cindices<-blocks[[b]];rgi<-rg[,cindices];gm<-minfi::preprocessRaw(rgi)
  meth<-as.matrix(minfi::getMeth(gm));unmeth<-as.matrix(minfi::getUnmeth(gm))
  class(meth) <- class(unmeth) <- "numeric"
  lindex <- list(1:num.assays, cindices[1]:cindices[length(cindices)])
  if(verbose){message("Writing first data blocks...")}
  rhdf5::h5write(meth, file = h5dbn, name = "meth", index = lindex)
  rhdf5::h5write(unmeth, file = h5dbn, name = "unmeth", index = lindex)
  if(verbose){message("Writing rows and cols.")}
  rhdf5::h5write(colnames(rg), file = h5dbn, name = "colnames", index = list(1:ncol(rg)))
  rhdf5::h5write(rownames(gm), file = h5dbn, name = "rownames", index = list(1:num.assays))
  message("Writing next data blocks...");t1 <- Sys.time()
  for(b in 2:length(blocks)){
    cindices <- blocks[[b]]; gm <- minfi::preprocessRaw(rg[,cindices])
    meth<-as.matrix(minfi::getMeth(gm));unmeth<-as.matrix(minfi::getUnmeth(gm))
    write.indices <- cindices[1]:cindices[length(cindices)]
    lindex <- list(1:num.assays, write.indices)
    rhdf5::h5write(meth, index = lindex, file = h5dbn, name = "meth")
    rhdf5::h5write(unmeth, index = lindex, file = h5dbn, name = "unmeth")
    if(verbose){message("Finished block ", b, ", time = ", Sys.time() - t1)}
  };message("Finished writing data blocks. Returning h5 name.")
  return(h5dbn)
}


#' Make an HDF5 h5 file with GenomicRatioSet data, from an RGChannelSet 
#' h5se file
#'
#' Makes an HDF5 h5 file containing noob-normalized Beta-values, including
#' rows and columns as separate datasets.
#' @param dbn Path to the RGChannelSet H5SE data directory.
#' @param version Version for new h5 file.
#' @param ts Timestamp for new h5 file.
#' @param num.samp Number of samples in HDF5 dataset.
#' @param blocksize Size of blocks/number of samples to process at a time 
#' (default 65).
#' @param platform Array platform (either "hm450k" or "epic").
#' @param newfnstem New filename stem for h5 file.
#' @param verbose Whether to return verbose messages (default TRUE).
#' @param replace.opt Whether to replace h5 files with duplicate filename 
#' (default TRUE).
#' @returns Path to new h5 object.
#' @export
make_h5_gr <- function(dbn, version, ts, num.samp, blocksize = 65, verbose=TRUE,
  platform=c("hm450k","epic"),newfnstem="remethdb_h5-gr",replace.opt=TRUE){
  if(platform == "hm450k"){num.assays = 485512} else if(platform == "epic"){
  num.assays = 866836} else{stop("Error, didn't recognize platform.")}
  rg <- HDF5Array::loadHDF5SummarizedExperiment(dbn)
  h5dbn <- paste(newfnstem, platform, gsub("\\.", "-", version), ts, sep = "_")
  message("Making h5 db: ", h5dbn);h5dbn <- paste0(c(h5dbn, "h5"), collapse=".")
  if(!file.exists(h5dbn)){rhdf5::h5createFile(h5dbn)}
  if(verbose){message("Making main noobbeta data table...")}
  rhdf5::h5createDataset(h5dbn, "noobbeta", dims = list(num.assays, num.samp), 
    maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()), 
    storage.mode = "double", level = 5, chunk = c(1000, 50))
  if(verbose){message("Making colnames data table...")}
  rhdf5::h5createDataset(h5dbn, "colnames", dims = list(num.samp), 
  maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character", size = 256, 
  level = 5, chunk = c(1000))
  if(verbose){message("Making, writing rownames data table...")}
  rhdf5::h5createDataset(h5dbn, "rownames", dims = list(num.assays), 
    maxdims = c(rhdf5::H5Sunlimited()), storage.mode = "character", size = 256, 
    level = 5, chunk = c(1000))
  rhdf5::h5write(rownames(gr),file=h5dbn,name="rownames",
          index=list(1:num.assays))
  blocks <- getblocks(ncol(rg), bsize = blocksize)
  message("Writing data blocks.");t1 <- Sys.time()
  for(b in 1:length(blocks)){cindices <- blocks[[b]];
      if(verbose){message("Beginning block ", b, "/", length(blocks), "...")} 
      lindex <- list(row.indices, cindices);rgi <- rg[,cindices]
      tryCatch(gr <- minfi::preprocessNoob(rgi),
          error = function(e){
              message("Error with preprocessNoob, trying samples in block...")
              for(si in seq(length(cindices))){
                  message("Attempting to write data for sample ", si, "...")
                  lindex <- list(row.indices, cindices[si])
                  tryCatch(gr <- minfi::preprocessNoob(rgi[,si]),
                      error = function(e){message("Error at ",si,", skipping...")}
                  )
                  if(class(gr) == "MethylSet"){
                      nbeta <- as.matrix(minfi::getBeta(gr));class(nbeta)<-"numeric"
                      rhdf5::h5write(nbeta,file=h5dbn,name="noobbeta",index=lindex)
                      message("Wrote data for sample ", si, ", continuing...")
                  }
              }
          }
      )
      if(class(gr) == "MethylSet"){
          nbeta <- as.matrix(minfi::getBeta(gr));class(nbeta) <- "numeric"
          rhdf5::h5write(nbeta, file=h5dbn, name="noobbeta", index=lindex)
      }
      if(verbose){message("Finished block ", b, "/", length(blocks), 
          ", time = ", Sys.time() - t1)}
  };message("Finished writing data blocks. Returning h5 name.");return(h5dbn)
}

#' Make an H5SE MethySet object from a MethySet HDF5 database file
#' 
#' @param dbn Name of the h5 file to read data from.
#' @param version Version for new filenames.
#' @param ts NTP timestamp.
#' @param platform Array platform (either "hm450k" or "epic").
#' @param replaceopt Whether to overwrite existing file of same name as new file
#' (default = TRUE).
#' @param verbose Whether to return verbose messages (default TRUE).
#' @param add.metadata Whether to add samples metadata (boolean, FALSE).
#' @param pdata Metadata object obtained from pData(se).
#' @param newdnstem Stem for new h5se object file.
#' @return NULL, generates a new h5se file.
#' @export
make_h5se_gm <- function(dbn, version, ts, platform = c("hm450k", "epic"),   
  replaceopt = TRUE, verbose = TRUE, add.metadata = FALSE, pdata = NULL, 
  newdnstem = "remethdb_h5se-gm", 
  semd=list("title"="GenomicMethylSet HDF5-SummarizedExperiment object",
    "preprocessing"="raw")){
  newfn <- paste(newdnstem, platform, gsub("\\.", "-", version), ts, sep="_")
  if(verbose){message("Making new H5SE database: ", newfn)}
  if(verbose){message("Getting granges...")}
  if(platform == "hm450k"){require(minfiData)
    ms <- minfi::mapToGenome(get(data("MsetEx")))
    } else if(platform == "epic"){require(minfiDataEPIC)
      ms <- minfi::mapToGenome(get(data("MsetEPIC")))
      } else{stop("Error, didn't recognize platform.")}
  anno <- minfi::annotation(ms);gr <- GenomicRanges::granges(ms)
  message("Reading datasets.")
  meth <- HDF5Array::HDF5Array(dbn, "meth")
  unmeth <- HDF5Array::HDF5Array(dbn, "unmeth")
  cn <- rhdf5::h5read(dbn, "colnames");rn <- rhdf5::h5read(dbn, "rownames")
  rownames(meth) <- rownames(unmeth) <- as.character(rn)
  colnames(meth) <- colnames(unmeth) <- as.character(cn)
  gr <- gr[order(match(names(gr), rownames(meth)))]
  icond <- identical(rownames(meth), names(gr))
  if(!icond){stop("Problem matching gr names to h5 meth rows")}
  icond <- identical(rownames(unmeth), names(gr))
  if(!icond){stop("Problem matching gr names to h5 unmeth rows")}
  if(add.metadata & !is.null(pdata)){
    icond <- identical(colnames(meth), rownames(pdata))
    if(!icond){stop("Problem matching pdata rows to h5 meth cols")}
    icond <- identical(colnames(unmeth), rownames(pdata))
    if(!icond){stop("Problem matching pdata rows to h5 unmeth cols")}
  }
  if(verbose){message("Defining new SummarizedExperiment object...")}
  gm <- minfi::GenomicMethylSet(gr = gr, Meth = meth, Unmeth = unmeth, 
                          annotation = anno)
  S4Vectors::metadata(gm) <- semd # H5SE object metadata
  if(add.metadata & !is.null(pdata)){message("Adding metadata...")
    if(is.character(pdata)){pdata <- get(load(pdata))};pData(gm) <- pdata}
  message("Writing data to new h5se object. This may take awhile.")
  HDF5Array::saveHDF5SummarizedExperiment(gm, dir = newfn, replace = replaceopt)
  message("Finished writing h5se data.")
  return(NULL)
}

#' Make an h5se object from an h5se gr file
#' 
#' @param version Version for new filenames.
#' @param ts NTP timestamp.
#' @param dbn Path to HDF5 database file to read data from.
#' @param pdata Metadata object obtained from pData(se). 
#' @param replaceopt Whether to overwrite existing file of same name as new file
#' (default = TRUE).
#' @param verbose Whether to return verbose messages (default TRUE).
#' @param add.metadata Whether to add samples metadata (boolean, FALSE).
#' @param pdata Samples metadata as a DataFrame, path to file, or NULL.
#' @param newdnstem Character stme for new h5se object file 
#' ("remethdb_h5se-gr").
#' @return NULL, generates a new h5se file.
#' @export
make_h5se_gr <- function(dbn, version, ts, platform = c("hm450k", "epic"), 
  replaceopt = TRUE, verbose = TRUE, add.metadata = FALSE, pdata = NULL, 
  newdnstem = "remethdb_h5se-gr", 
  semd=list("title"="GenomicMethylSet HDF5-SummarizedExperiment object",
    "preprocessing"="Normalization with out-of-band signal (noob)")){
  newfn <- paste(newdnstem, platform, gsub("\\.", "-", version), ts, sep = "_")
  if(verbose){message("Making new H5SE database: ", newfn)}
  if(platform == "hm450k"){require(minfiData)
    ms <- minfi::mapToGenome(get(data("MsetEx")))
    } else if(platform == "epic"){require(minfiDataEPIC)
      ms <- minfi::mapToGenome(get(data("MsetEPIC")))
      } else{stop("Error, didn't recognize platform.")}
  anno <- minfi::annotation(ms);gr <- GenomicRanges::granges(ms)
  message("Reading datasets.");nbeta <- HDF5Array::HDF5Array(dbn, "noobbeta")
  cn<-rhdf5::h5read(dbn,"colnames");rn<-rhdf5::h5read(dbn,"rownames")
  rownames(nbeta) <- as.character(rn);colnames(nbeta) <- as.character(cn)
  gr <- gr[order(match(names(gr), rownames(nbeta)))]
  icond <- identical(rownames(nbeta), names(gr))
  if(!icond){stop("Problem matching gr names to h5 nbeta rows")}  
  gmi <- minfi::GenomicRatioSet(gr = gr, Beta = nbeta, annotation = anno)
  S4Vectors::metadata(gmi) <- semd # H5SE object metadata
  if(add.metadata & !is.null(pdata)){message("Checking and adding metadata...")
    icond <- identical(colnames(nbeta), rownames(pdata))
    if(!icond){stop("Problem matching pdata rows to h5 nbeta cols")}
    if(is.character(pdata)){pdata <- get(load(pdata))};pData(gmi) <- pdata
  }
  message("Writing new data to h5se object. This may take awhile.")
  HDF5Array::saveHDF5SummarizedExperiment(gmi, dir = newfn, replace = replaceopt)
  message("Finished writing h5se data.");return(NULL)
}

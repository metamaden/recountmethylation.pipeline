#--------------
# OLD FUNCTIONS
#--------------
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
                      rmoldh5 = TRUE, newtables = FALSE,
                      addmd = TRUE, mdpath = NULL,
                      nr.inc = 10, ngsm.block = 50, verbose = TRUE,
                      dsnl = c("redsignal", "greensignal", "mdpost")){
  require(rhdf5)
  # make new h5 db filename
  fn <- paste(dbfnstem, ts, gsub("\\.", "-", version), sep = "_")
  dbn <- paste(fn, "h5", sep = ".")
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
  #if(newtables){
  #  if(verbose){message("Adding new data sets from h5 signal tables...")}
  #  h5_newtables(dbn, verbose = verbose)
  #}
  
  # optionally add sample metadata
  # note: FINISH THIS PART!
  #if(addmd & !is.null(mdpath)){
  #  if(verbose){message("Adding sample metadata to db...")}
  #  h5_addmd(dbn, mdpath, verbose = verbose)
  #} else{
  #  if(verbose){message("Failed adding metadata.",
  #                      "Option addmd is TRUE, but no mdpath specified!",
  #                      "Continuing..")}
  #}
  
  # finally, close open connections
  rhdf5::h5closeAll()
  
  if(verbose){message("Finished all processes. Returning.")}
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
#!/usr/bin/env R

# Author: Sean Maden
# Functions for instance snakemake workflow.

#------------------
# instance metadata
#------------------

#' Get new instance metadata
#'
#' Generates the new metadata, including version and timestamp, for the
#' `recountmethylation` instance. This is called by rule `new_instance_md`.
#' 
#' @param instdir Path to directory to contain all instance metadata.
#' @return NULL, produces the instance metadata file as side effect.
#' @export
new_instance_md <- function(instdir = file.path("recount-methylation-files", "metadata")){
  if(!dir.exists(instdir)){message("Making instdir: ", instdir)
    dir.create(instdir)}
  message("Provide version:"); version <- readLines("stdin", n = 1)
  message("Using ", version, " as the instance version...")
  message("Using `", instdir, "` as the instance directory...")
  if(!dir.exists(instdir)){
    message("Making dir ", instdir, "...");dir.create(instdir)}
  md <- get_metadata("compilation_metadata", version = version)
  # automatically detect platform
  message("Detecting platform from `settings.py`...")
  sett.path <- file.path("recountmethylation_server", "src", "settings.py")
  settrl <- gsub(".* |'", "", readLines(sett.path, n = 31)[31])
  md$accessionID <- settrl
  md$platform <- ifelse(settrl == "GPL13534", "hm450k",
                        ifelse(settrl == "GPL21145", "epic-hm850k",
                               ifelse(settrl == "GPL8490", "hm27k", "NA")))
  message("Saving new metadata...")
  md.fn <- paste0("metadata_v", gsub("\\.", "-", md$version), 
                  "_", md$timestamp, ".rda")
  md.path <- file.path(instdir, md.fn)
  message("Saving instance metadata to ", file.path(md.path))
  save(md, file = md.path); return(NULL)
}

#' Retrieve existing instance metadata
#'
#' User dialogue to handle instance metadata previously generated 
#' with `new_instance_md()`.
#'
#' @param instdir Path to directory to contain all instance metadata.
#' @return Instance metadata as a list containing the version and timestamp.
#' @export
get_data_md() <- function(instdir = "rmp_instance"){
  message("Getting metadata..."); sdv <- as.numeric(list.files(instdir))
  sdv.max <- as.character(sdv[sdv == max(sdv)])
  md.fn <- list.files(file.path(instdir, sdv.max))
  md.fn <- md.fn[grepl("metadata.*", md.fn)[1]]
  message("(Y/N) Use detected md file ", md.fn, "?")
  mdopt <- readLines("stdin", n = 1)
  if(mdopt %in% c("y", "Y", "yes", "Yes", "YES")){
    mdpath <- file.path("rmp_instance", sdv.max, md.fn)
  } else if(mdopt %in% c("n", "N", "No", "no", "NO")){
    message("Provide mdpath: "); mdpath <- readLines("stdin", n = 1)
  } else{stop("Error,invalid input")}
  if(!file.exists(mdpath)){stop("Error, mdpath doesn't exist: ", mdpath)}
  message("Loading metadata from file...")
  md <- get(load(mdpath)); version <- md$version; ts <- md$timestamp
  return(md)
}

#' Check DNAm data contained in an HDF5 database
#'
#' Checks DNAm data stored in an HDF5 database. Uses a GEO lookup to grab 
#' data for data validation.
#' 
#' @param comppath Path to compilation files for the `recountmethylation` instance.
#' @return
#' @export
check_h5_database <- function(comppath = file.path("recount-methylation-files", 
                                                 "compilations")){
  checkopt <- readLines("stdin", n = 1)
  if(mdopt %in% c("y", "Y", "yes", "Yes", "YES")){
    message("Checking h5 files..."); dbn <- cvf; rn <- h5read(dbn, "greensignal.rownames")
    rnf <- rn[!rn == "" & grepl("GSM", rn)] 
    rnff <- rnf[c(1, 2, length(rnf) - 1,length(rnf))]; gsmv <- gsub("\\..*", "", rnff)
    message("Getting IDATs from GEO..."); rg.gds <- gds_idat2rg(gsmv)
    message("Comparing h5 data to downloaded IDATs...")
    gs <- HDF5Array(dbn, "greensignal");rs <- HDF5Array(dbn, "redsignal")
  } else if(mdopt %in% c("n", "N", "No", "no", "NO")){
    message("Skipping h5 file check...")} else{stop("Error,invalid input")}
  return(NULL)
}

#-------------------------
# "run_pipeline" functions
#-------------------------

#' Get compiled data tables of red and green signals
#'
#' Generates the data tables compiling DNAm red and green signals 
#' from downloaded, and validated, IDATs. Makes two tables, one for 
#' each channel type. This is called by the rules `get_rg_compilations` 
#' and `run_dnam_pipeline`.
#'
#' @param files.dname Name of the instance files directory.
#' @param comp.dname Name of the instance compilations directory.
#' @return NULL, generates data tables as side effect
#' @export
get_rg_dtables <- function(files.dname = "recount-methylation-files",
                           comp.dname = "compilations"){
  comp.dpath <- file.path(files.dname, comp.dname)
  if(!dir.exists(comp.dpath)){stop("Error, didn't find compilations dir ", comp.dpath)}
  message("Handling metadata options...")
  md <- rmp_handle_metadata(); if(is.null(md)){stop("Couldn't get metadata...")}
  message("Getting platform info..."); accinfo <- rmp_handle_platform()
  platform <- accinfo[["platform_name"]]; 
  version <- md[["version"]]; ts <- md[["timestamp"]]
  dtables_rg(version, platform = platform, ts, destpath = comp.dpath)
  return(NULL)
}

#' Get HDF5 database of red/green signals, from compilation data tables
#'
#' Generates and HDF5 database of compiled red and green signals. Uses 
#' data tables previously generated using `get_rg_dtables` and `dtables_rg`.
#' This is called by the rules `get_h5db_rg` and `run_dnam_pipeline`.
#' 
#' @param files.dpath Path to instance files directory.
#' @param comp.dname Name of compilations directory.
#' @param ngsm.block Number of GSM IDs per processed data block.
#' @return NULL, generates an h5 rg database file as side effect.
#' @export
get_h5db_rg <- function(files.dpath = "recount-methylation-files", 
                        comp.dname = "compilations", ngsm.block = 50){
  comp.dpath <- file.path(files.dpath, comp.dname)
  if(!dir.exists(comp.dpath)){stop("Error, didn't find compilations dir ", comp.dpath)}
  message("Handling metadata options...")
  md <- rmp_handle_metadata(); if(is.null(md)){stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); accinfo <- rmp_handle_platform()
  message("Checking for signal compilation tables...")
  vform <- gsub("\\.", "-", version)
  lfv <- list.files(comp.dpath);
  lfv <- lfv[grepl(vform, lfv) & grepl(ts, lfv)] # filter on instance metadata
  lfv <- lfv[grepl("\\.mdat\\.compilation$", lfv)] # filter compilations
  lfv.grn <- lfv[grepl("^greensignal.*", lfv)]
  lfv.red <- lfv[grepl("^redsignal.*", lfv)]
  if(length(lfv.grn) == 0){stop("Error, couldn't find green compilation ",
                                "file at ",comp.dpath)}
  if(length(lfv.red) == 0){stop("Error, couldn't find red compilation ",
                                "file at ",comp.dpath)}
  read.path <- write.path <- comp.dpath; fnl <- c(lfv.grn[1], lfv.red[1])
  message("Making HDF5 database file from compilations tables...")
  make_h5db_rg(dbfnstem = "remethdb", dbpath = write.path, version = version, 
               ts = ts, mdpath = NULL, platform = platform, fnpath = read.path, 
               fnl = fnl, ngsm.block = ngsm.block);return(NULL)
}

#' Get HDF5-SummarizedExperiment of red/green signals, from h5 rg data
#'
#' @param files.dpath Path to instance files directory.
#' @param comp.dname Name of compilations directory.
#' @param ngsm.block Number of GSM IDs per processed data block.
#' @return NULL, generates an h5se rg database directory as side effect.
#' @export
get_h5se_rg <- function(files.dpath = "recount-methylation-files", 
                     comp.dname = "compilations", ngsm.block = 50){
  comp.dpath <- file.path(files.dpath, comp.dname)
  if(!dir.exists(comp.dpath)){stop("Error, didn't find compilations dir ", comp.dpath)}
  message("Handling metadata options...")
  md <- rmp_handle_metadata(); if(is.null(md)){stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); accinfo <- rmp_handle_platform()
  message("Checking for valid h5se RGChannelSet file...")
  vform <- gsub("\\.", "-", version);lfv <- list.files(comp.dpath)
  lfv <- lfv[grepl(vform, lfv) & grepl(ts, lfv)] # filter on instance metadata
  lfv <- lfv[grepl(".*h5se_rg.*", lfv)]; 
  if(length(lfv) == 0){
    stop("Couldn't find HDF5 rg database file at: ", comp.dpath, ".\n",
         "Try running rule `get_h5db_rg` first.")
  } else{message("Using HDF5 rg file: ", lfv[1])}
  make_h5se_rg(platform = accinfo[["platform"]], version = version, ts = ts, 
               dbn = file.path(comp.dpath, lfv[1])); return(NULL)
}

#' Get HDF5 database of methylated/unmethylated signals, from h5se rg data
#'
#' @param files.dpath Path to instance files directory.
#' @param comp.dname Name of compilations directory.
#' @param ngsm.block Number of GSM IDs per processed data block.
#' @return NULL, generates an h5 gm database file as side effect.
#' @export
get_h5db_gm <- function(files.dpath = "recount-methylation-files", 
                        comp.dname = "compilations", ngsm.block = 50){
  comp.dpath <- file.path(files.dpath, comp.dname)
  if(!dir.exists(comp.dpath)){
    stop("Error, didn't find compilations dir ", comp.dpath)}
  message("Handling metadata options...")
  md <- rmp_handle_metadata(); if(is.null(md)){
    stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); accinfo <- rmp_handle_platform()
  message("Checking for h5se RGChannelSet database file...")
  vform <- gsub("\\.", "-", version);lfv <- list.files(comp.dpath)
  lfv <- lfv[grepl(vform, lfv) & grepl(ts, lfv)]
  lfv <- lfv[grepl(".*h5se_rg.*", lfv)]
  if(length(lfv) == 0){
    stop("Couldn't find HDF5 rg database file at: ", comp.dpath, ".\n",
         "Try running rule `get_h5se_rg` first.")
  } else{message("Using HDF5 rg file: ", lfv[1])}
  dbn.dpath <- file.path(comp.dpath, lfv[1])
  message("Making HDF5 gm database file from h5se RGChannelSet...")
  make_h5_gm(dbn = dbn.dpath, version = version, ts = ts, 
             blocksize = ngsm.block, platform = platform);return(NULL)
}

#' Get HDF5-SummarizedExperiment of methylated/unmethylated signals, 
#' from h5se rg data
#'
#' @param files.dpath Path to instance files directory.
#' @param comp.dname Name of compilations directory.
#' @param ngsm.block Number of GSM IDs per processed data block.
#' @return NULL, generates an h5se gm database directory as side effect.
#' @export
get_h5se_gm <- function(files.dpath = "recount-methylation-files", 
                        comp.dname = "compilations"){
  comp.dpath <- file.path(files.dpath, comp.dname)
  if(!dir.exists(comp.dpath)){stop("Error, didn't find compilations dir ", 
                                   comp.dpath)}
  message("Handling metadata options...");md <- rmp_handle_metadata()
  if(is.null(md)){stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); accinfo <- rmp_handle_platform()
  message("Checking for valid HDF5 MethylSet ('gm') database file...")
  vform <- gsub("\\.", "-", version);lfv <- list.files(comp.dpath)
  lfv <- lfv[grepl(vform, lfv) & grepl(ts, lfv)] # filter on instance metadata
  lfv <- lfv[grepl(".*gm.*", lfv) & grepl(".*\\.h5$", lfv)];
  if(length(lfv) == 0){
    stop("Couldn't find HDF5 gm database file at: ", comp.dpath, ".\n",
         "Try running rule `get_h5_gm` first.")
  } else{message("Using HDF5 gm file: ", lfv[1])}
  dbn.dpath <- file.path(comp.dpath, lfv[1])
  message("Making h5se gm database file...")
  make_h5se_gm(dbn = dbn.dpath, version = version, ts = ts, platform=platform)
  return(NULL)
}

#' Get HDF5 database of DNAm fractions, from h5se rg data
#'
#' @param files.dpath Path to instance files directory.
#' @param comp.dname Name of compilations directory.
#' @param ngsm.block Number of GSM IDs per processed data block.
#' @return NULL, generates an h5 gr database file as side effect.
#' @export
get_h5db_gr <- function(files.dpath = "recount-methylation-files", 
                        comp.dname = "compilations", ngsm.block = 50){
  comp.dpath <- file.path(files.dpath, comp.dname)
  if(!dir.exists(comp.dpath)){
    stop("Error, didn't find compilations dir ", comp.dpath)}
  message("Handling metadata options...")
  md <- rmp_handle_metadata(); if(is.null(md)){
    stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); accinfo <- rmp_handle_platform()
  message("Checking for h5se RGChannelSet database file...")
  vform <- gsub("\\.", "-", version);lfv <- list.files(comp.dpath)
  lfv <- lfv[grepl(vform, lfv) & grepl(ts, lfv)]
  lfv <- lfv[grepl(".*h5se_rg.*", lfv)]
  if(length(lfv) == 0){
    stop("Couldn't find HDF5 rg database file at: ", comp.dpath, ".\n",
         "Try running rule `get_h5se_rg` first.")
  } else{message("Using HDF5 rg file: ", lfv[1])}
  dbn.dpath <- file.path(comp.dpath, lfv[1])
  message("Making HDF5 gr database file from h5se RGChannelSet...")
  make_h5_gr(dbn = dbn.dpath, version = version, ts = ts, 
             blocksize = ngsm.block, platform = platform);return(NULL)
}

#' Get HDF5-SummarizedExperiment of DNAm fractions, from h5se rg data
#'
#' @param files.dpath Path to instance files directory.
#' @param comp.dname Name of compilations directory.
#' @param ngsm.block Number of GSM IDs per processed data block.
#' @return NULL, generates an h5se gr database directory as side effect.
#' @export
get_h5se_gr <- function(files.dpath = "recount-methylation-files", 
                        comp.dname = "compilations"){
  comp.dpath <- file.path(files.dpath, comp.dname)
  if(!dir.exists(comp.dpath)){stop("Error, didn't find compilations dir ", 
                                   comp.dpath)}
  message("Handling metadata options...");md <- rmp_handle_metadata()
  if(is.null(md)){stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); accinfo <- rmp_handle_platform()
  message("Checking for valid HDF5 `gr` set database file...")
  vform <- gsub("\\.", "-", version);lfv <- list.files(comp.dpath)
  lfv <- lfv[grepl(vform, lfv) & grepl(ts, lfv)] # filter on instance metadata
  lfv <- lfv[grepl(".*gr.*", lfv) & grepl(".*\\.h5$", lfv)];
  if(length(lfv) == 0){
    stop("Couldn't find HDF5 gr database file at: ", comp.dpath, ".\n",
         "Try running rule `get_h5_gr` first.")
  } else{message("Using HDF5 gr file: ", lfv[1])}
  dbn.dpath<-file.path(comp.dpath,lfv[1]);message("Making h5se gr db file...")
  make_h5se_gr(dbn = dbn.dpath, version = version, ts = ts, platform=platform)
  return(NULL)
}

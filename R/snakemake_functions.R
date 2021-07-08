#!/usr/bin/env R

# Author: Sean Maden
# Functions for instance snakemake workflow.

#------------------
# instance metadata
#------------------

# library(recountmethylation.pipeline)

#' Get new instance metadata
#'
#' Generates the new metadata, including version and timestamp, for the
#' `recountmethylation` instance. This is called by rule `new_instance_md`.
#' 
#' @param instdir Path to directory to contain all instance metadata.
#' @param md.dname Name of dir to contain the instance metadata.
#' @param sett.path Path to the "settings.py" settings file containing the 
#' platform accession ID.
#' @return NULL, produces the instance metadata file as side effect.
#' @export
new_instance_md <- function(files.dname = "recount-methylation-files", 
                            md.dname = "metadata",
                            sett.path = file.path("recountmethylation_server", 
                                                  "src", "settings.py")){
  instdir <- file.path(files.dname, md.dname)
  if(!dir.exists(instdir)){message("Making instdir: ", instdir)
    dir.create(instdir)}
  message("Provide version:"); version <- readLines("stdin", n = 1)
  message("Using ", version, " as the instance version...")
  message("Using `", instdir, "` as the instance directory...")
  if(!dir.exists(instdir)){
    message("Making dir ", instdir, "...");dir.create(instdir)}
  md <- get_metadata("compilation_metadata", version = version)
  message("Detecting platform from `settings.py`...")
  setting.lines <- readLines(sett.path)
  platform.line <- setting.lines[grepl(".*platformid =.*", setting.lines)]
  platform.acc<-gsub(".* |'","",platform.line);md$accessionID<-platform.acc
  md$platform <- ifelse(platform.acc == "GPL13534", "hm450k",
                        ifelse(platform.acc == "GPL21145", "epic-hm850k",
                               ifelse(platform.acc=="GPL8490","hm27k","NA")))
  message("Saving new metadata...")
  md.fn <- paste0("metadata_v", gsub("\\.", "-", md$version), 
                  "_", md$timestamp, ".rda")
  md.path <- file.path(instdir, md.fn)
  message("Saving instance metadata to ", file.path(md.path))
  save(md, file = md.path); return(NULL)
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
  if(!dir.exists(comp.dpath)){
    message("Making new compilations dir: ",comp.dpath);dir.create(comp.dpath)}
  message("Handling metadata options...")
  md <- rmp_handle_metadata(); if(is.null(md)){stop("Couldn't get metadata...")}
  platform <- md[["platform"]]; message("Using platform ", platform, "...")
  version <- md[["version"]]; ts <- md[["timestamp"]]
  dtables_rg(platform=platform, version=version, ts=ts, destpath=comp.dpath)
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
  if(!dir.exists(comp.dpath)){stop("Error, didn't find compilations dir ", 
                                   comp.dpath)}
  message("Handling metadata options...")
  md <- rmp_handle_metadata(); 
  if(md == "NA"){stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); platform <- md[["platform"]]
  message("Using platform: ", platform)
  message("Checking for signal compilation tables...")
  vform <- gsub("\\.", "-", version)
  lfv <- list.files(comp.dpath);lfv <- lfv[grepl(vform, lfv) & grepl(ts, lfv)]
  lfv <- lfv[grepl("\\.compilation$", lfv)]
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
  message("Getting platform info...");
  platform <- md[["platform"]];message("Using platform: ", platform)
  message("Checking for valid HDF5 RGChannelSet file...")
  vform <- gsub("\\.", "-", version);lfv <- list.files(comp.dpath)
  lfv<-lfv[grepl(vform,lfv)&grepl(ts,lfv)];lfv<-lfv[grepl(".*h5_rg.*",lfv)]
  if(length(lfv) == 0){
    stop("Couldn't find HDF5 rg database file at: ", comp.dpath, ".\n",
         "Try running rule `get_h5db_rg` first.")
  } else{message("Using HDF5 rg file: ", lfv[1])}
  make_h5se_rg(newfnstem = "remethdb", platform = platform, 
               version = version, ts = ts, 
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
  message("Getting platform info...");platform <- md[["platform"]];
  message("Using platform: ", platform)
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
  message("Getting platform info..."); platform <- md[["platform"]];
  message("Using platform: ", platform)
  message("Checking for valid HDF5 MethylSet ('gm') database file...")
  vform <- gsub("\\.", "-", version);lfv <- list.files(comp.dpath)
  lfv<-lfv[grepl(vform,lfv)&grepl(ts,lfv)];lfv<-lfv[grepl(".*h5_gm.*",lfv)]
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
  message("Getting platform info...");
  platform <- md[["platform"]];message("Using platform: ", platform)
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
  message("Getting platform info...");platform <- md[["platform"]];
  message("Using platform: ", platform)
  message("Checking for valid HDF5 `gr` set database file...")
  vform <- gsub("\\.", "-", version);lfv <- list.files(comp.dpath)
  lfv <-lfv[grepl(vform,lfv)&grepl(ts,lfv)];lfv<-lfv[grepl(".*h5_gr.*",lfv)]
  if(length(lfv) == 0){
    stop("Couldn't find HDF5 gr database file at: ", comp.dpath, ".\n",
         "Try running rule `get_h5_gr` first.")
  } else{message("Using HDF5 gr file: ", lfv[1])}
  dbn.dpath<-file.path(comp.dpath,lfv[1]);message("Making h5se gr db file...")
  make_h5se_gr(dbn = dbn.dpath, version = version, ts = ts, platform=platform)
  return(NULL)
}

#------------------------
# handle metadata options
#------------------------


#' Get mapped metadata
#'
#' Snakemake function to get mapped metadata for a recountmethylation instance.
#' 
#' @param files.dname Files dir name for instance ("recount-methylation-files")
#' @param md.dname Metadata dir name ("metadata").
#' @return NULL, stores new mapped metadata.
#' @export
get_mdmap <- function(files.dname = "recount-methylation-files", 
                      md.dname = "metadata"){
  message("Handling metadata options...");md <- rmp_handle_metadata()
  if(is.null(md)){stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info...");platform <- md[["platform"]];
  message("Using platform: ", platform)
  message("Mining sample/GSM titles...");get_jsontitle(ts = ts)
  message("Getting study annotation tables from JSON files...")
  suppressMessages(get_atables(ts = ts))
  message("Performing metadata preprocessing...");md_preprocess(ts = ts)
  md.dpath <- file.path(files.dname, md.dname); md.lf <- list.files(md.dpath)
  mdpre.cond <- grepl(paste0(".*",ts,".*"), md.lf) & 
    grepl(paste0("md_preprocess.*"), md.lf);mdpre.fname <- md.lf[mdpre.cond][1]
  if(length(mdpre.fname) == 0){
    stop("Preprocessed metadata not found at ", md.dpath, "...")
  };mdpre <- get(load(file.path(md.dpath, mdpre.fname)))
  message("Performing metadata postprocessing (term mapping)...")
  suppressMessages(md_postprocess(ts = ts, mdpre = mdpre));
  message("Completed metadata mapping. Completed files are ",
          "located at: ", md.dpath);return(NULL)
}

#' Get DNAm-based metadata
#'
#' Gets DNAm-based metadata, including model-based predictions and quality 
#' metrics.
#'
#' @param files.dname Files dir name for instance ("recount-methylation-files")
#' @param md.dname Metadata dir name ("metadata"), contained at files.dname.
#' @param comp.dname Compilations dir name ("compilations"), contained at 
#' files.dname.
#' @return NULL, stores new DNAm-based metadata file at md.dpath.
#' @export
get_mddnam <- function(files.dname = "recount-methylation-files", 
                       md.dname = "metadata", comp.dname = "compilations"){
  message("Handling metadata options...");md <- rmp_handle_metadata()
  if(is.null(md)){stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); platform <- md[["platform"]];
  message("Using platform: ", platform)
  message("Detecting DNAm compilations...")
  comp.dpath <- file.path(files.dname, comp.dname)
  comp.lf <- list.files(comp.dpath)
  cond <- grepl(paste0(".*",ts,".*"), comp.lf) &
    grepl(paste0(".*",gsub("\\.", "-", version),".*"), comp.lf)
  comp.lf <- comp.lf[cond]
  h5se.rg.fn <- comp.lf[grepl(".*h5se_rg.*", comp.lf)][1]
  h5se.gr.fn <- comp.lf[grepl(".*h5se_gr.*", comp.lf)][1]
  h5se.gm.fn <- comp.lf[grepl(".*h5se_gm.*", comp.lf)][1]
  if(length(h5se.rg.fn[!is.na(h5se.rg.fn)]) == 0){
    stop("Couldn't find h5se rg dataset.")}
  if(length(h5se.gm.fn[!is.na(h5se.gm.fn)]) == 0){
    stop("Couldn't find h5se gm dataset.")}
  if(length(h5se.gr.fn[!is.na(h5se.gr.fn)]) == 0){
    stop("Couldn't find h5se gr dataset.")}
  message("Getting model predictions...")
  md_predictions(ts = ts, rgset.fname = h5se.rg.fn, grset.fname = h5se.gr.fn)
  message("Getting quality metrics...");get_qcmetrics(ts = ts, 
                                                      rgset.fname = h5se.rg.fn,
                                                      gmset.fname = h5se.gm.fn)
  message("Getting replicates info...");get_replicates(ts = ts, 
                                                       rg.fname = h5se.rg.fn)
  return(NULL)
}

#' Get available metadata as a single table
#'
#' Detects available metadata and aggregates tables into a single composite 
#' table.
#'
#' @param files.dname Files dir name for instance ("recount-methylation-files")
#' @param md.dname Metadata dir name ("metadata"), contained at files.dname.
#' @param comp.dname Compilations dir name ("compilations"), contained at 
#' files.dname.
#' @return NULL, stores new composite metadata table.
#' @export
#'
get_all_md <- function(files.dname = "recount-methylation-files", 
                       md.dname = "metadata", comp.dname = "compilations"){
  message("Handling metadata options...");md <- rmp_handle_metadata()
  if(is.null(md)){stop("Couldn't get metadata...")};ts <- md[["timestamp"]]
  message("Getting platform info...");platform <- md[["platform"]];
  message("Using platform: ", platform)
  message("Aggregating metadata...");md_agg(ts = ts, platform = platform)
  return(NULL)
}

#' Append metadata to compilation files
#' 
#' Append aggregate metadata to compilation files. This step follows 
#' formation of the compilation files (e.g. rule `get_rg_compilations`, etc.)
#' and metadata files (e.g. rule `do_mdmap`, etc.).
#'
#' @param files.dname Files dir name for instance ("recount-methylation-files")
#' @param md.dname Metadata dir name ("metadata"), contained at files.dname.
#' @param comp.dname Compilations dir name ("compilations"), contained at 
#' files.dname.
#' @return NULL, stores new composite metadata table.
#' @export
append_md <- function(files.dname = "recount-methylation-files", 
                      md.dname = "metadata", comp.dname = "compilations"){
  message("Handling metadata options...");md <- rmp_handle_metadata()
  if(is.null(md)){stop("Couldn't get metadata...")};ts <- md[["timestamp"]]
  message("Appending metadata...")
  append_md(ts = ts, files.dname = files.dname, md.dname = md.dname, 
            comp.dname = comp.dname); return(NULL)
}

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
new_instance_md <- function(instdir = file.path("recount-methylation-files",
                                                "metadata")){
  if(!dir.exists(instdir)){message("Making instdir: ", instdir)
    dir.create(instdir)}
  message("Provide version:"); version <- readLines("stdin", n = 1)
  message("Using ", version, " as the instance version...")
  message("Using `", instdir, "` as the instance directory...")
  if(!dir.exists(instdir)){
    message("Making dir ", instdir, "...");dir.create(instdir)}
  md <- get_metadata("compilation_metadata", version = version)
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

#------------------------
# handle metadata options
#------------------------


#' Get mapped metadata
#'
#' Snakemake function to get mapped metadata for a recountmethylation instance.
#' 
#' @param files.dname Files dir name for instance ("recount-methylation-files")
#' @param md.dname Metadata dir name ("metadata").
#' @return NULL, stored new mapped metadata.
#' @export
get_mdmap <- function(files.dname = "recount-methylation-files", 
                      md.dname = "metadata"){
  message("Handling metadata options...");md <- rmp_handle_metadata()
  if(is.null(md)){stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); accinfo <- rmp_handle_platform()
  platform <- accinfo[["platform_name"]];message("Using platform: ", platform)
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
#' @return NULL, stored new DNAm-based metadata file at md.dpath.
#' @export
get_mddnam <- function(files.dname = "recount-methylation-files", 
                       md.dname = "metadata", comp.dname = "compilations"){
  message("Handling metadata options...");md <- rmp_handle_metadata()
  if(is.null(md)){stop("Couldn't get metadata...")}
  version <- md[["version"]]; ts <- md[["timestamp"]]
  message("Getting platform info..."); accinfo <- rmp_handle_platform()
  platform <- accinfo[["platform_name"]];message("Using platform: ", platform)
  
  
  message("Detecting DNAm compilations...")
  comp.dpath <- file.path(files.dname, comp.dname)
  comp.lf <- list.files(comp.dpath)
  cond <- grepl(paste0(".*",ts,".*"), comp.lf) &
    grepl(paste0(".*",gsub("\\.", "-", version),".*"), comp.lf)
  comp.lf <- comp.lf[cond]
  h5se.rg.fn <- comp.lf[grepl(".*h5se_rg.*", comp.lf)][1]
  h5se.gr.fn <- comp.lf[grepl(".*h5se_gr.*", comp.lf)][1]
  h5se.gm.fn <- comp.lf[grepl(".*h5se_gm.*", comp.lf)][1]
  if(length(h5se.rg.fn) == 0){stop("Couldn't find h5se rg dataset.")}
  if(length(h5se.gm.fn) == 0){stop("Couldn't find h5se gm dataset.")}
  if(length(h5se.gr.fn) == 0){stop("Couldn't find h5se gr dataset.")}
  
  message("Getting model predictions...")
  md_predictions(ts = ts, rgset.fname = h5se.rg.fn, grset.fname = h5se.gr.fn)
  message("Getting quality metrics...");get_qcmetrics()
  message("Getting replicates info...");get_replicates()
  return(NULL)
}




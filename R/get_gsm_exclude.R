#!/usr/bin/env R

# Author: Sean Maden
# Get GSM IDs to exclude for instance.

#' Get GSM IDs to exclude for this instance
#'
#' Get the GSM IDs to exclude from target/download for this instance. The
#' IDs vector informs an updated edirect query filter object (e.g. file
#' named "gsequery_filt.*").
#' 
#' @param md Metadata table containing sample ID column (id.gsm.cname, 
#' e.g. "gsm"). If this is NULL, try loading an se object instead.
#' @param fn.str Name of sample IDs vector file to save.
#' @param fn.dpath Directory path to save samples ID vector.
#' @param se.path Path to a valid HDF5-SummarizedExperiment to check/load.
#' @param id.gsm.cname 
#' @param verbose
#' @return Vector of GSM IDs, with option to save.
#' @export
get_gsm_exclude <- function(md = NULL, se.path = NULL, fn.str = "gsm_exclude", 
                            fn.dpath = ".", id.gsm.cname = "gsm", verbose = TRUE){
  if(is.null(md) & is.null(se.path)){stop("Provide either md or se.path.")}
  if(!is.null(md)){
    if(!id.gsm.cname %in% colnames(md)){
      stop("Coudln't find sample ID column ",id.gsm.cname," in md...")}
  } else{
    if(!file.exists(rg.path)){
      stop("Path to se object ", se.path," doesn't exist.")}
    if(verbose){message("Loading the se object...")}
    se <- HDF5Array::loadHDF5SummarizedExperiment(se.path)
    md <- minfi::pData(se)
    if(!id.gsm.cname %in% colnames(md)){
      stop("Couldn't find sample ID column ",id.gsm.cname," in se metadata.")}}
  if(verbose){message("Writing sample ID vector...")}
  gsmv.path <- file.path(fn.dpath, fn.str)
  write(as.character(md[,id.gsm.cname]), sep = " ", file = gsmv.path)
  return(gsmv.path)}

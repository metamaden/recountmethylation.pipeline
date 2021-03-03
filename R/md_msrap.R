#!/usr/bin/env R

# Author: Sean Maden
# 
# Get the MetaSRA-pipeline output data as a flat table.
#

# get mapped terms from metasra-pipeline

#' Get the pipeline data files as a flat table
#'
#' Reads data in the sample/GSM-level metadata files produced by the pipeline 
#' and "map_msrap.py" script. Data is read and stored as a flat table ready to
#' be appended to other metadata for the instance.
#' 
#' @param ts NTP timestamp for the new metadata table (string).
#' @param files.dir Files dir for the instance ("recount-methylation-files").
#' @param md.dname Metadata dir name, located in the files.dir ("metadata").
#' @param msrap.dname Name of dir, located in the files.dir, containing 
#' pipeline outputs for samples/GSM IDs ("gsm_msrap_outfiles").
#' @param msrap.regex.str Regex pattern string for valid metadata pipeline 
#' output files ("^msrapout.*").
#' @param gsmid.fnindex File name index of the GSM ID for valid metadata 
#' pipeline output files (2).
#' @param gmap.fn File name stem of the new mapped metadata table 
#' ("md_msrapout").
#' @return NULL, outputs the mapped data table to the metadata directory.
#' @export
get_msrap <- function(ts, files.dir = "recount-methylation-files", md.dname = "metadata",
                      msrap.dname = "gsm_msrap_outfiles", msrap.regex.str = "^msrapout.*",
                      gsmid.fnindex = 2, gmap.fn = "md_msrapout"){
  rpath <- file.path(files.dir, msrap.dname)
  wpath <- file.path(files.dir, md.dname, paste(gmap.fn, ts, collapse = "."))
  rl <- list.files(rpath); rl <- rl[grepl(msrap.regex.str, rl)]
  gmap <- matrix(nrow = 0, ncol = 2);lfl.filt <- list.files(rpath)
  for(f in rl){
    gsmi <- unlist(strsplit(f,"\\."))[gsmid.fnindex]
    jfi <- jsonlite::fromJSON(txt = file.path(dpath, f), flatten=T)
    xi <- jfi[,!colnames(jfi) %in% 
                c("mapped ontology terms", "real-value properties")]
    pi <- paste0("'",names(xi),"':'",as.character(xi[1,]),"'",collapse=";")
    if("real-value properties" %in% names(jfi)){
      if(length(unlist(jfi[,"real-value properties"]))>0){
        pi <- paste0(pi,";'real-value properties':",
                     as.character(unlist(jfi[,"real-value properties"])),
                     collapse="")}}
    mot <- paste0(as.character(unlist(jfi[,"mapped ontology terms"])), collapse=";")
    jfi.mapped.flat<-paste0(pi,";",mot)
    gmapi<-matrix(c(gsmi,jfi.mapped.flat),nrow=1);gmap <- rbind(gmap, gmapi)}
  gmap <- as.data.frame(gmap, stringsAsFactors = FALSE)
  colnames(gmap) <- c("gsm", "msrapout");save(gmap, file = wpath);return(NULL)}

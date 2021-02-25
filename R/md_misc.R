#!/usr/bin/env Rscript

# Author: Sean Maden
#
# Preprocess available sample metadata
#

#------------------
# md misc functions
#------------------

# get annotation tables

#' Generate study annotation tables
#'
#' Generate a list of study annotation tables from filtered sample JSON files.
#' Sample/GSM JSON files should be filtered to remove study-level information,
#' retaining just sample-level metadata.
#' 
#' @param ts Timestamp for the preprocessed metadata table to output 
#' (integer or character).
#' @param atable.fn Name of the annotation tables file to output 
#' ("geo_gse-atables_list)
#' @param eq.dname Name of directory, in files.dname, containing the equery 
#' files produced by the server ("equery").
#' @param md.dname Name of directory, in files.dname, containing the instance
#' metadata files ("metadata).
#' @param json.filt.dname Name of directory, in files.dname, containing the
#' filtered GSM JSON files (e.g. containing only sample metadata, 
#' "gsm_json_filt")
#' @param files.dname Main recountmethylation instance files directory 
#' ("recount-methylation-files").
#' @param verbose Whether to show status messages (TRUE).
#' @return NULL, produces the tables list file.
#' @seealso md_preprocess(); md_postprocess()
#' @export
get_atables <- function(ts, atable.fn = "geo_gse-atables_list",
                        eq.dname = "equery", md.dname = "metadata",
                        json.filt.dname = "gsm_json_filt",
                        files.dname = "recount-methylation-files",
                        verbose = TRUE){
  md.dpath <- file.path(files.dname, md.dname)
  if(!dir.exists(md.dpath)){
    message("Making md dir ",md.dpath); dir.create(md.dpath)}
  message("Getting paths...")
  eqdpath <- file.path(files.dname, eq.dname) # equery dir path
  eqd.ts <- max(as.numeric(gsub(".*\\.", "", list.files(eqdpath))))
  eq.fpath <- recount.synth::get_eqfpath(ts = eqd.ts, eqdpath = eqdpath)
  if(length(eq.fpath) == 0){
    stop("Couldn't find an equery file with ts ", eqd.ts," at ", eqdpath, ".")}
  eq.fn <- gsub(".*/", "", eq.fpath) # equery file name
  message("Getting GSM IDs by GSE IDs, from equery file...")
  x <- scan(eq.fpath, what="", sep="\n");gsel <- list()
  for(i in 1:length(x)){
    ssi <- unlist(strsplit(x[[i]]," "));gsel[[ssi[1]]] <- ssi[2:length(ssi)]}
  message("Getting tables list from JSON files...")
  json.filt.dpath <- file.path(files.dname, json.filt.dname)
  lf.all <- list.files(json.filt.dpath); 
  lf.all <- lf.all[grepl(".*\\.json\\.filt$", lf.all)]; tgse.list <- list()
  for(g in seq(length(gsel))){
    lgse <- list();gseid <- as.character(names(gsel)[g]);ggse <- gsel[[g]]
    message("Parsing GSM JSON metadata...")
    for(j in seq(length(ggse))){
      ffilt <- grepl(ggse[j], lf.all) & grepl("\\.filt", lf.all)
      lff <- lf.all[ffilt][1]; cond <- length(lff)>0 & !(lff=="NA"|is.na(lff))
      if(cond){
        fnj <- file.path(json.filt.dpath, lff)
        lgse[[ggse[j]]] <- rjson::fromJSON(paste(readLines(fnj), collapse=""))}}
    if(length(lgse) > 0){
      message("Making table columns from unique GSM keys..."); tcols <- c()
      for(l in 1:length(lgse)){
        gsmid <- names(lgse)[l]; gsmval <- unlist(lgse[[l]])
        for(k in 1:length(gsmval)){
          if(grepl(":",gsmval[k])){
            kk <- as.character(gsub(":.*","",gsmval[k]))
            if(!kk=="" & !kk %in% tcols){tcols <- c(tcols, kk)}}}; message(l)
      };tcols <- c("gsm", "gse", tcols)
      tgse <- matrix(nrow=0,ncol=length(tcols));colnames(tgse) <- tcols
      message("Coercing GSM data to GSE/study-specific table")
      for(l in 1:length(lgse)){
        gsmid <- names(lgse)[l]; gsmval <- unlist(lgse[[l]])
        gvk <- c(gsmid, gseid);message("Looping on GSM values...")
        for(c in 3:ncol(tgse)){gvv <- "NA"; tc <- colnames(tgse)[c]
        for(i in 1:length(gsmval)){
          gsmdati <- as.character(gsub(".*:","",gsmval[i]))
          gsmlabi <- as.character(gsub(":.*","",gsmval[i]))
          if(tc %in% gsmlabi){gvv <- as.character(gsmdati)}
        }; gvk <- c(gvk,gvv)
        }; tgse <- rbind(tgse, gvk); message(l)
      }; tgse.list[[gseid]] <- tgse
    }; message("gse:", g)
  };atable.fpath <- file.path(md.dpath, paste0(atable.fn,"_", ts, ".rda"))
  message("Saving annotation tables data to ", atable.fpath, "...")
  save(tgse.list, file = atable.fpath); return(NULL)
}

# get gsm titles from json files

#' Get sample titles
#' 
#' Mines the sample titles from the instance sample/GSM JSON files. These
#' are included in the mapped metadata tables and used in term mappings.
#' 
#' @param ts Timestamp for the preprocessed metadata table to output 
#' (integer or character).
#' @param json.dname Name of directory, in files.dname, containing the instance 
#' filtered sample/GSM JSON files ("gsm_json_filt").
#' @param jsonext.regex Regex pattern to identify filtered JSON files 
#' (".*\\.json.filt$").
#' @param json.titlesdf.fn 
#' @param md.dname Name of directory, in files.dname, containing the instance 
#' metadata files ("metadata).
#' @param files.dname Main recountmethylation instance files directory 
#' ("recount-methylation-files").
#' @return NULL, produces a file containing the sample titles.
#' @seealso md_preprocess(); md_postprocess(); get_atables()
#' @export
get_jsontitle <- function(ts, json.dname = "gsm_json_filt",
                          jsonext.regex = ".*\\.json.filt$",
                          json.titlesdf.fn = "gsm_jsontitledf",
                          md.dname = "metadata", 
                          files.dname = "recount-methylation-files"){
  json.dpath <- file.path(files.dname, json.dname)
  md.dpath <- file.path(files.dname, md.dname)
  message("Looking for JSON files at path ", json.dpath, "...")
  lf<-list.files(json.dpath);ffilt<-grepl(jsonext.regex,lf);lff<-lf[ffilt]
  message("Getting JSON titles from ", length(lff), " files...")
  gsmtitledf <- matrix(nrow = 0, ncol = 2)
  for(jf in lff){
    jsym <- unlist(strsplit(jf, "\\.")); gsm.id <- jsym[2]
    json.lines <- readLines(file.path(json.dpath, jf))
    json.convert <- rjson::fromJSON(paste(json.lines,collapse="")) 
    st.catch <- as.character(unlist(json.convert)["!Sample_title"])
    stm <- matrix(c(gsm.id, st.catch), nrow=1)
    gsmtitledf <- rbind(gsmtitledf, stm)
  };colnames(gsmtitledf) <- c("gsm", "gsm_title")
  json.titlesdf.fpath <- file.path(md.dpath, 
                                   paste0(json.titlesdf.fn, "_", ts, ".rda"))
  message("Saving JSON titles to file ", json.titlesdf.fpath, "...")
  save(gsmtitledf, file = json.titlesdf.fpath); return(NULL)
}



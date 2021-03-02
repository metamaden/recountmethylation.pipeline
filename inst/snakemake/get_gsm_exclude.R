#!/usr/bin/env R

# Author: Sean Maden
# Get GSM IDs to exclude for instance.

get_gsm_exclude_path <- function(){
  # check for default file at current dir
  lf <- list.files(); lf <- lf[grepl("gsm_exclude", lf)][0]
  if(length(lf) > 0){
    message("Detected existing sample ID vector: ", lf)
    gsmv.path <- lf} else{gsmv.path <- get_gsm_exclude()}
  return(gsmv.path)}

get_gsm_exclude_path()



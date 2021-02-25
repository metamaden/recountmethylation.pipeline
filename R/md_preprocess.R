#!/usr/bin/env Rscript

# Author: Sean Maden
#
# Preprocess available sample metadata
#

#------------------------
# md preprocess functions
#------------------------

# define strings to search for variables

#' Terms vector for variable tissue
#' 
#' Terms to seed regex pattern matching for md_preprocess().
#' @seealso md_preprocess(); md_postprocess()
#' @return Terms vector
#' @export
mdpre_vars_tissue <- function(){
  cname.tissue <- c("tissue", "tissu", "sample type", "sample_type", 
                    "sample.type", "cell line", "cell_line", "cell.line",
                    "cell type", "cell_type", "cell.type", "tissue type", 
                    "tissue_type", "tissue.type", "tissue region", 
                    "tissue_region", "tissue.region", "histology", 
                    "region", "brnum", "location") 
  return(cname.tissue)
}

#' Terms vector for variable disease
#' 
#' Terms to seed regex pattern matching for md_preprocess().
#' @seealso md_preprocess(); md_postprocess()
#' @return Terms vector
#' @export
mdpre_vars_disease <- function(){
  cname.disease_state <- c("disease state", "disease_state", "disease.state",
                           "subject status", "subject_status", 
                           "subject.status", "sample group", "sample_group", 
                           "sample.group", "diagnosis", "group", "condition",
                           "disease", "subgroup", "status") 
  return(cname.disease_state)
}

#' Terms vector for variable info
#' 
#' Terms to seed regex pattern matching for md_preprocess().
#' @seealso md_preprocess(); md_postprocess()
#' @return Terms vector
#' @export
mdpre_vars_info <- function(){
  cname.info <- c("state", "passages", "race", "individual",
                  "race", "stage", "condition", "risk", "twinid", "smok",
                  "cigarette", "pack year", "pack_year", "pack.year", 
                  "pack yr", "pack_yr", "pack.yr", "drinks", "drink year", 
                  "drink_year", "drink.year", "drink yr", "drink_yr", 
                  "drink.yr", "drug use", "drug_use", "drug.use", 
                  "alcohol", "treatment", "material", "run", "batch", 
                  "plate", "developmental stage", "developmental_stage", 
                  "developmental.stage", "source", "storage", "zygosity", 
                  "family", "weight", "bmi", "drug", "intervention")
  return(cname.info)
}

#' Terms vector for variable sex
#' 
#' Terms to seed regex pattern matching for md_preprocess().
#' @seealso md_preprocess(); md_postprocess()
#' @return Terms vector
#' @export
mdpre_vars_sex <- function(){cname.sex <- c("sex", "gender");return(cname.sex)}

#' Terms vector for variable age
#' 
#' Terms to seed regex pattern matching for md_preprocess().
#' @seealso md_preprocess(); md_postprocess()
#' @return Terms vector
#' @export
mdpre_vars_age <- function(){
  cname.age <- c("age", "passage", "age (years)", "age_(years)", "age.(years)")
  return(cname.age)
}

# get list of search strings

#' Get the regex patterns for variable mappings
#'
#' Gets the regex patterns from vectors of seed terms for each variable mapped 
#' by md_preprocess().
#' @seealso md_preprocess(); md_postprocess()
#' @return Returns list containing regex patterns for each mapped variable.
#' @export
mdpre_vars <- function(){
  cl <- lapply(list(mdpre_vars_tissue(), mdpre_vars_disease(), 
                    mdpre_vars_sex(), mdpre_vars_age(), 
                    mdpre_vars_info()), get_pstr);
  names(cl) <- c("sample_type", "disease_state", "sex", "age", "info")
  return(cl)
}

#' Preprocess sample metadata
#'
#' Preprocess sample metadata by coercing JSON-formatted metadata files
#' into a flat matrix, with regex pattern matching to identify and 
#' categorize various variable types. Additional data such as the sample
#' titles contained in the file declared by the titles.fn argument.
#'
#' @param ts Timestamp for the preprocessed metadata table to output 
#' (integer or character).
#' @param mdpre.fn Name of preprocessed metadata table file to output 
#' ("md_preprocess").
#' @param md.dname Name of directory, in files.dname, containing the instance 
#' metadata files ("metadata).
#' @param titlesfn.str Name of sample/GSM titles file ("gsm_jsontitledf).
#' @param atablefn.str Name of study annotation tables file 
#' ("geo_gse-atables_list")
#' @param files.dname Main recountmethylation instance files directory 
#' ("recount-methylation-files").
#' @param verbose Whether to show status messages (TRUE).
#' @return NULL, produces the mdpre preprocessed metadata table.
#' @seealso md_postprocess()
#' @export
md_preprocess <- function(ts, mdpre.fn = "md_preprocess", 
                          md.dname = "metadata", 
                          titlesfn.str = "gsm_jsontitledf",
                          atablefn.str = "geo_gse-atables_list",
                          files.dname = "recount-methylation-files",
                          verbose = TRUE){
  if(verbose){message("Loading GSM data files...")}
  md.dpath <- file.path(files.dname, md.dname)
  md.lf<-list.files(md.dpath);md.lf<-md.lf[grepl(paste0(".*",ts,".*"),md.lf)]
  at.fname <- md.lf[grepl(paste0(atablefn.str, ".*"), md.lf)][1]
  titles.fname <- md.lf[grepl(paste0(titlesfn.str, ".*"), md.lf)][1]
  if(!file.exists(file.path(md.dpath, at.fname))){
    stop("Study anno. tables file not found at ", at.fpath)}
  if(!file.exists(file.path(md.dpath, titles.fname))){
    stop("GSM titles file not found at ",titles.fpath)}
  tls <- get(load(file.path(md.dpath, at.fname)))
  tdf <- get(load(file.path(md.dpath, titles.fname)))
  message("Processing GSE flat files...")
  colv <- c("gsm","gse","sample_type","disease_state","sex","age")
  gat.all <- matrix(nrow = 0, ncol = length(colv));colnames(gat.all) <- colv
  cl <- mdpre_vars(); message("Appending tables data...")
  for(gseid in names(tgse.list)){
    gsedat <- tgse.list[[gseid]]
    gati <- data.frame(gsm = gsedat[,1], gse=gsedat[,2],stringsAsFactors=FALSE)
    for(cn in names(cl)){
      cnvar <- rep("NA", nrow(gsedat))
      cn.dat <- rep("NA", nrow(gsedat));gsedat.rep <- gsedat;cn.cv <- cl[[cn]]
      cname.filt <- grepl(cn.cv, colnames(gsedat.rep))
      if(length(which(cname.filt)) > 0){
        if(cn %in% c("age", "info")){ # parse age mapping logic, excluding fp's
          cmv <- colnames(gsedat.rep)[cname.filt]
          if(verbose){message("Removing info matches (e.g. stage columns...)")}
          if(cn == "age"){
            cmv.info <- grepl(cl[["info"]], colnames(gsedat.rep))
            cmv <- colnames(gsedat.rep)[cname.filt & !cmv.info]}
          if(verbose){message("Appending colnames to row entries...")}
          for(colnamei in cmv){
            gsedat.rep[,colnamei]<-paste0(colnamei,":",gsedat.rep[,colnamei])}
        };gf <- gsedat.rep[, cname.filt, drop = FALSE]
        cnvar <- as.character(apply(gf,1,paste,collapse = ";"))
      };gati[,ncol(gati) + 1] <- cnvar; colnames(gati)[ncol(gati)] <- cn
    };gat.all <- rbind(gat.all, gati);message("Finished study: ", gseid)
  };gat.all <- gat.all[!duplicated(gat.all[,1]),]
  message("Appending GSM titles...")
  d1 <- gat.all; d2 <- gsmtitledf;gsm.all <- unique(c(d1[,1], d2[,1]))
  gsm1 <- gsm.all[!gsm.all %in% d1[,1]];gsm2 <- gsm.all[!gsm.all %in% d2[,1]]
  if(length(gsm1) > 0){
    nav <- rep(rep("NA", length(gsm1)), ncol(d1) - 1)
    mna <- matrix(c(gsm1, nav), nrow = length(gsm1), ncol = ncol(d1))
    colnames(mna) <- colnames(d1); d1 <- rbind(d1, mna)}
  if(length(gsm2) > 0){
    nav <- rep(rep("NA", length(gsm2)), ncol(d2) - 1)
    mna <- matrix(c(gsm2, nav), nrow = length(gsm2), ncol = ncol(d2))
    d2 <- rbind(d2, mna)}
  if(nrow(d2) > nrow(d1)){d2 <- d2[d2[,1] %in% d1[,1] & !duplicated(d2[,1]),]}
  match.gsm1 <- match(as.character(d1[,1]), as.character(d2[,1]))
  order.gsm1 <- order(match.gsm1);d1 <- d1[order.gsm1,]
  match.gsm2 <- match(as.character(d2[,1]), as.character(d1[,1]))
  order.gsm2 <- order(match.gsm2);d2 <- d2[order.gsm2,]
  cond <- identical(as.character(d2[,1]), as.character(d1[,1]))
  if(cond){
    d1 <- as.data.frame(d1, stringsAsFactors = FALSE)
    d1$gsm_title <- as.character(d2[,2])}
  mdpre.fpath <- file.path(md.dpath, paste0(mdpre.fn, "_", ts, ".rda"))
  message("Saving mdpre at ", mdpre.fpath, "...")
  mdpre <- d1; save(mdpre, file = mdpre.fpath); return(NULL)
}

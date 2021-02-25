#!/usr/bin/env Rscript

# Author: Sean Maden
#
# Preprocess available sample metadata
#

require(data.table); require(rjson)
library(wateRmelon) # for age predictions
library(minfi) # for cell type and sex predictions
library(HDF5Array)
library(minfi)
library(ewastools)


#--------------
# dnam metadata
#--------------
# model-based predictoins

#' Make model-based metadata predictions 
#'
#' Make model-based metadata predictions using DNAm assays. Predictions for 
#' age, sex, and blood cell fractions are produced. Sex and blood cell 
#' predictions use raw/unnormalized red/green signal from h5se rg data 
#' (rgset.fname argument), while age predictions use noob-normalized 
#' Beta-values from h5se gr data (grset.fname argument).
#' 
#' @param ts Timestamp for the preprocessed metadata table to output 
#' (integer or character).
#' @param rgset.fname Name of the compilation file containing red signals 
#' extracted from IDATs (character).
#' @param grset.fname Name of the compilation file containing green signals 
#' extracted from IDATs (character).
#' @param mdmod.fname Name of the table of model-based predictions produced 
#' ("mdmod_dnam-predictions").
#' @param files.dname Main recountmethylation instance files directory 
#' ("recount-methylation-files").
#' @param md.dname Name of directory, in files.dname, containing the instance
#' metadata files ("metadata).
#' @param nsamp.block Number of samples per data block processed (integer, 50).
#' @param rgset.path Path 
#' @param comp.dname Name of directory, in files.dname, containing the
#' compilation files with red and green signals for prediction calculations 
#' ("compilations")
#' @param verbose Whether to show status messages (TRUE).
#' @return NULL, produces table of model-based metadata predictions.
#' @seealso get_qcmetrics(); md_preprocess(); md_postprocess()
#' @export
md_predictions <- function(ts, rgset.fname, grset.fname,
                           mdmod.fname = "mdmod_dnam-predictions", 
                           files.dname = "recount-methylation-files",
                           md.dname = "metadata", nsamp.block = 50,
                           comp.dname = "compilations", verbose = TRUE){
  md.dpath <- file.path(files.dname, md.dname)
  rgset.path <- file.path(files.dname, comp.dname, rgset.fname)
  grset.path <- file.path(files.dname, comp.dname, grset.fname)
  if(!file.exists(rgset.path)){stop("Couldn't locate red signals compilation ",
                                    "file at:\n",rgset.path)}
  if(!file.exists(grset.path)){stop("Couldn't locate red signals compilation ",
                                    "file at:\n",grset.path)}
  if(verbose){message("Loading h5se datasets...")}
  rgset <- HDF5Array::loadHDF5SummarizedExperiment(rgset.path)
  grset <- HDF5Array::loadHDF5SummarizedExperiment(grset.path)
  if(verbose){message("Getting model predictions for blocks of ",
                      nsamp.block," samples...")}
  if(!ncol(rgset) == ncol(grset)){
    if(verbose){message("Warning, ncol not identical for rgset and grset. ",
                        "Choosing lowest column count for sample indexing...")}
    col.tot <- ifelse(ncol(rgset) > ncol(grset), ncol(grset), ncol(rgset))
  } else{ncol.tot <- ncol(rgset)};blocks <- getblocks(ncol.tot, nsamp.block)
  if(length(blocks[[length(blocks)]]) == 1 & length(blocks) == 1){
    stop("Detected 1 block of 1 sample. Sample blocks should have > 1 sample.")
  } else if(length(blocks[[length(blocks)]]) == 1 & length(blocks) > 1){
    if(verbose){message("Aggregating final block with 1 sample...")}
    blocks[length(blocks)-1] <- list(c(blocks[[length(blocks)-1]], 
                                           blocks[[length(blocks)]]))
    blocks <- blocks[1:(length(blocks)-1)]
  } else{if(verbose){message("No singleton blocks found. Continuing...")}}
  if(verbose){message("Getting predictions...")};mdmod<-matrix(nrow=0,ncol=8)
  for(b in blocks){
    rgf <- rgset[, unlist(b)] # sex and cell preds use raw red/grn signals
    rgf.matrix <- minfi::RGChannelSet(Green = as.matrix(minfi::getGreen(rgf)),
                                      Red = as.matrix(minfi::getRed(rgf)),
                                      annotation=BiocGenerics::annotation(rgf))
    celltypepred <- minfi::estimateCellCounts(rgf.matrix)
    msf <- minfi::mapToGenome(minfi::preprocessRaw(rgf))
    sexpred <- minfi::getSex(msf);grf <- grset[,colnames(rgf)]
    predage<-wateRmelon::agep(minfi::getBeta(grf)) # age pred from noob'd bvals
    mdf<-cbind(predage,cbind(sexpred[,3],celltypepred));mdmod<-rbind(mdmod,mdf)
    if(verbose){message("Completed index ", max(unlist(b)))}}
  if(verbose){message("Formatting mdmod...")}
  mdmod <- as.data.frame(mdmod, stringsAsFactors = FALSE)
  colnames(mdmod) <- c("predage", "predsex", 
                       paste0("predcell.", colnames(mdmod)[3:8]))
  mdmod$gsm<-gsub("\\..*", "",rownames(mdmod))
  mdmod<-mdmod[!duplicated(mdmod$gsm),]
  mdmod.fpath <- file.path(md.dpath, paste0(mdmod.fname,"_", ts, ".rda"))
  if(verbose){message("Saving mdmod to ", mdmod.fpath, "...")}
  save(mdmod, file = mdmod.fpath); return(NULL)
}

# get qc metrics from red/grn signals

#' Get quality metrics from DNAm assays
#'
#' Calculates quality metrics from DNAm assays contained in the compilation 
#' files. These include signals for 17 BeadArray controls, 2 controls for the 
#' methylated and unmethylated signals, and predicted genotypes by sample
#' using methods in the ewastools package.
#' 
#' @param ts Timestamp for the preprocessed metadata table to output 
#' (integer or character).
#' @param mdqc.fname Name of the quality metrics table output ("mdqc")
#' @param athresh Similarity threshold (percent similarity) for the ewastools 
#' predicted genotype (decimal, 0.1).
#' @param nsamp.block Samples per data block processed (integer, 50).
#' @param md.dname Name of directory, in files.dname, containing the instance
#' metadata files ("metadata).
#' @param files.dname Main recountmethylation instance files directory 
#' ("recount-methylation-files").
#' @param verbose Whether to show status messages (TRUE). 
#' @return NULL, produces a table of quality metrics.
#' @seealso md_predictions(); md_preprocess(); md_postprocess()
#' @export
get_qcmetrics <- function(ts, rgset.fname, gmset.fname, mdqc.fname = "mdqc",
                          athresh = 0.1, nsamp.block = 50, 
                          md.dname = "metadata", comp.dname = "compilations",
                          files.dname = "recount-methylation-files",
                          verbose = TRUE){
  if(verbose){message("Loading h5se data...")}
  rgset.path <- file.path(files.dname, comp.dname, rgset.fname)
  gmset.path <- file.path(files.dname, comp.dname, gmset.fname)
  if(!file.exists(rgset.path)){stop("Couldn't find h5se rg at: ",rgset.path)}
  if(!file.exists(rgset.path)){stop("Couldn't find h5se gm at: ",gmset.path)}
  rgset <- HDF5Array::loadHDF5SummarizedExperiment(rgset.path)
  gmset <- HDF5Array::loadHDF5SummarizedExperiment(gmset.path)
  mdqc.fpath <- file.path(md.dpath, paste0(mdqc.fname, "_", ts, ".rda"))
  message("Working on file: ", mdqc.fpath, "...")
  md.dpath <- file.path(files.dname, md.dname)
  blocks <- getblocks(slength = ncol(rgset), bsize = nsamp.block)
  ms <- matrix(nrow = 0, ncol = 20)
  cdf <- as.data.frame(getProbeInfo(rgset, type = "Control"))
  if(verbose){message("Getting model predictions for blocks of ",
                      nsamp.block," samples...")}
  if(!ncol(rgset) == ncol(gmset)){
    if(verbose){message("Warning, ncol not identical for rgset and grset. ",
                        "Choosing lowest column count for sample indexing...")}
    col.tot <- ifelse(ncol(rgset) > ncol(gmset), ncol(gmset), ncol(rgset))
  } else{ncol.tot <- ncol(rgset)};blocks <- getblocks(ncol.tot, nsamp.block)
  if(length(blocks[[length(blocks)]]) == 1 & length(blocks) == 1){
    stop("Detected 1 block of 1 sample. Sample blocks should have > 1 sample.")
  } else if(length(blocks[[length(blocks)]]) == 1 & length(blocks) > 1){
    if(verbose){message("Aggregating final block with 1 sample...")}
    blocks[length(blocks)-1] <- list(c(blocks[[length(blocks)-1]], 
                                       blocks[[length(blocks)]]))
    blocks <- blocks[1:(length(blocks)-1)]
  } else{if(verbose){message("No singleton blocks found. Continuing...")}}
  dft <- df_bathresh() # BeadArray control thresholds
  for(bi in seq(length(blocks))){
    b <- blocks[[bi]];rgf <- rgset[, b];gmf <- gmset[, b]
    colnames(rgf) <- gsub("\\..*", "", colnames(rgf))
    basignals <- suppressMessages(get_bactrl(rgset = rgf))
    bathresh <- get_bathresh(basignals)
    bathresh <- bathresh[,c(2:ncol(bathresh))]
    colnames(bathresh) <- paste0(colnames(bathresh),".thresh",as.numeric(dft[1,]))
    mqc <- cbind(basignals, bathresh)
    if(verbose){message("Getting log2 median meth/unmeth from gm set...")}
    ms <- as.matrix(getMeth(gmf)); us <- as.matrix(getUnmeth(gmf))
    meth.l2med <- apply(ms, 2, function(x){log2(median(x))})
    unmeth.l2med <- apply(us, 2, function(x){log2(median(x))})
    dfs <- data.frame(meth.l2med = meth.l2med, unmeth.l2med = unmeth.l2med,
                      stringsAsFactors = FALSE)
    mi <- cbind(mqc, dfs);ms <- rbind(ms, mi);
    if(verbose){message("Finished block ", bi)}}
  save(mdqc, file = mdqc.fpath);return(NULL)
}

#' Get likely sample replicates
#'
#' Use inferred genotypes to predict sample replicates by study. This function
#' uses the genotype prediction methods in the ewastools package to predict
#' likely sample replicates for each included study.
#' 
#' @param ts Timestamp for the metadata output (integer or character).
#' @param rgset An RGChannelSet object.
#' @param athresh Similarity threshold for replicate identification (0.1).
#' @param md Preprocessed metadata, rsheet, or other table containing columns
#' "gsm" and "gse". If NULL, attempt to find mdpre at md.dname (NULL).
#' @param md.dname Name of metadata dir ("metadata"), contained at files.dname.
#' @param comp.dname Name of compilations dir ("compilations"), contained at 
#' files.dname.
#' @param files.dname Name of files dir for the instance 
#' ("recount-methylation-files").
#' @param verbose Whether to show status messages (TRUE).
#' @return NULL, stores a table of replicate information.
#' @export
get_replicates <- function(ts, rg.fname, athresh = 0.1, md = NULL,
                           md.dname = "metadata", comp.dname = "compilations",
                           files.dname = "recount-methylation-files",
                           verbose = TRUE){
  if(verbose){message("Loading h5se rg set...")}
  comp.dpath <- file.path(files.dname, comp.dname);
  comp.lf <- list.files(comp.dpath)
  rgset.fname <- comp.lf[grepl(paste0(".*",ts,".*"), comp.lf) & 
                       grepl(".*h5se_rg.*", comp.lf)][1]
  rgset.fpath <- file.path(comp.dpath, rgset.fname)
  if(!file.exists(rgset.fpath)){stop("Couldn't load h5se rg at ",rgset.fpath)}
  rgset <- HDF5Array::loadHDF5SummarizedExperiment(rgset.fpath)
  if(is.null(md)){
    md.dpath <- file.path(files.dname, md.dname);md.ls <- list.files(md.dpath)
    md.fname <- md.ls[grepl(paste0(".*",ts,".*"), md.ls) & 
                           grepl(paste0(".*md_preprocess.*"), md.ls)][1]
    if(length(md.fname)==0){stop("Didn't find md file at: ",md.dpath)}
  };md <- get(load(file.path(md.dpath, md.fname)))
  if(!"gsm" %in% colnames(md) | !"gse" %in% colnames(md)){
    stop("Metadata table requies 'gsm' and 'gse' ID columns")}
  snp1.info <- minfi::getProbeInfo(rgset, type = "SnpI")
  snp2.info <- minfi::getProbeInfo(rgset, type = "SnpII")
  snp.addr <- c(snp1.info$AddressA, snp1.info$AddressB, snp2.info$AddressA)
  beta.snp <- minfi::getSnpBeta(rgset)
  colnames(beta.snp) <- gsub("\\..*", "", colnames(beta.snp))
  ama<-matrix(nrow=0,ncol=3);colnames(ama)<-c("gsm","gseid","cgsnp.gsm.geno")
  for(id in unique(md$gse)){
    gsmv <- md[md$gse == id,]$gsm
    sbi <- beta.snp[, colnames(beta.snp) %in% gsmv, drop = FALSE]
    if(ncol(sbi) > 1){
      sbi <- as.matrix(sbi); class(sbi) <- "numeric"
      sgenoi <- ewastools::call_genotypes(sbi, learn = FALSE)
      sagree <- ewastools::check_snp_agreement(sgenoi, 
                                               donor_ids = colnames(sbi),
                                               sample_ids = colnames(sbi))
      for(sai in sagree){
        am <- data.frame(gsm = unique(c(sai$sample1, sai$sample2)), 
                         stringsAsFactors = FALSE)
        am$gse<-id;am$gsm.geno<-"NA"
        for(gsm in am$gsm){
          sai.cond <- (sai$sample1==gsm | sai$sample2==gsm) & 
            sai$agreement > athresh
          sai.gsm <- sai[sai.cond,]
          sai.id <- unique(c(sai.gsm$sample1, sai.gsm$sample2))
          am[am$gsm == gsm,]$gsm.geno <- paste(sai.id, collapse = ";")
        };ama <- rbind(ama, am)
      }
    } else{if(verbose){message("Not enough samples with SNP data for study ",
                               id,". Skipping")}}
    if(verbose){message("Finished with study ", id)}}
  ama$num.shared <- unlist(lapply(ama$gsm.geno, function(x){
    length(unique(unlist(strsplit(x, ";"))))}));return(ama)
}

#!/usr/bin/env R

# Functions for the Recount Methylation Pipeline and `recountmethylation` R package.

#----------
# Utilities
#----------

# Generic utilities

#' Get list of index blocks
#'
#' Get list of index blocks allowing for remainders.
#'
#' @param slength Total length of index vector
#' @param bsize Size of index blocks along length
#' @return List of index blocks of min length `slength`/`bsize`
#' @export
getblocks <- function(slength, bsize){
  sc <- 1; ec <- sc + bsize - 1
  iv <- list()
  nblocks <- slength %/% bsize
  for(b in 1:nblocks){
    iv[[b]] <- seq(sc, ec, 1)
    sc <- ec + 1; ec <- ec + bsize
  }
  # add final indices
  if(nblocks < (slength/bsize)){
    iv[[length(iv) + 1]] <- seq(sc, slength, 1)
  }
  return(iv)
}

# Make the base h5 db from data tables

#' Append metadata to HDF5 or SE object
#' 
#' Passes version info and timestamp from Python to object metadata
#' @param version Numeric version to be passed, should conform to ##.##.## nomenclature
#' @return Metadata content for the object
#' @export 


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
h5_newtables <- function(dbn = "remethdb2.h5", dsn.nb = "noobbeta",
                         dsn.meth = "methylated_signal",
                         dsn.unmeth = "unmethylated_signal",
                         dsn.red = "redsignal", dsn.grn = "greensignal",
                         verbose = TRUE, ngsm.block = 50,
                         ncol.chunk = 5000){
  # Generate noobbeta and meth/unmeth signal tables
  # get dimensions from red and grn signal data
  rs.rn <- rhdf5::h5read(dbn, paste0(dsn.red, ".rownames"))
  rs.cn <- rhdf5::h5read(dbn, paste0(dsn.red, ".colnames"))
  gs.rn <- rhdf5::h5read(dbn, paste0(dsn.grn, ".rownames"))
  gs.cn <- rhdf5::h5read(dbn, paste0(dsn.grn, ".colnames"))
  sbv <- getblocks(length(rs.rn), ngsm.chunk)
  # get new cg  dimensions
  anno.name = "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  man = eval(parse(text = paste(anno.name, "Manifest", sep = "::")))
  ncg = nrow(man)
  # new h5 data params
  newdims <- c(length(rs.rn), ncg)
  chunkvars <- c(5, ncol.chunk)
  # make new tables
  rhdf5::h5createDataset(dbn, "unmethylated_signal", dims = newdims, 
                         maxdims = c(H5Sunlimited(), H5Sunlimited()), 
                         storage.mode = "double", level = 5, chunk = chunkvars)
  rhdf5::h5createDataset(dbn, "methylated_signal", dims = newdims, 
                         maxdims = c(H5Sunlimited(), H5Sunlimited()), 
                         storage.mode = "double", level = 5, chunk = chunkvars)
  rhdf5::h5createDataset(dbn, "noobbeta", dims = newdims, 
                         maxdims = c(H5Sunlimited(), H5Sunlimited()), 
                         storage.mode = "double", level = 5, chunk = chunkvars)
  # append new data
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
    se.nb <- preprocessNoob(se.rgi)
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
      message("iter: ", i," of ", 
              length(sbv),", time: ", 
              Sys.time() - t1)
    }
  }
  if(verbose){message("Completed addition of new tables!")}
  return(TRUE)
}

# H5 database to SEH5 dir objects with DelayedArray

#' Use DelayedArray function to store H5 SummarizedExperiment directory
#'
#' Use DelayedArray function to store H5 SummarizedExperiment directory from h5 db. Handles 3 classes of Summarized Experiment objects.
#'
#' @param newfn Filename of new H5 SE directory.
#' @param dsn.data1 Name of first dataset in h5 db (required).
#' @param dsn.data2 Name of second dataset in h5 db (required if se either rg or gm).
#' @param mdpost Sample metadata file (columns are fields, rows are samples).
#' @param dsn.rn Name of object in h5 db with the rownames (GSM or sample basenames).
#' @param se Type of SummarizedExperiment to parse (either rg for RGChannelSet, gr for GenomicRatioSet, or gm for GenomicMethylSet).
#' @param verbose Whether to show verbose status messages.
#' @param dsn.cn Name of h5 db object with column names (cgids or assays, optional).
#' @param dbn Name or path to h5 db. 
#' @param replace.opt Whether to replace/overwrite any existing H5 SE directories of the same name as `newfn`.
#' @return List of index blocks of min length `slength`/`bsize`
#' @export
make.h5se <- function(newfn, dsn.data1, dsn.data2 = NULL, 
                      mdpost, dsn.rn = "redsignal.rownames", 
                      se = c("rg", "gr", "gm"), verbose = TRUE,
                      dsn.cn = FALSE, dbn = "remethdb2.h5", replace.opt = TRUE){
  # Sets up SummarizedExperiment creation from h5 file
  # Stores SE object to SE H5 file with DelayedArray processing
  
  # anno for se sets
  anno = c("IlluminaHumanMethylation450k", "ilmn12.hg19")
  names(anno) = c("array", "annotation")
  # get granges object
  if(se %in% c("gr", "gm")){
    require(minfiData)
    data(MsetEx)
    mrset <- minfi::mapToGenome(MsetEx)
    grcg <- GenomicRanges::granges(mrset)
  }
  # load data table
  if(verbose){message("Getting dsn.data1...")}
  ldat <- list()
  nb <- HDF5Array::HDF5Array(dbn, dsn.data1)
  rn <- rhdf5::h5read(dbn, dsn.rn)
  rownames(nb) <- as.character(rn)
  nb <- t(nb)
  ldat[[1]] <- nb
  # sanity checks and parse data2
  if(!is.null(dsn.data2)){
    if(verbose){message("Getting dsn.data2...")}
    nb <- HDF5Array::HDF5Array(dbn, dsn.data2)
    rn <- rhdf5::h5read(dbn, dsn.rn)
    rownames(nb) <- as.character(rn)
    nb <- t(nb)
    if(!identical(nrow(nb), nrow(ldat[[1]])) | 
       !identical(ncol(nb), ncol(ldat[[1]]))){
      stop("Matrix dsn.data2 not similar dim to dsn.data1!")
    }
    ldat[[2]] <- nb
  } else if (se %in% c("rg", "gm")){
    stop("Must provide dsn.data2 for se as rg or gm!")
  }
  # get probe ids or addresses
  if(!dsn.cn){
    man = eval(parse(text = paste("IlluminaHumanMethylation450kanno.ilmn12.hg19", 
                                  "Manifest", sep = "::")))
    cgrn <- rownames(man)
    ldat <- lapply(ldat, function(x){
      rownames(x) <- cgrn
      return(x)
    })
    if(!se == "rg"){
      ordergr <- order(match(names(grcg), 
                             rownames(ldat[[1]])))
      grcg <- grcg[ordergr]
      if(!identical(rownames(ldat[[1]]), names(grcg)) |
         !identical(rownames(ldat[[2]]), names(grcg))){
        stop("Couldn't match grcg names to rows in at least one data matrix.")
      }
    }
  } else{
    # append cgids (stored cnames) to rows/assays
    cn <- rhdf5::h5read(dbn, dsn.rn)
    rownames(nb) <- as.character(cn)
  }
  # make the new se set
  if(verbose){message("Making the new se object...")}
  if(se == "rg"){
    if(verbose){message("Making RGChannelSet...")}
    gri <- minfi::RGChannelSet(Red = ldat[[1]], 
                               Green = ldat[[2]], 
                               anno = anno)
  } else if (se == "gr"){
    if(verbose){message("Making GenomicRatioSet...")}
    gri <- minfi::GenomicRatioSet(gr = grcg,
                                  Beta = ldat[[1]],
                                  anno = anno)
  } else{
    if(verbose){message("Making Making GenomicMethylSet...")}
    gri <- minfi::GenomicMethylSet(gr = grcg,
                                   Meth = ldat[[1]], 
                                   Unmeth = ldat[[2]],
                                   anno = anno)
  }
  # append pheno
  gri.pheno <- se.addpheno(mdp = mdpost, se = gri)
  # save seh5 object and finally run process
  t1 <- Sys.time()
  if(verbose){message("Saving seh5 object...")}
  HDF5Array::saveHDF5SummarizedExperiment(gri.pheno, 
                                          dir = "rmseh5_grnoob", 
                                          replace = replace.opt)
  if(verbose){
    message("Save complete! Process finished at time ", 
            Sys.time() - t1)
  }
  return(TRUE)
}

#' Append phenotype data to a SummarizedExperiment object
#'
#' Append phenotype data to a SummarizedExperiment object
#'
#' @param mdp Metadata object (columns are fields, rows are GSMs/samples).
#' @param se SummarizedExperiment object.
#' @return SummarizedExperiment object with appended phenotype data.
#' @export
se.addpheno <- function(mdp, se){
  # Adds pheno data to a SummarizedExperiment objects
  mdp <- mdp[mdp$basename %in% colnames(se),]
  bnv <- colnames(se)
  gsmv <- gsub("\\..*", "", bnv)
  gsmfilt <- !gsmv %in% mdp$gsm
  gsm.ov <- gsmv[gsmfilt]
  bnv.ov <- bnv[gsmfilt]
  # add na values
  nacol <- rep("NA", length(gsm.ov))
  mdp.ov <- matrix(c(gsm.ov,
                     rep(nacol, 6),
                     bnv.ov,
                     rep(nacol, 11)), 
                   nrow = length(gsm.ov))
  colnames(mdp.ov) <- colnames(mdp)
  mdp <- rbind(mdp, mdp.ov)
  mdp <- mdp[order(match(mdp$gsm, gsmv)),]
  identical(mdp$gsm, gsmv)
  # add valid basenames
  mdp$basename <- colnames(se)
  rownames(mdp) <- colnames(se)
  identical(rownames(mdp), colnames(se))
  # append pheno to se object
  minfi::pData(se) <- DataFrame(mdp)
  return(se)
}
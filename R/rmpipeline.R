#!/usr/bin/env R

# Functions for the Recount Methylation Pipeline and `recountmethylation` R package.

#----------
# Utilities
#----------

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

#' Append metadata to HDF5 or SE object
#' 
#' Passes version info and timestamp from Python to object metadata
#' @param title Object title
#' @param version Numeric version to be passed, should conform to ##.##.## nomenclature
#' @param pname Name of pipeline package (default: "rmpipeline")
#' @param sname Name of Python script (default: "get_timestamp.y")
#' @return Metadata content for the object
#' @export 
get.metadata <- function(title, version, pname = "rmpipeline", 
                         sname = "get_timestamp.py"){
  mdl <- list(title = title, version = version)
  # get timestamp from package python script
  path <- paste(system.file(package = pname), 
                sname, sep="/")
  command <- paste("python", path, sname, 
                   "to", sep = " ")
  mdl[["timestamp"]] <- system(command, intern = T)
  return(mdl)
}

#----------------------------
# Make data tables from IDATs
#----------------------------

#' Make data tables from IDATs
#' 
#' Make data tables from IDATs, including red and green channel signals and noob-normalized Beta-values. 
#' This function should be run from the "recount-methylation" base directory.
#' 
#' @param version Version of the run for data table filenames.
#' @param verbose Whether to return verbose notifications.
#' @param gsmint Number of GSMs to process at a time, typically runs best near 50 samples.
#' @param fnstem Filename stem for data tables.
#' @param sepval Separator symbol for data being written.
#' @param getnb Whether to get noob-normalized Beta-vlaues (default: FALSE).
#' @return Creates new HDF5 database from DNAm signal data.
#' @export

h5db.fromsignal <- function(version, verbose = TRUE, gsmint = 60,
                            fnstem = "mdat.compilation", sepval = " ", getnb = FALSE,
                            idats.path = paste("recount-methylation-files", "idats", sep = "/"),
                            dest.path = paste("recount-methylation-analysis", 
                                              "files", "mdata", "compilations", sep = "/")){
  
  # get run metadata
  runmd <- get.metadata(version, "notitle")
  nts = runmd[["timestamp"]]
  
  # get valid gsms idats dir
  idats.lf = list.files(idats.path)
  which.valid1 = grepl("\\.idat$", substr(idats.lf, nchar(idats.lf) - 4, nchar(idats.lf)))
  which.valid2 = grepl(".*hlink.*", idats.lf)
  idats.valid = idats.lf[which.valid1 & which.valid2]
  gsmu = gsub("\\..*", "", idats.valid)
  gsmu = gsmu[grepl("^GSM.*", gsmu)]
  gsmu = unique(gsmu)
  
  # determine target idats
  # note: add max timestamps filter
  gpath = c()
  for(i in 1:length(gsmu)){
    cont.cond  <- FALSE
    g = gsmu[i]
    gfilt = grepl(paste0(".*", g, ".*"), idats.valid)
    ig = idats.valid[gfilt]
    igv.red = ig[grepl(".*_Red.*", ig)][1]
    # check for matched green channel file
    while(length(ig.red) > 0 & !cont.cond){
      ig.red <- igv.red[1]
      fchan.filt = grepl(gsub("_Red.*", "", ig.red), ig) & grepl(".*_Grn.*", ig)
      ig.grn = ig[fchan.filt][1]
      cont.cond = length(ig.red) == 1 & length(ig.grn) == 1
      igv.red[!igv.red == ig.red]
    }
    if(cont.cond){
      gpath = c(gpath, gsub("_Red.*", "", ig.red))
    }
    if(verbose){message(
      "Finished validating idat fn's for sample ", i)
    }
  }
  gsmii = seq(1, length(gsmu), gsmint)
  
  # make new data paths
  reds.fn <- paste(paste("redsignal", nts, version, sep = "_"), 
                   fnstem, sep = ".")
  grns.fn <- paste(paste("greensignal", nts, version, sep = "_"), 
                   fnstem, sep = ".")
  reds.path = paste(destpath, reds.fn, sep = "/")
  grns.path = paste(destpath, grns.fn, sep = "/")
  if(getnb){
    nb.fn <- paste(paste("noobbeta", nts, version, sep = "_"), 
                   fnstem, sep = ".")
    nb.path = paste(destpath, nb.fn, sep = "/")
  }
  
  # instantiate new empty data tables
  cn = c("gsmi")
  rgi = minfi::read.metharray(c(paste(idats.path, gpath[1:2], sep = "/")))
  nbi = minfi::preprocessNoob(rgi)
  rgcni = colnames(t(getRed(rgi)))
  grcni = colnames(t(getGreen(rgi)))
  rgcni == grcni
  rgcn = matrix(c(cn, colnames(t(getRed(rgi)))), nrow = 1)
  nbcn = matrix(c(cn, colnames(t(getBeta(nbi)))), nrow = 1)
  fwrite(rgcn, reds.path, sep = sepval, append = FALSE)
  fwrite(rgcn, grns.path, sep = sepval, append = FALSE)
  if(getnb){
    fwrite(nbcn, nb.path, sep = sepval, append = FALSE)
  }
  
  # append new methdata
  tt = Sys.time()
  for(i in 1:length(gsmii)){
    gi = gsmii[i]
    # read in new data
    pathl = paste(idats.path, gpath[gi:(gi + gsmint - 1)], sep = "/")
    rgi = minfi::read.metharray(c(pathl))
    
    # get data matrices
    redi = matrix(c(colnames(rgi), t(getRed(rgi))), ncol = nrow(rgi) + 1)
    grni = matrix(c(colnames(rgi), t(getGreen(rgi))), ncol = nrow(rgi) + 1)
    
    # append new data
    data.table::fwrite(redi, reds.path, sep = sepval, append = TRUE)
    data.table::fwrite(grni, grns.path, sep = sepval, append = TRUE)
    
    # parse noob-normalized data option
    if(getnb){
      gsi = preprocessNoob(rgi)
      nbi = matrix(c(colnames(gsi), t(getBeta(rgi))), ncol = nrow(gsi) + 1)
      data.table::fwrite(nbi, grns.path, sep = sepval, append = TRUE)
    }
    
    if(verbose){
      message("Finished gsmi ", gi, " to ", gi + gsmint - 1, 
              ", time : ", Sys.time() - tt)
    }
  }
  return(NULL)
}

#------------------------------------
# Make and populate the HDF5 database
#------------------------------------

#' Append metadata to HDF5 or SE object
#' 
#' Add signal data (red and green channel) to the HDF5 database.
#' @param dbn Name of the HDF5 database targeted.
#' @param fnl List of signal tables containing data to be added.
#' @param dsnl List of data set names in HDF5 database to be populated. 
#' @param rmax Total rows to append to data sets, reflecting total samples.
#' @param cmax Total columns to append to data sets, reflecting total assays or probes.
#' @param verbose Whether to print verbose progress messages.
#' @param nr.inc Number of samples to append at a time (default: 10).
#' @return Populates the HDF5 database
#' @export 
h5.addtables = function(dbn, fnl, dsnl, rmax, cmax, 
                    verbose = TRUE, nr.inc = 10){
  for(di in 1:length(dsnl)){
    fnread = fnl[di]; dsn = dsnl[di]
    h5createDataset(dbn, dsn, dims = c(rmax, cmax), 
                    maxdims = c(H5Sunlimited(), H5Sunlimited()), 
                    storage.mode = "double", level = 5, chunk = c(1, 5))
    rn = cn = c()
    con <- file(fnread, "r")
    cn = unlist(strsplit(readLines(con, n = 1), " "))
    cn = cn[2:length(cn)] # filt first value 
    cn = gsub("\n", "",gsub('\"', '', cn[1:cmax])) # grab the max indexed value
    nri = seq(1, rmax, nr.inc)
    tt <- Sys.time()
    j = 1 # j tracks last line written to hdf5d
    for(i in nri){
      dati = unlist(strsplit(readLines(con, n = nr.inc), " "))
      wdi = which(grepl(".*GSM.*", dati))
      dff = matrix(nrow = 0, ncol = cmax)
      ngsm = gsub("\n", "", gsub('\"', '', gsub("\\..*", "", dati[wdi]))) # all read gsm ids
      wgsm = c() # new gsm ids to write
      for(wi in 1:length(wdi)){
        # filter redundant gsm ids
        if(!ngsm[wi] %in% rn){
          wadd = wdi[wi] + 1
          dff = rbind(dff, matrix(dati[wadd:(wadd + cmax - 1)], nrow = 1))
          wgsm = c(wgsm, ngsm[wi])
        }
      }
      rn = c(rn, wgsm) # gsm ids for new data
      class(dff) = "numeric"
      h5write(dff, file = dbn, name = dsn, 
              index = list(j:(j + nrow(dff) - 1), 1:cmax))
      h5closeAll()
      j = j + nrow(dff) # start index for next dff
      message("For ds ", dsn,", finished reading index ", i, " to ", i + (nr.inc - 1), 
              ", time; ", Sys.time() - tt)
    }
    message("Adding row and column names for ds ", dsn)
    cnn = paste0(dsn, ".colnames"); rnn = paste0(dsn, ".rownames");
    h5createDataset(dbn, cnn, dims = length(cn), maxdims = c(H5Sunlimited()), 
                    storage.mode = "character", level = 5, # compression level, 1-9
                    chunk = c(20), size = 256) # chunk dims
    message("Added colnames...")
    h5createDataset(dbn, rnn, dims = length(rn), maxdims = c(H5Sunlimited()), 
                    storage.mode = "character", level = 5, # compression level, 1-9
                    chunk = c(20), size = 256) # chunk dims
    message("Added rownames...")
    h5write(cn, file = dbn, name = cnn, index = list(1:length(cn)))
    h5write(rn, file = dbn, name = rnn, index = list(1:length(rn)))
    h5closeAll()
    message("Completed writing data, rownames, colnames for ds, ", dsn)
  }
  message("Completed hdf5 write for all ds's in list. Returning")
}

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
h5.newtables <- function(dbn, dsn.nb = "noobbeta",
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

#' Make and populate a new HDF5 database
#' 
#' Add signal data (red and green channel) to the HDF5 database.
#' @param dbn Name of the HDF5 database targeted.
#' @param fnl List of signal tables containing data to be added.
#' @param dsnl List of data set names in HDF5 database to be populated. 
#' @param rmax Total rows to append to data sets, reflecting total samples.
#' @param cmax Total columns to append to data sets, reflecting total assays or probes.
#' @param newtables Whether to also add new data tables (noob-norm. Beta-values, meth. and unmeth. signal).
#' @return Populates the HDF5 database
#' @export 
make.h5db <- function(dbn, fnl, rmax = 35500; cmax = 622399, newtables = TRUE,
                 dsnl = c("redsignal", "greensignal", "noobbeta")){
  try(rhdf5::h5createFile(dbn))
  # remove old data if present
  if(verbose){message("Removing any existing old data.")}
  for(d in dsnl){
    rhdf5::h5delete(dbn, d)
    rhdf5::h5delete(dbn, paste0(d, ".colnames"))
    rhdf5::h5delete(dbn, paste0(d, ".rownames"))
  }
  if(verbose){message("Adding and populating data tables to HDF5 database")}
  h5ds.add(fnl = fnl, dsnl = dsnl, rmax = rmax, cmax = cmax, nr.inc = nr.inc)
  if(verbose){message("Finished adding red and green channel data to h5 databse.")}
  
  if(newtables){
    if(verbose){message("Adding new data tables.")}
    h5.newtables(dbn)
  }
  
  if(verbose){message("Finished all processes. Returning.")}
  return(NULL)
}

#-----------------------
# Make the SE-H5 objects
#-----------------------

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
    metadata(gri) <- get_metadata(paste0("Recount Methylation ", class(gri)[1]),
                                  version)
  } else if (se == "gr"){
    if(verbose){message("Making GenomicRatioSet...")}
    gri <- minfi::GenomicRatioSet(gr = grcg,
                                  Beta = ldat[[1]],
                                  anno = anno)
    metadata(gri) <- get_metadata(paste0("Recount Methylation ", class(gri)[1]),
                                  version)
  } else{
    if(verbose){message("Making GenomicMethylSet...")}
    gri <- minfi::GenomicMethylSet(gr = grcg,
                                   Meth = ldat[[1]], 
                                   Unmeth = ldat[[2]],
                                   anno = anno)
    metadata(gri) <- get_metadata(paste0("Recount Methylation ", class(gri)[1]),
                                  version)
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
    message("Save complete. Process duration was ", 
            Sys.time() - t1)
  }
  return(TRUE)
}

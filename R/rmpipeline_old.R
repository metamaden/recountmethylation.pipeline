#--------------
# OLD FUNCTIONS
#--------------
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
dtables_fromsignal <- function(version, timestamp = NULL, verbose = TRUE, gsmint = 60,
                               fnstem = "mdat.compilation", sepval = " ", getnb = FALSE,
                               idatspath = paste("recount-methylation-files", "idats", sep = "/"),
                               destpath = paste("recount-methylation-analysis",
                                                "files", "mdata", "compilations", sep = "/")){
  # get run metadata
  if(is.null(timestamp)){
    runmd <- get_metadata(version, "notitle")
    nts = runmd[["timestamp"]]
  } else{
    nts <- timestamp
  }
  # get valid gsms idats dir
  if(verbose){message("Getting valid GSM IDs from IDAT filenames...")}
  idats.lf = list.files(idatspath)
  which.valid1 = grepl("\\.idat$", substr(idats.lf,
                                          nchar(idats.lf) - 4,
                                          nchar(idats.lf))) # idat ext
  which.valid2 = grepl(".*hlink.*", idats.lf) # hlink
  which.valid3 = length(gsub("\\..*", "",
                             gsub(".*\\.", "", idats.lf))) > 0 # timestamp
  idats.valid = idats.lf[which.valid1 &
                           which.valid2 &
                           which.valid3]
  gsmu = gsub("\\..*", "", idats.valid)
  gsmu = gsmu[grepl("^GSM.*", gsmu)]
  gsmu = unique(gsmu)
  # red/grn channel files
  if(verbose){message("Finding paired red and grn channel files...")}
  gstr <- gsub(".*hlink\\.", "", gsub("(_Red.idat$|_Grn.idat$)", "", idats.valid))
  rsub <- idats.valid[grepl(".*_Red.idat$", idats.valid)]
  gsub <- idats.valid[grepl(".*_Grn.idat$", idats.valid)]
  # intersect strings between ".*hlink\\." and "(_Red.idat|_Grn.idat)"
  gstr.rg <- intersect(gsub(".*hlink\\.", "", gsub("_Red.idat$", "", rsub)),
                       gsub(".*hlink\\.", "", gsub("_Grn.idat$", "", gsub)))
  idats.validf <- idats.valid[gstr %in% gstr.rg]
  gpath <- unique(gsub("(_Red.idat|_Grn.idat)", "", idats.validf))
  # filter on timestamps
  if(verbose){message("Filtering files on latest timestamps...")}
  gpath.ts <- c()
  gpath.gid <- c()
  gidv <- gsub("\\..*", "", gpath)
  tidv <- gsub(".*\\.", "", gsub("\\.hlink.*", "", gpath))
  # iterate on timestamps
  tsu <- unique(tidv)
  tsu <- tsu[rev(order(tsu))] # sort so that latest comes first
  for(t in tsu){
    gpathf <- gpath[tidv == t & !gidv %in% gpath.gid]
    # rm redundant gsm ids in this iter
    gpathf.gid <- gsub("\\..*", "", gpathf)
    gpathf <- gpathf[!duplicated(gpathf.gid)]
    # append new fn
    gpath.ts <- c(gpath.ts, gpathf)
    # remake gid list
    gpath.gid <- gsub("\\..*", "", gpath.ts)
  }
  gpath <- gpath.ts
  gsmu <- gsub("\\..*", "", gpath)
  if(verbose){message("Finished filtering filenames. ",
                      "Found ", length(gsmu),
                      " valid file pairs.")
  }
  # get GSM ID indices by interval, as list
  gsmii = getblocks(length(gsmu), gsmint)
  # make new data paths
  if(verbose){message("Making and instantiating new data table files...")}
  reds.fn <- paste(paste("redsignal", nts,
                         gsub("\\.", "-", version),
                         sep = "_"),
                   fnstem, sep = ".")
  grns.fn <- paste(paste("greensignal", nts,
                         gsub("\\.", "-", version),
                         sep = "_"),
                   fnstem, sep = ".")
  reds.path = paste(destpath, reds.fn, sep = "/")
  grns.path = paste(destpath, grns.fn, sep = "/")
  if(getnb){
    nb.fn <- paste(paste("noobbeta", nts,
                         gsub("\\.", "-", version),
                         sep = "_"),
                   fnstem, sep = ".")
    nb.path = paste(destpath, nb.fn, sep = "/")
  }
  # instantiate new empty data tables with probes as colnames
  cn = c("gsmi")
  rgi = minfi::read.metharray(c(paste(idatspath, gpath[1:2], sep = "/")))
  rgcni = colnames(t(minfi::getRed(rgi))) # grcni = colnames(t(getGreen(rgi)))
  # rgcni == grcni
  rgcn = matrix(c(cn, rgcni), nrow = 1)
  data.table::fwrite(rgcn, reds.path, sep = sepval, append = FALSE, col.names = F)
  data.table::fwrite(rgcn, grns.path, sep = sepval, append = FALSE, col.names = F)
  if(getnb){
    if(verbose){message("Instantiating noob beta table...")}
    nbi = minfi::preprocessNoob(rgi)
    nbcn = matrix(c(cn, colnames(t(getBeta(nbi)))), nrow = 1)
    data.table::fwrite(nbcn, nb.path, sep = sepval, append = FALSE, col.names = F)
  }
  # append new methdata
  tt = Sys.time()
  for(i in 1:length(gsmii)){
    gi = gsmii[[i]]
    # read in new data
    pathl = paste(idatspath, gpath[gi], sep = "/")
    rgi = minfi::read.metharray(c(pathl))
    
    # get data matrices
    redi = matrix(c(colnames(rgi), t(minfi::getRed(rgi))), ncol = nrow(rgi) + 1)
    grni = matrix(c(colnames(rgi), t(minfi::getGreen(rgi))), ncol = nrow(rgi) + 1)
    
    # append new data
    data.table::fwrite(redi, reds.path, sep = sepval, append = TRUE)
    data.table::fwrite(grni, grns.path, sep = sepval, append = TRUE)
    
    # parse noob-normalized data option
    if(getnb){
      gsi = minfi::preprocessNoob(rgi)
      nbi = matrix(c(colnames(gsi), t(minfi::getBeta(rgi))), ncol = nrow(gsi) + 1)
      data.table::fwrite(nbi, nb.path, sep = sepval, append = TRUE)
    }
    
    if(verbose){
      message("Finished gsmi ", gi[1],
              " to ", gi[length(gi)],
              ", time elapsed: ", Sys.time() - tt)
    }
  }
  return(NULL)
}

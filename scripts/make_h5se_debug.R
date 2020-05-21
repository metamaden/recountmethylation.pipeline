require(rmpipeline)
require(rhdf5)

# args
version <- "0.0.1"
ts <- 1590090412
dbn <- "remethdb_1590090412_0-0-1.h5"
fnstem <- "remethdb_h5se-rg"
dsnv <- c("redsignal", "greensignal")
#dsn.data1 = "redsignal"
#dsn.data2 = "greensignal"
se = "rg"
verbose = TRUE
dsn.rnv = c(paste0(dsnv[1], ".rownames"),
           paste0(dsnv[2], ".rownames"))
dsn.cnv = c(paste0(dsnv[1], ".colnames"),
           paste0(dsnv[2], ".colnames"))
addpheno = TRUE
phenopath = NULL
dsn.md <- "mdpost"
semd = list("title" = "Recount Methylation H5-SE Object",
            "version" = version,
            "timestamp" = ts,
            "preprocessing" = "raw")

# test function
newfn <- paste(newfnstem, gsub("\\.", "-", version), ts, sep = "_")
if(verbose){message("Setting annotation info...")}
anno = c("IlluminaHumanMethylation450k", "ilmn12.hg19")
names(anno) = c("array", "annotation")
if(verbose){message("Getting dsn.data1...")}
ldat <- list()
for(i in 1:length(dsnv)){
  nb <- HDF5Array::HDF5Array(dbn, dsnv[i])
  rn <- rhdf5::h5read(dbn, dsn.rnv[i])
  cn <- rhdf5::h5read(dbn, dsn.cnv[i])
  nb <- nb[c(1:length(rn)),c(1:length(cn))]
  rownames(nb) <- as.character(rn)
  colnames(nb) <- as.character(cn)
  nb <- t(nb)
  ldat[[dsnv[i]]] <- nb
}
# make the new H5-SE set(s)
if(verbose){message("Making RGChannelSet...")}
gri <- minfi::RGChannelSet(Red = ldat[[1]],
                           Green = ldat[[2]],
                           anno = anno)
S4Vectors::metadata(gri) <- semd
# parse md file options
if(addmd){
  if(verbose){message("Adding samples metadata...")}
  if(is.null(mdpath)){
    if(dsn.md %in% h5ls(dbn)$name &
       paste0(dsn.md, ".colnames") %in% h5ls(dbn)$name){
      if(verbose){message("Adding metadata from dbn...")}
      mdp <- recountmethylation::data_mdpost(dbn, dsn.md)
      gri <- se_addpheno(pdat = mdp, se = gri)
    }
    message("Couldn't add metadata. Set phenopath or check dsn.md in file dbn.")
  } else{
    if(verbose){message("Adding metadata from phenopath...")}
    gri <- se_addpheno(mdpath, se = gri)
  }
}
# start to instantiate the H5SE object
t1 <- Sys.time()
if(verbose){message("Starting process to make new file ", newfn, "...")}
HDF5Array::saveHDF5SummarizedExperiment(gri, dir = newfn, replace = replace.opt)
if(verbose){message("Save complete, time elapsed:", Sys.time() - t1)}
return(NULL)





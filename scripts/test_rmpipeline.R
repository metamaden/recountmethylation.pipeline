#!/usr/bin/env R

# Run the Recount Methylation Pipeline to generate data objects for the R package.

require(rmpipeline)

# make data tables from IDAT files
# rmpipeline::h5db.fromsignal(version = commandArgs(T)[1])

vers <- "00.00.01"
rmd <- get_metadata("newrun", vers)
ts <- rmd[["timestamp"]]

dtables_fromsignal(version = "0.0.1", timestamp = ts,
                   idatspath = "idats", destpath = "compilations")

{
  version = "0.0.1"
  timestamp = "123"
  idatspath = "idats"
  destpath = ""
  verbose = TRUE
  gsmint = 60

  if(is.null(timestamp)){
    runmd <- get_metadata(version, "notitle")
    nts = runmd[["timestamp"]]
  } else{
    nts <- timestamp
  }

  # get valid GSM IDs after checking IDAT fn strings
  if(verbose){message("Getting valid GSM IDs from IDAT filenames...")}
  idats.lf = list.files(idatspath)
  which.valid1 = grepl("\\.idat$", substr(idats.lf,
                                          nchar(idats.lf) - 4,
                                          nchar(idats.lf))) # idat ext
  which.valid2 = grepl(".*hlink.*", idats.lf) # hlink
  which.valid3 = length(gsub("\\..*", "", gsub(".*\\.", "", idats.lf))) > 0 # timestamp
  idats.valid = idats.lf[which.valid1 & which.valid2]
  gsmu = gsub("\\..*", "", idats.valid)
  gsmu = gsmu[grepl("^GSM.*", gsmu)]
  gsmu = unique(gsmu)

  # determine target idats
  # note: add max timestamps filter
  if(verbose){message("Determining valid IDAT file targets from GSM ID list...")}
  gpath = c()
  for(i in 1:length(gsmu)){
    cont.cond  <- FALSE
    g = gsmu[i]
    gfilt = grepl(paste0(".*", g, ".*"), idats.valid)
    ig = idats.valid[gfilt]
    # order files by their timestamp, newest to oldest
    ig <- ig[order(gsub("\\..*", "", gsub(".*\\.", "", ig)))]
    igv.red = ig[grepl(".*_Red.*", ig)][1]
    # check for matched green channel file
    ig.red <- igv.red[1]
    while(length(ig.red) > 0 & !cont.cond){
      fchan.filt = grepl(gsub("_Red.*", "", ig.red), ig) & grepl(".*_Grn.*", ig)
      ig.grn = ig[fchan.filt][1]
      cont.cond = length(ig.red) == 1 & length(ig.grn) == 1
      igv.red[!igv.red == ig.red]
      ig.red <- igv.red[1]
    }
    if(cont.cond){
      gpath = c(gpath, gsub("_Red.*", "", ig.red))
    }
    if(verbose){message("Finished validating idat fn's for sample ", i)}
  }

  # get GSM ID indices by interval, as list
  gsmii = getblocks(length(gsmu), gsmint)

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
  rgcni = colnames(t(getRed(rgi))) # grcni = colnames(t(getGreen(rgi)))
  # rgcni == grcni
  rgcn = matrix(c(cn, colnames(t(getRed(rgi)))), nrow = 1)
  data.table::fwrite(rgcn, reds.path, sep = sepval, append = FALSE)
  data.table::fwrite(rgcn, grns.path, sep = sepval, append = FALSE)
  if(getnb){
    if(verbose){message("Instantiating noob beta table...")}
    nbi = minfi::preprocessNoob(rgi)
    nbcn = matrix(c(cn, colnames(t(getBeta(nbi)))), nrow = 1)
    data.table::fwrite(nbcn, nb.path, sep = sepval, append = FALSE)
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
}

# make an HDF5 databse from data tables
make_h5db(dbfnstem = "remethdb",
          version = "0.0.1", ts = ts,
          fnl = list.files("compilations"),
          addmd = TRUE, mdpath = "mdpost_all-gsm-md.rda",
          fnpath = "compilations", rmax = 2)

{

  version = "0.0.1"
  ts = "1123"
  dbfnstem = "remethdb"
  fnl = list.files("compilations")
  rmoldh5 = TRUE
  rmax = 35500
  cmax = 622399
  newtables = TRUE
  verbose = TRUE
  dsnl = c("redsignal", "greensignal", "noobbeta")

  dbn <- paste(paste(dbfnstem, ts,
                       gsub("\\.", "-", version),
                       sep = "_"), "h5",
                 sep = ".")
  try(rhdf5::h5createFile(dbn))
  # remove old data if present
  if(rmoldh5){
    if(verbose){message("Removing any existing old data.")}
    for(d in dsnl){
      try(rhdf5::h5delete(dbn, d))
      try(rhdf5::h5delete(dbn, paste0(d, ".colnames")))
      try(rhdf5::h5delete(dbn, paste0(d, ".rownames")))
    }
  }

  if(verbose){message("Adding and populating data tables to HDF5 database")}
  h5_addtables(fnl = fnl, dsnl = dsnl, rmax = rmax, cmax = cmax, nr.inc = nr.inc)

  {
    dbn = dbn
    fnl = fnl
    dsnl = dsnl
    rmax = 2
    cmax = cmax
    verbose = TRUE
    nr.inc = 10
    fnpath = "compilations"

    for(di in 1:length(dsnl)){
      fnread = fnl[di]; dsn = dsnl[di]
      rhdf5::h5createDataset(dbn, dsn, dims = c(rmax, cmax),
                      maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                      storage.mode = "double", level = 5, chunk = c(1, 5))
      rn = cn = c()
      con <- file(paste(fnpath, fnread, sep = "/"), "r") # make new connection object
      cn = unlist(strsplit(readLines(con, n = 1), " ")) # read first line (colnames)
      cn = cn[2:length(cn)] # filt first value
      cn = gsub("\n", "",gsub('\"', '', cn[1:cmax])) # grab the max indexed value
      nri = getblocks(rmax, nr.inc)
      tt <- Sys.time()
      for(ni in nri){
        i = ni[1]
        # read new lines; note this is contextual/resumes at next new line each iter
        dati = unlist(strsplit(readLines(con, n = length(ni)), " "))
        wdi = which(grepl(".*GSM.*", dati))
        dff = matrix(nrow = 0, ncol = cmax)
        # get stripped GSM IDs from rownames (first col each line)
        ngsm = gsub("\n", "",
                    gsub('\"', '',
                         gsub("\\..*", "", dati[wdi])))
        wgsm = c() # new gsm ids to write
        for(wi in 1:length(wdi)){
          # filter redundant gsm ids
          if(!ngsm[wi] %in% rn){
            wadd = wdi[wi] + 1
            dff = rbind(dff,
                        matrix(dati[wadd:(wadd + cmax - 1)],
                               nrow = 1))
            wgsm = c(wgsm, ngsm[wi])
          }
        }
        rn = c(rn, wgsm) # gsm ids for new data
        class(dff) = "numeric"
        rhdf5::h5write(dff, file = dbn, name = dsn,
                index = list(ni[1]:ni[length(ni)], 1:cmax))
        rhdf5::h5closeAll()
        message("For ds ", dsn,", finished reading index ", i,
                " to ", i + (nr.inc - 1),
                ", time; ", Sys.time() - tt)
      }
      message("Adding row and column names for ds ", dsn)
      cnn = paste0(dsn, ".colnames"); rnn = paste0(dsn, ".rownames");
      rhdf5::h5createDataset(dbn, cnn, dims = length(cn), maxdims = c(rhdf5::H5Sunlimited()),
                      storage.mode = "character", level = 5, # compression level, 1-9
                      chunk = c(20), size = 256) # chunk dims
      message("Added colnames...")
      rhdf5::h5createDataset(dbn, rnn, dims = length(rn), maxdims = c(rhdf5::H5Sunlimited()),
                      storage.mode = "character", level = 5, # compression level, 1-9
                      chunk = c(20), size = 256) # chunk dims
      message("Added rownames...")
      rhdf5::h5write(cn, file = dbn, name = cnn, index = list(1:length(cn)))
      rhdf5::h5write(rn, file = dbn, name = rnn, index = list(1:length(rn)))
      rhdf5::h5closeAll()
      message("Completed writing data, rownames, colnames for ds, ", dsn)
    }
    message("Completed hdf5 write for all ds's in list. Returning")
  }

  if(verbose){message("Finished adding red and green channel data to h5 databse.")}

  if(newtables){
    if(verbose){message("Adding new data tables.")}
    h5_newtables(dbn)


    {
      h5_newtables <- function(dbn, dsn.nb = "noobbeta",
                               dsn.meth = "methylated_signal",
                               dsn.unmeth = "unmethylated_signal",
                               dsn.red = "redsignal", dsn.grn = "greensignal",
                               verbose = TRUE, ngsm.block = 50,
                               ncol.chunk = 5000)

      dbn = "remethdb_1123_0-0-1.h5"
      dsn.nb = "noobbeta"
      dsn.meth = "methylated_signal"
      dsn.unmeth = "unmethylated_signal"
      dsn.red = "redsignal"
      dsn.grn = "greensignal"
      verbose = TRUE
      ngsm.block = 50
      ncol.chunk = 5000


      # Generate noobbeta and meth/unmeth signal tables
      # get dimensions from red and grn signal data
      if(verbose){message("Getting red and green signal h5 object info...")}
      rs.rn <- rhdf5::h5read(dbn, paste0(dsn.red, ".rownames"))
      rs.cn <- rhdf5::h5read(dbn, paste0(dsn.red, ".colnames"))
      gs.rn <- rhdf5::h5read(dbn, paste0(dsn.grn, ".rownames"))
      gs.cn <- rhdf5::h5read(dbn, paste0(dsn.grn, ".colnames"))
      if(verbose){message("Getting blocked indices of row data to process...")}
      sbv <- getblocks(length(rs.rn), ngsm.block)
      # get new cg  dimensions
      if(verbose){message("Getting manifest and setting new colname dim...")}
      anno.name = "IlluminaHumanMethylation450kanno.ilmn12.hg19"
      man = eval(parse(text = paste(anno.name, "Manifest", sep = "::")))
      ncg = nrow(man)
      # new h5 data params
      newdims <- c(length(rs.rn), ncg)
      chunkvars <- c(5, ncol.chunk)
      # make new tables
      rhdf5::h5createDataset(dbn, "unmethylated_signal", dims = newdims,
                             maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                             storage.mode = "double", level = 5, chunk = chunkvars)
      rhdf5::h5createDataset(dbn, "methylated_signal", dims = newdims,
                             maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                             storage.mode = "double", level = 5, chunk = chunkvars)
      rhdf5::h5createDataset(dbn, "noobbeta", dims = newdims,
                             maxdims = c(rhdf5::H5Sunlimited(), rhdf5::H5Sunlimited()),
                             storage.mode = "double", level = 5, chunk = chunkvars)


      # append new data
      if(verbose){message("Appending new data to h5 datasets...")}
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
        se.nb <- minfi::preprocessNoob(se.rgi)
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
          message("finished block ", i," of ",
                  length(sbv),", time elapsed: ",
                  Sys.time() - t1)
        }
      }


    }


  }

  if(verbose){message("Finished all processes. Returning.")}
}

# make HDF5-SummarizedExperiment objects
rmpipeline::make.h5se()


make_h5se("remethdb-seh5", "0.0.1", "1123", se = "rg",
          dbn = "remethdb_1123_0-0-1.h5",
          dsn.data1 = "redsignal", dsn.data2 = "greensignal",
          dsn.rn = "redsignal.rownames",
          dsn.cn = "redsignal.colnames")

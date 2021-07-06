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

#' Aggregate available metadata tables
#'
#' Detects valid metadata at md.dpath, then aggregates these into a single
#' composite table.
#' 
#' @param ts Timestamp.
#' @param platform The DNAm array platform used.
#' @param lfn Optional list of filenames. If NULL, automatically detect 
#' valid metadata files (NULL).
#' @param id.cname Column name for sample IDs, which should be present in 
#' each table to be compiled ("gsm").
#' @param mda.fn Filename stem of new aggregate metadata table to save 
#' ("mdall").
#' @param md.fnv Vector of regex patterns to detect valid metadata tables.
#' @param md.dpath Path to metadata directory.
#' @param verbose Whether to show status messages (TRUE).
#' @return NULL, stores aggregate/composite metadata table.
#' @export
md_agg <- function(ts = NULL, platform = NULL, lfn = NULL, 
                   id.cname = "gsm", mda.fn = "mdall",
                   md.fnv = c("^md_postprocess.*","^mdmod_dnam-predictions_.*",
                              "^mdqc_.*", "^mdrep_.*", "^md_msrapout_.*"),
                   md.dpath = file.path("recount-methylation-files",
                                        "metadata"), verbose = TRUE){
  if(is.null(ts) | is.null(platform)){
    message("Handling metadata options...");md <- rmp_handle_metadata()
    if(is.null(md)){stop("Couldn't get metadata...")}
    version <- md[["version"]]; ts <- md[["timestamp"]]
    message("Getting platform info..."); accinfo <- rmp_handle_platform()
    platform<-accinfo[["platform_name"]];message("Using platform: ",platform)}
  if(is.null(lfn)){
    if(verbose){message("Provided lfn is NULL, detecting metadata files...")}
    lfn <- list.files(md.dpath);lfn <- lfn[grepl(paste0(".*",ts,".*"),lfn)]
    lfn <- lfn[grepl(paste(md.fnv, collapse = "|"), lfn)]}
  if(verbose){message("Loading metadata files...")};ldat <- list()
  for(ii in seq(length(lfn))){
    mdi <- get(load(file.path(md.dpath, lfn[ii])));cnv <- colnames(mdi)
    if(length(cnv[grepl(id.cname, colnames(mdi))]) > 0){
      which.cnid <- which(grepl(id.cname, cnv))
      colnames(mdi)[which.cnid]<-id.cname;ldat[[lfn[ii]]] <- mdi} else{
      message("Couldn't find sample ID column in file ",lfn[ii],
        ". Skipping...")}}
  if(verbose){message("Getting all sample IDs...")}
  idv <- unique(unlist(lapply(ldat, function(x){
    return(as.character(x[,id.cname]))})))
  mda <- as.matrix(data.frame(id = idv, stringsAsFactors = FALSE))
  colnames(mda) <- id.cname; rownames(mda) <- mda[,1]
  if(verbose){
    message("Detected ",length(idv)," unique sample IDs. Coercing tables...")}
  for(ii in seq(length(ldat))){
    tname <- names(ldat)[ii]
    mdt <- as.matrix(ldat[[ii]]); mdt.gsmv <- as.character(mdt[,id.cname])
    idout <- idv[!idv %in% mdt.gsmv]
    if(length(idout) > 0){
      mna <- matrix(rep(rep("NA", ncol(mdt)), length(idout)),
                    nrow=length(idout)); colnames(mna) <- colnames(mdt)
      mna[,id.cname] <- idout; mdt <- rbind(mdt, mna)
      if(verbose){message("Bound ",nrow(mna)," rows to data ",tname)}}
    if(verbose){message("Matching tables on sample IDs...")}
    mdt <- mdt[mdt[,id.cname] %in% mda[,id.cname],]
    mdt <- mdt[order(match(mdt[,id.cname], mda[,id.cname])),]
    rownames(mdt) <- mdt[,id.cname]; 
    mdt <- mdt[,!colnames(mdt) %in% colnames(mda)] # filter redundant vars
    cond <- identical(as.character(rownames(mdt)), 
                      as.character(rownames(mda)))
    if(cond){mda <- cbind(mda, mdt)} else{
      if(verbose){
        message("Couldn't match sample IDs. Skipping table ",tname)}}
    if(verbose){message("Finished with table ",tname)}}
  mda <- as.data.frame(mda, stringsAsFactors = FALSE)
  if(!is.null(platform)){mda$platform <- platform}
  mda.fpath <- file.path(md.dpath, paste0(mda.fn, "_", ts, ".rda"))
  if(verbose){message("Saving aggregate md table to ", mda.fpath)}
  save(mda, file = mda.fpath); return(NULL)
}

#' Append metadata to DNAm data compilations
#'
#' Append metadata to DNAm data compilations. Handles metadata addition for 
#' either HDF5 or HDF5-SummarizedExperiment compilation files.
#' 
#' @param ts Timestamp.
#' @param files.dname Files dir name for the instance 
#' ("recount-methylation-files")
#' @param md.dname Metadata dir name for the instance ("metadata"), contained 
#' in files.dname.
#' @param comp.dname Compilations dir name for instance ("compilations"), contained
#' in files.dname.
#' @param comp.fnv Vector of regular expression patterns to filter compilation 
#' files list.
#' @param mdall.fnstr Filename string stem to detect aggregate metadata table at 
#' md.dname.
#' @param verbose Whether to show status messages (TRUE).
#' @return NULL, if all appends successful, otherwise a list of the compilation 
#' files for which appends weren't successful.
#' @export
append_md <- function(ts, files.dname = "recount-methylation-files", 
                      md.dname = "metadata", comp.dname = "compilations", 
                      comp.fnv = c(".*_h5_.*", ".*_h5se_.*"),
                      mdall.fnstr = "mdall_*", verbose = TRUE){
  comp.dpath <- file.path(files.dname, comp.dname)
  md.dpath <- file.path(files.dname, md.dname);tsstr <- paste0(".*",ts,".*")
  comp.lf <- list.files(comp.dpath); comp.lf <- comp.lf[grepl(tsstr, comp.lf)]
  comp.lf <- comp.lf[grepl(paste(comp.fnv, collapse = "|"), comp.lf)]
  if(length(comp.lf) == 0){
    stop("Didn't find any compilation files at ",comp.dpath,".")}
  md.lf <- list.files(md.dpath); md.lf <- md.lf[grepl(tsstr, md.lf)]
  mdall.fname <- md.lf[grepl(mdall.fnstr, md.lf)][1]
  if(length(mdall.fname) == 0){
    stop("Couldn't find mdall file at ",mdall.fpath,
         ". Try running rule make_md_all first.")}
  mdall <- get(load(file.path(md.dpath, mdall.fname)))
  if(verbose){message("Appending mdall to ",length(comp.lf),
                      " compilation files...")}
  ltry <- list()
  for(comp.fn in comp.lf){
    comp.fpath <- file.path(comp.dpath, comp.fn)
    if(grepl(comp.fnv[1], comp.fn)){
      if(verbose){message("Appending metadata to h5db ", comp.fn, "...")}
      ltry[[comp.fn]] <- try(
        suppressMessages(append_md_h5db(ts = ts, mdall = mdall, 
                                        comp.fpath = comp.fpath)))
    } else{
      if(verbose){message("Appending metadata to h5se file ", comp.fn, "...")}
      ltry[[comp.fn]] <- try(
        suppressMessages(append_md_h5se(ts = ts, mdall = mdall, 
                                        comp.fpath = comp.fpath)))}
    if(verbose){message("Finished compilation file ", comp.fn)}}
  if(length(ltry) == 0){
    if(verbose){"Append success for all compilations."}} else{
      message("Couldn't complete md append for certain files. ",
              "Returning status list...");return(ltry)
    };return(NULL)
}

#' Append metadata to HDF5 DNAm data compilation file
#'
#' Append metadata to DNAm HDF5 data compilations.
#' 
#' @param ts Timestamp.
#' @param mdall Valid metadata table containing column id.cname.
#' @param comp.fpath File path to a valid h5 compilation file.
#' @param id.cname Sample ID column name in mdall ("gsm").
#' @param dsn Name stem of the metadata table in h5db object (e.g. writes table
#' named dsn and column names to paset0(dsn".colnames") to h5db)
#' @return NULL, saves h5 object with appended metadata.
#' @export
append_md_h5db <- function(ts, mdall, comp.fpath, id.cname = "gsm", 
                           dsn = "mdall", overwrite = TRUE){
  dbn <- comp.fpath;tnames <- as.character(rhdf5::h5ls(dbn)[,2])
  if(dsn %in% tnames){
    if(!overwrite){
      stop("Metadata table '",dsn,"' already exists in h5db.")
    } else{
      if(verbose){message("Removing detected metadata table...")}
      rhdf5::h5delete(dbn, dsn);rhdf5::h5delete(dbn, paste0(dsn,".colnames"))}}
  cnstr <- ifelse(grepl(".*h5_rg.*", gsub(".*\\/", "", comp.fpath)),
                  ".*rownames.*", ".*colnames.*")
  tname.sampid <- tnames[grepl(cnstr, tnames)][1]
  idv <- gsub("\\..*", "", as.character(rhdf5::h5read(dbn, tname.sampid)))
  idint <- intersect(idv, mdall[,id.cname])
  if(length(idint) == 0){
    stop("No sample IDs overlap compilation and mdall.")
  } else{
    if(verbose){
      message(length(idint)," samples overlap compilation and mdall.")}}
  mdf <- mdall[mdall[,id.cname] %in% idint, ,drop = FALSE]
  idout <- mdf[!mdf[,id.cname %in% idv],id.cname]
  if(length(idout) > 0){
    if(verbose){
      message("Appending null metadata for ",length(idout)," samples...")}
    mna <- matrix(rep(rep("NA", ncol(mdf)), length(idout)), ncol = ncol(mdf))
    colnames(mna) <- colnames(mdf)
    mna[, which(colnames(mdf) == id.cname)] <- rownames(mna) <- idout
    mdf <- rbind(mdf, mna)}
  mdf <- mdf[order(match(mdf[,id.cname], idv)),]
  cond <- identical(mdf[, id.cname], idv)
  if(cond){
    if(verbose){message("Making metadata objects...")}
    rhdf5::h5createDataset(dbn, dsn, dims = c(nrow(mdall), ncol(mdall)),
                           maxdims = c(rhdf5::H5Sunlimited(), 
                                       rhdf5::H5Sunlimited()),
                           storage.mode = "character", level = 5, 
                           chunk = c(10, 16), size = 256)
    rhdf5::h5createDataset(dbn, paste0(dsn, ".colnames"),level = 5, 
                           dims=ncol(mdall),maxdims=c(rhdf5::H5Sunlimited()),
                           storage.mode = "character",chunk = c(5), size = 256)
    if(verbose){message("Populating new HDF5 entities...")}
    mdf.cnv <- as.character(colnames(mdf))
    mdf <- as.matrix(mdf); class(mdf) <- "character"
    rhdf5::h5write(mdf, file = dbn, name = dsn,
                   index = list(1:nrow(mdf), 1:ncol(mdf)))
    if(verbose){message("Appending metadata...")}
    rhdf5::h5write(mdf.cnv, file = dbn, name = paste0(dsn, ".colnames"),
                   index = list(1:length(mdf.cnv)))};rhdf5::h5closeAll()
  if(verbose){message("Finished appending metadata.")};return(NULL)
}

#' Append metadata to HDF5-SummarizedExperiment DNAm data compilation file
#'
#' Append metadata to DNAm HDF5-SummarizedExperiment data compilations.
#' 
#' @param ts Timestamp.
#' @param mdall Valid metadata table containing column id.cname.
#' @param comp.fpath File path to a valid h5se compilation file.
#' @param id.cname Sample ID column name in mdall ("gsm").
#' @param overwrite Whether to overwrite existing metadata (TRUE).
#' @return NULL, saves h5se object with appended metadata.
#' @export
append_md_h5se <- function(ts, mdall, comp.fpath, id.cname = "gsm", 
                           overwrite = TRUE){
  se <- HDF5Array::loadHDF5SummarizedExperiment(comp.fpath);
  pdat <- pData(se); if(ncol(pdat) > 0){
    if(verbose){message("Overwriting detected metadata")}}
  cnv <- colnames(se); idv <- gsub("\\..*", "", cnv)
  idint <- intersect(idv, mdall[,id.cname])
  if(length(idint) == 0){
    stop("No sample IDs overlap compilation and mdall.")
  } else{
      if(verbose){
        message(length(idint)," samples overlap compilation and mdall.")}}
  mdf <- mdall[mdall[,id.cname] %in% idint, ,drop = FALSE]
  idout <- mdf[!mdf[,id.cname %in% idv],id.cname]
  if(length(idout) > 0){
    if(verbose){
      message("Appending null metadata for ",length(idout)," samples...")}
    mna <- matrix(rep(rep("NA", ncol(mdf)), length(idout)), ncol = ncol(mdf))
    colnames(mna) <- colnames(mdf)
    mna[, which(colnames(mdf) == id.cname)] <- rownames(mna) <- idout
    mdf <- rbind(mdf, mna)}
  mdf <- mdf[order(match(mdf[,id.cname], idv)),]
  cond <- identical(mdf[, id.cname], idv)
  if(cond){
    if(verbose){message("Appending metadata...")}
    rownames(mdf) <- colnames(se);pData(se) <- DataFrame(mdf)
    HDF5Array::quickResaveHDF5SummarizedExperiment(se)
  } else{stop("Couldn't match sample IDs in compilation and mdall.")}
  return(NULL)
}



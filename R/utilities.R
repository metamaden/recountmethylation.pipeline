#!/usr/bin/env R

# Author: Sean Maden
# Utilities for rmpipeline and snakemake workflow.

#------------------
# instance metadata
#------------------
#' Append metadata to HDF5 or SE object
#'
#' Passes version info and timestamp from Python to object metadata
#' @param title Object title
#' @param version Numeric version to be passed, should conform to ##.##.## 
#' nomenclature
#' @param ts Timestamp. If NULL, get new timestamp using Python function.
#' @param pname Name of pipeline package (default: "rmpipeline")
#' @param sname Name of Python script (default: "get_timestamp.y")
#' @return Metadata content for the object
#' @export
get_metadata <- function(title, version, ts = NULL, pname = "rmpipeline",
                         sname = "get_timestamp.py"){
  mdl <- list(title = title, version = version)
  path <- paste(system.file(package = pname), sname, sep="/")
  if(is.null(ts)){
    ts <- system(paste("python", path, sep = " "),intern = TRUE, wait = FALSE)}
  mdl[["timestamp"]] <- ts; return(mdl)
}

#---------------------
# sample block indices
#---------------------

#' Get list of index blocks
#'
#' Get list of index blocks, allowing for remainders.
#'
#' @param slength Total length of index vector
#' @param bsize Size of index blocks along length
#' @return List of index blocks of min length `slength`/`bsize`
#' @export
getblocks <- function(slength, bsize){
  iv <- list()
  if(slength < bsize){
    iv[[1]] <- seq(1, slength, 1)
  } else{
    sc <- 1; ec <- sc + bsize - 1
    nblocks <- slength %/% bsize
    for(b in 1:nblocks){
      iv[[b]] <- seq(sc, ec, 1)
      sc <- ec + 1; ec <- ec + bsize
    }
    # add final indices
    if(nblocks < (slength/bsize)){
      iv[[length(iv) + 1]] <- seq(sc, slength, 1)
    }
  }
  return(iv)
}

#-------------------------------
# get path to edirect query file
#-------------------------------

#' Get an edirect query file path
#'
#' Lookup and return an edirect query file matching pattern `fn.str`.
#'
#'@param fn.str File name string matching the queried file.
#'@param eqdpath Path to directory containing the equery files
#'@return Path to the queried edirect file. 
#'@export
get_eqfpath <- function(fn.str = "gsequery_filt.*",
  eqdpath = file.path("recount-methylation-files", "equery")){
  message("Getting files at eqdpath ", eqdpath, "...")
  eqfl <- list.files(eqdpath);eqfn <- eqfl[grepl(fn.str, eqfl)]
  if(length(eqfn) > 1){
    eqts <- gsub(".*\\.", "", eqfn); eqfn <- eqfn[which(eqts == max(eqts))]
  };eqfpath <- file.path(eqdpath, eqfn);return(eqfpath)
}

#--------------
# Regex helpers
#--------------
# Helper functions for regular expression term matching and mapping

#' Get negative regex string mapping
#'
#' Do negative regex term lookup using patterns mapped by pstr().
#' 
#' @param pstr The mapped regex patterns output from a pstr() lookup.
#' @returns Mapped negative lookup strings.
#' @seealso pstr, get_filt, appendvar
#' @export
get_pstr_neg <- function(pstr){
  pstrg <- gsub("\\.\\*", "_", pstr);pstrg <- gsub("\\|", "", pstrg)
  uv <- unlist(strsplit(pstrg, "_"));uv <- uv[!uv == ""]
  for (ui in 1:length(uv)) {
    s <- uv[ui]
    if(ui == 1){
      ns <- paste(paste0(".*non-", s, ".*"), 
                  paste0(".*Non-", s, ".*"), 
                  paste0(".*non ", s, ".*"), 
                  paste0(".*Non ", s, ".*"), 
                  paste0(".*not ", s, ".*"), 
                  paste0(".*Not ", s, ".*"),
                  paste0(".*Never ", s, ".*"), 
                  paste0(".*never ", s, ".*"), 
                  paste0(".*", s, "-free.*"), 
                  paste0(".*", s, "-Free.*"), 
                  paste0(".*", s, " free.*"), 
                  paste0(".*", s, " Free.*"), 
                  sep = "|")
    } else{
      ns <- paste(ns, 
                  paste0(".*non-", s, ".*"), 
                  paste0(".*Non-", s, ".*"), 
                  paste0(".*non ", s, ".*"), 
                  paste0(".*Non ", s, ".*"), 
                  paste0(".*not ", s, ".*"), 
                  paste0(".*Not ", s, ".*"),
                  paste0(".*Never ", s, ".*"), 
                  paste0(".*never ", s, ".*"),
                  paste0(".*", s, "-free.*"), 
                  paste0(".*", s, "-Free.*"), 
                  paste0(".*", s, " free.*"), 
                  paste0(".*", s, " Free.*"), 
                  sep = "|")
    }
  }
  return(ns)
}

#' Get regex patterns for a string, or the affirmative match for the 
#' given string
#'
#' Gets the regex pattern matching the affirmative of a given string `v`.
#' Does progressive capitalization on values separated by spaces
#' for each value in v, appends flanking '.*' (matches any char)
#' for each value in v, appends "|" OR conditional separator
#'
#' @param v Character string or vector of such character strings. 
#' @returns Regex patterns for detecting matching strings from metadata.
#' @seealso get_pstr_neg, get_filt, appendvar
#' @export
get_pstr <- function(v){
  rs <- ""
  for(ci in seq(length(v))){
    c <- v[ci]
    if(ci == 1){rs = paste0(".*", c, ".*")} else{
      rs = paste(c(rs, paste0(".*", c, ".*")), collapse = "|")
    }
    uv <- unlist(strsplit(c, " ")) # num space-sep units
    # for each unit use lower- and uppercase
    uvstr <- c(paste(uv, collapse = " ")); uvl <- list(uv)
    for(i in seq(length(uv))){
      uvi <- c()
      for(ui in seq(length(uv))){
        chari <- uv[ui]
        if(ui <= i){
          if(nchar(chari)>1){
            ssi <- paste0(toupper(substr(chari, 1, 1)),
                          substr(chari, 2, nchar(chari)))
          } else{
            ssi <- paste0(toupper(substr(chari, 1, 1)))
          }
        }
        else{ssi <- chari}
        uvi <- c(uvi, ssi)
      }; uvl[[length(uvl)+1]] <- uvi
    }
    # append to new str
    for(si in 1:length(uvl)){
      s <- uvl[[si]]
      if(length(uv) > 1){
        if(!si==1){
          # space sep
          rs <- paste(c(rs, paste0(".*", paste(s, collapse = " "), ".*")), collapse = "|")
        }
        # underline sep
        rs <- paste(c(rs, paste0(".*", paste(s, collapse = "_"), ".*")), collapse = "|")
        # dash sep
        rs <- paste(c(rs, paste0(".*", paste(s, collapse = "-"), ".*")), collapse = "|")
      } else{
        if(!si==1){
          rs <- paste(c(rs, paste0(".*", s, ".*")), collapse = "|")
        }
      }
    }
  }
  return(rs)
}

#' Get the negation regex patterns for a given regex pattern
#'
#' Does progressive capitalization on values separated by spaces
#' for each value in v, appends flanking '.*' (matches any char)
#' for each value in v, appends "|" OR conditional separator
#'
#' @param pstr Regex patterns for matching, or output from `get_pstr()`.
#' @returns Regex for negation of indicated patterns.
#' @seealso get_pstr, get_filt, appendvar
#' @export
get_pstr_neg <- function(pstr){
  pstrg = gsub("\\.\\*", "_", pstr); pstrg = gsub("\\|", "", pstrg)
  uv = unlist(strsplit(pstrg, "_")); uv = uv[!uv==""]
  for(ui in 1:length(uv)){
    s = uv[ui]
    if(ui == 1){
      ns = paste(paste0(".*non-", s, ".*"), paste0(".*Non-", s, ".*"), paste0(".*non ", s, ".*"),
                 paste0(".*Non ", s, ".*"), paste0(".*not ", s, ".*"), paste0(".*Not ", s, ".*"),
                 paste0(".*", s, "-free.*"), paste0(".*", s, "-Free.*"), paste0(".*", s, " free.*"),
                 paste0(".*", s, " Free.*"), sep = "|")
    } else{
      ns = paste(ns, paste0(".*non-", s, ".*"), paste0(".*Non-", s, ".*"), paste0(".*non ", s, ".*"),
                 paste0(".*Non ", s, ".*"), paste0(".*not ", s, ".*"), paste0(".*Not ", s, ".*"),
                 paste0(".*", s, "-free.*"), paste0(".*", s, "-Free.*"), paste0(".*", s, " free.*"),
                 paste0(".*", s, " Free.*"), sep = "|")
    }
  }
  return(ns)
}

#' Get the outcome of pattern matching a string with regex patterns
#'
#' Does pattern matching for regex patterns constructed from a character
#' string or vector of such character strings. 
#'
#' @param v Character string or vector of such strings.
#' @param m Preprocessed metadata used for pattern matching/lookups (data.frame, mdpre).
#' @param filtrel Logical symbol joining each regex pattern (default "|").
#' @param ntfilt Regex pattern corresponding to negative lookup filter (default NULL).
#' @param ptfilt Regex pattern corresponding to positive lookup filter (default NULL).
#' @returns The result of assessing a regex match on a metadata variable.
#' @seealso appendvar
#' @export
get_filt <- function(v, m = mdpre, filtrel = "|", ntfilt = NULL, ptfilt = NULL,
                     varl = c("gsm_title", "sample_type", "disease_state", 
                              "anatomic_location", "misc")){
  # positive lookup filter
  filtl <- grepl(v, m[,varl[1]])
  if(!is.null(ptfilt)){
    message("Using positive lookup filter with ptfilt: ", 
            paste(ptfilt, collapse = "; "))
    filtl <- filtl & grepl(get_pstr(ptfilt), m[,varl[1]])
  }
  # negative lookup filter
  if(!is.null(ntfilt)){
    message("Using negative lookup filter with ntfilt: ", 
            paste(ntfilt, collapse = "; "))
    nfiltv <- get_pstr_neg(ntfilt)
    filtl <- filtl & !grepl(nfiltv, m[,varl[1]])
  }
  # proceed if additional vars specified
  if(length(varl)>1){
    if(!filtrel %in% c("|", "&")){
      message("Please provide a valid filter relation symbol.")
      return(NULL)
    } else{
      for(vi in varl[2:length(varl)]){
        if(filtrel == "|"){
          filtl <- filtl | grepl(v, m[,vi])
          if(!is.null(ntfilt)){filtl <- filtl & !grepl(nfiltv, m[,vi])}
        } else if(filtrel == "&"){
          filtl <- filtl & grepl(v, m[,vi])
          if(!is.null(ntfilt)){filtl <- filtl & !grepl(nfiltv, m[,vi])}
        }
      }
    }
  }
  return(filtl)
}

#' Append mapped terms to a metadata variable
#'
#' Appends new variable data to a metadata variable, preserving current
#' variable terms.
#'
#' @param var Variable in m for which to append new terms.
#' @param val Character string to append to matched samples.
#' @param filtv Vector of boolean values identifying samples for which to
#'  append the new term.
#' @param m The postprocessed metadata for which to append the new terms
#'  (data.frame, mdpost).
#' @returns The result of appending new terms to specified var.
#' @seealso get_filt, get_pstr_neg, get_pstr
#' @export
appendvar <- function(var, val, filtv, m = mdpost){
  varr <- m[, var]
  # get composite filter
  filti <- !grepl(paste0("(^|;)", val, "(;|$)"), varr)
  compfilt <- filti & filtv
  # assess filter results
  if(length(compfilt[compfilt]) == 0){
    message("No unique values to append. Returning var unchanged.")
    return(varr)
  } else{
    varr[compfilt] <- ifelse(varr[compfilt] == "NA", val,
                             paste(varr[compfilt], val, sep = ";")
    )
    message("Appended n = ", length(varr[compfilt]), " values")
    return(varr)
  }
  return(NULL)
}

#' Match and combine 2 datasets
#'
#' Orders 2 datasets on corresponding variables, then appends an NA
#' matrix as needed before binding on columns. This is mainly used to 
#' combine different types of sample metadata on common/shared GSM IDs.
#' Note columns should be matchable after calling `as.character()` for
#' each.
#'
#' @param d1 First set to combine (matrix).
#' @param d2 Second set to combine (matrix).
#' @param ci1 Column index in first dataset to match on
#' @param ci2 Column index in second dataset to match on
#' @returns Union of the 2 datasets, including NA values where appropriate.
#' @export
match1to2 <- function(d1, d2, ci1 = 2, ci2 = 1){
  id.all <- unique(c(d1[, ci1], d2[, ci2]))
  id1 <- id.all[!id.all %in% d1[, ci1]]
  id2 <- id.all[!id.all %in% d2[, ci2]]
  # append na slices as necessary
  if(length(id1) > 0){
    nav <- rep(rep("NA", length(id1)), ncol(d1) - 1)
    mna <- matrix(c(id1, nav), nrow = length(id1), ncol = ncol(d1))
    d1 <- rbind(d1, mna)
  }
  if(length(id2) > 0){
    nav <- rep(rep("NA", length(id2)), ncol(d2) - 1)
    mna <- matrix(c(id2, nav), nrow = length(id2), ncol = ncol(d2))
    d2 <- rbind(d2, mna)
  }
  # reorder and assign title var
  match.id1 <- match(as.character(d1[, ci1]), as.character(d2[, ci2]))
  order.id1 <- order(match.id1); d1 <- d1[order.id1,]
  match.id2 <- match(as.character(d2[, ci2]), as.character(d1[, ci1]))
  order.id2 <- order(match.id2); d2 <- d2[order.id2,]
  cond <- identical(as.character(d2[, ci2]), as.character(d1[, ci1]))
  if(cond){return(cbind(d1, d2))} else{
    stop("there was an issue matching the final datasets...")
  }
  return(NULL)
}

#------------------------
# handle instance options
#------------------------

#' Get instance metadata info
#'
#' Looks up and loads previously generated instance metadata 
#' (e.g. using rule `new_instance_md`), which includes the version and 
#' timestamp for the `recountmethylation` instance.
#'
#' @param files.dname
#' @param md.dname
#' @return Metadata object containing the instance version and timestamp
#' @export
rmp_handle_metadata <- function(files.dname = "recount-methylation-files",
                                md.dname = "metadata"){
  md.dpath <- file.path(files.dname, md.dname)
  if(!dir.exists(md.dpath)){stop("Error, could find dir ", md.dpath,"...")}
  lfmd <- list.files(md.dpath); lfmd <- lfmd[grepl("metadata", lfmd)]
  if(length(lfmd) > 0){
    ts <- as.numeric(gsub(".*_|\\.rda", "", lfmd)); ts.max <- max(ts)
    which.md <- which(ts == ts.max)[1]
    md.fpath <- file.path(md.dpath, lfmd[which.md])
    message("Using metadata at ", md.fpath);md <- get(load(md.fpath))
    return(md)
  } else{
    stop("Error, couldn't find metadata for this instance.", 
         " Try running the `new_inst_md` rule first.")}
  return(NULL)
}
#' Get instance platform info
#'
#' Retrieves the previously established instance platform info (e.g. using 
#' rule `set_acc`) from the settings script `settings.py` .
#'
#' @param settings.fname Name of the settings script containing platform info.
#' @param src.path Path the `src` dir in the `recountmethylation_server` repo.
#' @return Platform info as list contianing the accession ID and name.
#' @export
rmp_handle_platform <- function(settings.fname = "settings.py",
                                src.path=file.path("recountmethylation_server",
                                                     "src")){
  settings.fpath = file.path(src.path, settings.fname)
  if(!file.exists(settings.fpath)){stop("Error, couldn't find ",settings.fn)}
  accid <- gsub(" |.* '|'", "", readLines(settings.fpath, n = 25)[25])
  pname <- ifelse(accid == "GPL13534", "hm450k",
                  ifelse(accid == "GPL21145", "epic-hm850k",
                         ifelse(accid == "GPL8490", "hm27k", "NA")))
  return(list(accid = accid, platform_name = pname))
}
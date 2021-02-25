#!/usr/bin/env R

# BeadArray controls functions formatted for R package checks.

#------------------
# Control functions
#------------------

#' Get the background address from control probe data
#'
#' @param cdf Control probe information (data.frame).
#' @param verbose Whether to show status messages (TRUE).
#' @examples 
#' cdf <- system.file("extdata", "controldf", "cgcontroldf.rda", 
#' package = "recountmethylation")
#' get_background(cdf)
#' @export
get_background <- function(cdf, verbose = TRUE){
  if(verbose){message("Getting background addresses...")};ct <- "EXTENSION"
  ci <- cdf[cdf$Type==ct,];which.bkg <- grepl("(A)|(T)", ci[,"ExtendedType"])
  return(ci[which.bkg,]$Address)
}

#' Restoration metric
#'
#' Single metric, using grn channel
#' 
#'
#' @param bg.addr Probe address of background control probe.
#' @param rm Matrix of control metric signals.
#' @param cname Column name of the control probe.
#' @param gs Green signal data (data.frame, columns are probes, rows are samples, 
#' column names are addresses, rownames are samples/GSM IDs).
#' @param cdf Control probe annotations (data.frame, cols = properties, rows = probes).
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @param verbose Whether to show status messages (TRUE).
#' @return Control restoration signal.
#' @examples 
#' cdf <- system.file("extdata", "controldf", "cgcontroldf.rda", 
#' package = "recountmethylation")
#' @export
control_restoration <- function(bg.addr, rm, cname = "Restore", 
                                gs, cdf, baseline = 3000, biotin.baseline = 1,
                                verbose = TRUE){
  if(verbose){message("Calculating Restore metric...")}
  ct <- "RESTORATION"; ci <- cdf[cdf$Type == ct,]
  which.rest <- which(colnames(gs) %in% ci$Address)
  m1 <- apply(gs, 1, function(x){
    x[which.rest]/(max(x[bg.addr]) + baseline)})
  rm <- cbind(rm, m1);colnames(rm)[ncol(rm)] <- cname;return(rm)
}

#' Biotin staining metrics
#'
#' Two metrics, 1 per color channel
#' 
#' @param rm Matrix of control metric signals.
#' @param rs Red signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param gs Green signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param cdf Control probe annotations (data.frame, cols = properties, 
#' rows = probes).
#' @param cnames Vector of control probe column names for the metric.
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @param verbose Whether to show status messages (TRUE).
#' @return Control biotin staining signal.
#' @examples 
#' cdf <- system.file("extdata", "controldf", "cgcontroldf.rda", 
#' package = "recountmethylation")
#' 
control_biotinstaining <- function(rm, rs, gs, cdf, 
                                   cnames = c("Biotin_red", "Biotin_green"), 
                                   baseline = 3000, biotin.baseline = 1, 
                                   verbose = TRUE){
  ci <- cdf[grepl("Biotin|DNP", cdf$ExtendedType),]
  if(verbose){message("Calculating Biotin Staining from red signal...")}
  which.stain <- which(colnames(rs) %in% 
                         ci[grepl("DNP \\(High", ci$ExtendedType),]$Address)
  which.bkg <- which(colnames(rs) %in% 
                       ci[grepl("DNP \\(Bkg", ci$ExtendedType),]$Address)
  m1 <- apply(rs,1,function(x){x[which.stain]/(x[which.bkg]+biotin.baseline)})
  rm <- cbind(rm, m1); colnames(rm)[ncol(rm)] <- cnames[1]
  if(verbose){message("Calculating Biotin Staining from green signal...")}
  which.stain <- which(colnames(gs) %in% 
                         ci[grepl("Biotin \\(High", ci$ExtendedType),]$Address)
  which.bkg <- which(colnames(gs) %in% 
                       ci[grepl("Biotin \\(Bkg", ci$ExtendedType),]$Address)
  m2 <- apply(gs,1,function(x){x[which.stain]/(x[which.bkg]+biotin.baseline)})
  rm <- cbind(rm, m2); colnames(rm)[ncol(rm)] <- cnames[2]; return(rm)
}

#' Specificity metrics
#'
#' Three metrics (two for I and one for type II), for each color channel and 
#' BeadArray probe type.
#'
#' @param rm Matrix of control metric signals.
#' @param rs Red signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param gs Green signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param cdf Control probe annotations (data.frame, cols = properties, 
#' rows = probes).
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param cnames Vector of control probe column names for the metric.
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @param verbose Whether to show status messages (TRUE).
#' @return Control specificity metrics.
#' @examples 
#' cdf <- system.file("extdata", "controldf", "cgcontroldf.rda", 
#' package = "recountmethylation")
#' @export
control_specificity <- function(rm, rs, gs, cdf, baseline = 3000, 
                                cnames = c("Specificity_I_red", 
                                           "Specificity_I_grn", 
                                           "Specificity_II"),
                                biotin.baseline = 1, verbose = TRUE){
  ct <- "SPECIFICITY I"; ci1 <- cdf[cdf$Type == ct,]
  if(verbose){message("Calculating Specificity I from red signal...")}
  ci <- ci1[grepl("Mismatch (4|5|6)", ci1$ExtendedType),]
  addr.mm.index <- which(colnames(rs) %in% 
                           ci[grepl("MM", ci$ExtendedType),]$Address)
  addr.pm.index <- which(colnames(rs) %in% 
                           ci[grepl("PM", ci$ExtendedType),]$Address)
  m1 <- apply(rs, 1, function(x){min(x[addr.pm.index])/max(x[addr.mm.index])})
  rm <- cbind(rm, m1); colnames(rm)[ncol(rm)] <- cnames[1]
  if(verbose){message("Calculating Specificity I from green signal...")}
  ci <- ci1[grepl("Mismatch (1|2|3)", ci1$ExtendedType),]
  addr.mm.index <- which(colnames(gs) %in% 
                           ci[grepl("MM", ci$ExtendedType),]$Address)
  addr.pm.index <- which(colnames(gs) %in% 
                           ci[grepl("PM", ci$ExtendedType),]$Address)
  m2 <- apply(gs, 1, function(x){min(x[addr.pm.index])/max(x[addr.mm.index])})
  rm <- cbind(rm, m2); colnames(rm)[ncol(rm)] <- cnames[2]
  if(verbose){message("Calculating SpecificityII metric...")}
  ct <- "SPECIFICITY II"; ci <- cdf[cdf$Type == ct,]
  which.addr.red <- which(colnames(rs) %in% ci$Address)
  which.addr.grn <- which(colnames(gs) %in% ci$Address)
  m0.1 <- apply(rs, 1, function(x){min(x[which.addr.red])})
  m0.2 <- apply(gs, 1, function(x){max(x[which.addr.grn])})
  m1 <- m0.1/m0.2; rm <- cbind(rm,m1); colnames(rm)[ncol(rm)] <- cnames[3]
  return(rm)
}

#' Extension metrics
#'
#' Two metrics, one per color channel.
#' 
#' @param rm Matrix of control metric signals.
#' @param rs Red signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param gs Green signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param cdf Control probe annotations (data.frame, cols = properties, 
#' rows = probes).
#' @param cnames Vector of control probe column names for the metric.
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @param verbose Whether to show status messages (TRUE).
#' @return Control extension metrics.
#' @examples 
#' cdf <- system.file("extdata", "controldf", "cgcontroldf.rda", 
#' package = "recountmethylation")
#' @export
control_extension <- function(rm, rs, gs, cdf, verbose = TRUE,
                              cnames = c("Extension_red","Extension_green"),
                              baseline = 3000, biotin.baseline = 1){
  ct <- "EXTENSION"; ci <- cdf[cdf$Type == ct,]
  addr.ext.cg <- ci$Address[grepl("(C)|(G)", ci$ExtendedType)]
  addr.ext.at <- ci$Address[grepl("(A)|(T)", ci$ExtendedType)]
  if(verbose){message("Calculating Extension metric from red signal...")}
  which.cg <- which(colnames(rs) %in% addr.ext.cg)
  which.at <- which(colnames(rs) %in% addr.ext.at)
  m1 <- apply(rs, 1, function(x){min(x[which.at])/max(x[which.cg])})
  if(verbose){message("Calculating Extension metric from green signal...")}
  which.cg <- which(colnames(gs) %in% addr.ext.cg)
  which.at <- which(colnames(gs) %in% addr.ext.at)
  m2 <- apply(gs, 1, function(x){min(x[which.cg])/max(x[which.at])})
  rm <- cbind(rm, m1); colnames(rm)[ncol(rm)] <- cnames[2]
  rm <- cbind(rm, m2); colnames(rm)[ncol(rm)] <- cnames[2]
  return(rm)
}

#' Hybridization metrics
#'
#' Two metrics, one comparing high vs. medium level, the other comparing medium 
#' vs. low level. Both metrics calculated from Green color channel signals.
#'
#' @param rm Matrix of control metric signals.
#' @param rs Red signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param gs Green signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param cdf Control probe annotations (data.frame, cols = properties, 
#' rows = probes).
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param cnames Vector of control probe column names for the metric.
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @param verbose Whether to show status messages (TRUE).
#' @return Control hybridization metrics.
#' @examples 
#' cdf <- system.file("extdata", "controldf", "cgcontroldf.rda", 
#' package = "recountmethylation")
control_hybridization <- function(rm, rs, gs, cdf, baseline = 3000, 
                                  cnames = c("Hybridization_medium-vs-high",
                                             "Hybridization_low-vs-medium"),
                                  biotin.baseline = 1, verbose = TRUE){
  ct <- "HYBRIDIZATION"; ci <- cdf[cdf$Type == ct,]
  which.hi <- which(colnames(gs) %in% 
                      ci$Address[grepl("(High)", ci$ExtendedType)])
  which.med <- which(colnames(gs) %in% 
                       ci$Address[grepl("(Medium)", ci$ExtendedType)])
  which.low <- which(colnames(gs) %in% 
                       ci$Address[grepl("(Low)", ci$ExtendedType)])
  if(verbose){
    message("Calculating Hybridization metrics for high versus medium...")}
  m1 <- apply(gs, 1, function(x){x[which.hi]/x[which.med]})
  if(verbose){
    message("Calculating Hybridization metrics for medium versus low...")}
  m2 <- apply(gs, 1, function(x){x[which.med]/x[which.low]})
  rm <- cbind(rm, m1); colnames(rm)[ncol(rm)] <- cnames[1]
  rm <- cbind(rm, m2); colnames(rm)[ncol(rm)] <- cnames[2]; return(rm)
}

#' Target removal metrics
#'
#' Two metrics, one for each BeadArray probe type, each calculated from green
#'  color channel signals.
#'
#' @param rm Matrix of control metric signals.
#' @param rs Red signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param gs Green signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param cdf Control probe annotations (data.frame, cols = properties, 
#' rows = probes).
#' @param cnames Vector of control probe column names for the metric.
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @param verbose Whether to show status messages (TRUE).
#' @return Control target removal metrics.
#' @examples 
#' cdf <- system.file("extdata", "controldf", "cgcontroldf.rda", 
#' package = "recountmethylation")
#' @export
control_targetremoval <- function(rm, rs, gs, cdf, verbose = TRUE,
                                  cnames = c("Target_removal_1",
                                             "Target_removal_2"),
                                  baseline = 3000, biotin.baseline = 1){
  ct <- "EXTENSION"; ci <- cdf[cdf$Type == ct,]
  which.bkg <- which(colnames(gs) %in% 
                       ci[grepl("(A)|(T)", ci$ExtendedType), ]$Address)
  ct <- "TARGET REMOVAL"; ci <- cdf[cdf$Type == ct,]
  if(verbose){message("Calculating Target Removal 1 metrics...")}
  which.t1 <- which(colnames(gs) %in% 
                      ci[grepl("Removal 1", ci$ExtendedType),]$Address)
  if(verbose){message("Calculating Target Removal 2 metrics...")}
  which.t2 <- which(colnames(gs) %in% 
                      ci[grepl("Removal 2", ci$ExtendedType),]$Address)
  m1 <- apply(gs, 1, function(x){(max(x[which.bkg]) + baseline)/x[which.t1]})
  m2 <- apply(gs, 1, function(x){(max(x[which.bkg]) + baseline)/x[which.t2]})
  rm <- cbind(rm, m1); colnames(rm)[ncol(rm)] <- cnames[1]
  rm <- cbind(rm, m2); colnames(rm)[ncol(rm)] <- cnames[2]; return(rm)
}

#' Bisulfite conversion I and II
#'
#' Three metrics. For Bisulfite conversion I, calculates one metric from 
#' controls C1-3 vs. U1-3 in green channel, the other from controls C4-6 vs. 
#' U4-6 in red channel. For Bisulfite conversion II, calculates one metric.
#'
#' @param rm Matrix of control metric signals.
#' @param rs Red signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param gs Green signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param cdf Control probe annotations (data.frame, cols = properties, 
#' rows = probes).
#' @param cnames Vector of control probe column names for the metric.
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @param verbose Whether to show status messages (TRUE).
#' @return Control bisulfite conversion metrics.
#' @examples 
#' cdf <- system.file("extdata", "controldf", "cgcontroldf.rda", 
#' package = "recountmethylation")
#' @export
control_bsconversion <- function(rm, rs, gs, cdf, baseline = 3000,
                                 cnames = c("Bisulfite_conversion_I_red",
                                            "Bisulfite_conversion_I_green",
                                            "Bisulfite_conversion_II"),
                                 biotin.baseline = 1, verbose = TRUE){
  if(verbose){message("Setting up Bisulfite conversion calculations...")}
  ct <- "BISULFITE CONVERSION I"; ci <- cdf[cdf$Type == ct,]
  which.c123 <- which(colnames(gs) %in% 
                        ci$Address[grepl("C1|C2|C3", ci$ExtendedType)])
  which.u123 <- which(colnames(gs) %in% 
                        ci$Address[grepl("U1|U2|U3", ci$ExtendedType)])
  which.c456 <- which(colnames(rs) %in% 
                        ci$Address[grepl("C4|C5|C6", ci$ExtendedType)])
  which.u456 <- which(colnames(rs) %in% 
                        ci$Address[grepl("U4|U5|U6", ci$ExtendedType)])
  if(verbose){
    message("Calculating Bisulfite Conversion I metrics, red channel...")}
  m1 <- apply(rs, 1, function(x){min(x[which.c456])/max(x[which.u456])})
  if(verbose){
    message("Calculating Bisulfite Conversion I metrics, green channel...")}
  m2 <- apply(gs, 1, function(x){min(x[which.c123])/max(x[which.u123])})
  rm <- cbind(rm, m1); colnames(rm)[ncol(rm)] <- cnames[1]
  rm <- cbind(rm, m2); colnames(rm)[ncol(rm)] <- cnames[2]
  if(verbose){message("Calculating Bisulfite Conversion II metric...")}
  ct <- "BISULFITE CONVERSION II"; ci <- cdf[cdf$Type == ct,]
  which.ci.red <- which(colnames(rs) %in% ci$Address)
  which.ci.grn <- which(colnames(gs) %in% ci$Address)
  m0.1 <- apply(rs, 1, function(x){min(x[which.ci.red])})
  m0.2 <- apply(gs, 1, function(x){max(x[which.ci.grn])})
  m1 <- m0.1/m0.2; rm <- cbind(rm,m1); colnames(rm)[ncol(rm)] <- cnames[3]
  return(rm)
}

#' Non-polymorphic metric
#'
#' Two metrics, one per color channel.
#'
#' @param rm Matrix of control metric signals.
#' @param rs Red signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param gs Green signal data (data.frame, columns are probes, rows are 
#' samples, column names are addresses, rownames are samples/GSM IDs).
#' @param cdf Control probe annotations (data.frame, cols = properties, 
#' rows = probes).
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param cnames Vector of control probe column names for the metric.
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @param verbose Whether to show status messages (TRUE).
#' @return Control non-polymorphic metrics.
#' @examples 
#' cdf <- system.file("extdata", "controldf", "cgcontroldf.rda", 
#' package = "recountmethylation")
#' @export
control_nonpolymorphic <- function(rm, rs, gs, cdf, baseline = 3000, 
                                   cnames = c("Non-polymorphic_red",
                                              "Non-polymorphic_green"),
                                   biotin.baseline = 1,verbose = TRUE){
  ct <- "NON-POLYMORPHIC"; ci <- cdf[cdf$Type == ct,]
  if(verbose){
    message("Calculating Non-polymorphic metric for red channel...")}
  which.cg <- which(colnames(rs) %in% ci$Address[grepl("(C)|(G)", 
                                                       ci$ExtendedType)])
  which.at <- which(colnames(rs) %in% ci$Address[grepl("(A)|(T)", 
                                                       ci$ExtendedType)])
  m1 <- apply(rs, 1, function(x){min(x[which.at])/max(x[which.cg])})
  if(verbose){
    message("Calculating Non-polymorphic metric for green channel...")}
  which.cg <- which(colnames(gs) %in% ci$Address[grepl("(C)|(G)", 
                                                       ci$ExtendedType)])
  which.at <- which(colnames(gs) %in% ci$Address[grepl("(A)|(T)", 
                                                       ci$ExtendedType)])
  m2 <- apply(gs, 1, function(x){min(x[which.cg])/max(x[which.at])})
  rm <- cbind(rm, m1); colnames(rm)[ncol(rm)] <- cnames[1]
  rm <- cbind(rm, m2); colnames(rm)[ncol(rm)] <- cnames[2];return(rm)
}



#---------------
# main functions
#---------------
# Main functions coordinating control signal calculations, assessments

#' Get control probe metrics from matrices of red/grn signal 
#'
#' @description
#'This resource takes guidance from 3 key sources: the Illumina GenomeStudio Methylation 
#'Module (v1.8, source 1), the BeadArray Controls Reporter Software Guide (v00, source 2), 
#'and the ewastools resource (v1.5, source 3). Notes on how function relates to these sources 
#'follows: 
#'
#' * Use C and U 1-3 and 4-6 for Bisulfite Conversion I per source 2, source 3 only uses 
#'    1-2/4-5, and may have string matching error
#'    
#' * Use probe address "34648333" and "43603326" (DNP, Biotin subtype) for Biotin 
#'    'Staining Background where appropriate, per sources 2 and 3; 
#' 
#' * To offset many Inf/-Inf values due to 0 background signal, set a special 
#'    Biotin specific baseline offset (default: 1) in denom of Biotin red/grn 
#'    Bkg value. No offset specified in source 2;
#' 
#' * For Specificity I, use PM, MM 1-3 for grn, 4-6 for red (per source 3, 
#'    source 2 ambiguous);
#' 
#' * For Specificity II, use just probes S1-3, as probe S4 unavailable in 
#'    control probe annotation (per source 2, corroborated in source 3)
#' 
#' * Note definition for background is ambiguous (can be either Red or Grn 
#'    channel, from extension), per source 2; source 3 uses extension Grn 
#'    A/T probes for system background
#' 
#' * Note additional source 2 metrics using background (e.g. specificity 
#'    background/U1-3) currently not calculated by `bactrl()`.
#'  
#' @param rgset RGChannelSet objet containing control probe signals.
#' @param baseline Baseline measure for signals (integer, 3000).
#' @param biotin.baseline Baseline to use for biotin controls (integer, 1).
#' @param verbose Whether to show status messages (TRUE, boolean).
#' @param cnv Vector of control assay names for columns in returned table.
#' @returns Table of BeadArray control assay signals.
#' @seealso bathresh
#' @export
get_bactrl <- function(rgset,baseline=3000,biotin.baseline=1,verbose=TRUE,
                       cnv = c("restoration.grn","biotin.stain.red", 
                               "biotin.stain.grn", "specificityI.red", 
                               "specificityI.grn", "specificityII",
                               "extension.red", "extension.grn", 
                               "hyb.hi.med", "hyb.med.low", "target.removal.1",
                               "target.removal.2", "bisulfite.conv.I.red",
                               "bisulfite.conv.I.grn", "bisulfite.conv.II", 
                               "nonpolymorphic.red", "nonpolymorphic.grn")){
  if(verbose){message("Getting control probe info and red/green signals...")}
  cdf <- as.data.frame(minfi::getProbeInfo(rgset, type = "Control"))
  rs <- t(minfi::getRed(rgset)); gs <- t(minfi::getGreen(rgset))
  rm <- data.frame(gsm = rownames(rs)); lctrl <- list()
  if(verbose){message("Getting background probe addresses...")}
  bg <- get_background(cdf = cdf, verbose = verbose)
  if(verbose){message("Getting control assay signals...")}
  rm <- control_restoration(bg.addr = bg, rm = rm, gs = gs, cdf = cdf,
                            baseline = baseline, verbose = verbose,
                            biotin.baseline = biotin.baseline)
  rm <- control_biotinstaining(rm=rm,rs=rs,gs=gs,cdf=cdf,baseline=baseline, 
                               verbose=verbose,biotin.baseline=biotin.baseline)
  rm <- control_specificity(rm=rm,rs=rs,gs=gs,cdf=cdf,baseline=baseline, 
                            verbose = verbose, biotin.baseline=biotin.baseline)
  rm <- control_extension(rm=rm,rs=rs,gs=gs,cdf=cdf,baseline=baseline, 
                          verbose = verbose, biotin.baseline = biotin.baseline)
  rm <- control_hybridization(rm=rm,rs=rs,gs=gs,cdf=cdf,baseline=baseline,
                              verbose=verbose,biotin.baseline=biotin.baseline)
  rm <- control_targetremoval(rm=rm,rs=rs,gs=gs,cdf=cdf,baseline=baseline, 
                              verbose=verbose,biotin.baseline=biotin.baseline)
  rm <- control_bsconversion(rm=rm,rs=rs,gs=gs,cdf=cdf,baseline=baseline, 
                             verbose=verbose,biotin.baseline=biotin.baseline)
  rm <- control_nonpolymorphic(rm=rm,rs=rs,gs=gs,cdf=cdf,baseline=baseline, 
                               verbose=verbose,biotin.baseline=biotin.baseline)
  colnames(rm) <- c("gsmid", cnv); return(rm)
}

#' Get minimum thresholds for BeadArray controls
#' 
#' Get a table of the minimum thresholds for BeadArray control assessments. 
#' These are called by get_bathresh() to do the quality assessments (e.g. 
#' if signal < threshold -> fail, else pass).
#'
#' @return Data frame containing minimum thresholds for BeadArray controls.
#' @seealso get_bathresh(); get_bactrl()
#' @export
df_bathresh <- function(){
  dft=data.frame(restoration.grn=0, biotin.stain.red=5, biotin.stain.grn=5, 
                 specificityI.red=1, specificityI.grn=1, specificityII=1, 
                 extension.red=5, extension.grn=5, hyb.hi.med=1, 
                 hyb.med.low=1, target.removal.1=1, target.removal.2=1, 
                 bisulfite.conv.I.red=1, bisulfite.conv.I.grn=1, 
                 bisulfite.conv.II=1, nonpolymorphic.red=5, 
                 nonpolymorphic.grn=5, stringsAsFactors=F);return(dft)
}

#' Get BeadArray control outcomes from a matrix of metric signals
#'
#' Apply BeadArray thresholds to determine passing and failing samples. 
#' This function takes a matrix of signals and performs controls using 
#' minimum signal thresholds from the BeadArray Controls Reporter 
#' Software Guide ewastools.
#'
#' @param rm BeadArray signals returned by `bactrl()` (matrix, rows = samples, 
#' cols = metrics).
#' @param pass.lab Label for passing samples ("PASS", chracter).
#' @param fail.lab Label for failing samples ("FAIL", chracter).
#' @returns Matrix of threshold assessments, either 'FAIL' or 'PASS', 
#' of same dimensions as input matrix.
#' @seealso bactrl
#' @export
get_bathresh <- function(rm, pass.lab = "PASS", fail.lab = "FAIL"){
  dft <- df_bathresh()
  for(c in colnames(rm)[2:ncol(rm)]){
    rm[,c]<-ifelse(rm[,c]<dft[,c],fail.lab,pass.lab)};return(rm)}

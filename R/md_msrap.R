
# get mapped terms from metasra-pipeline
get_msrap <- function(files.dir = "recount-methylation-files", 
                      eq.dname = "equery", msrap.dname = "gsm_msrap_outfiles"){
  # get paths
  eqdpath = file.path(files.dir, eq.dname) # equery dir path
  eq.dpath <- get_eqfpath(eqdpath) # equery file path
  eq.fn <- gsub(".*/", "", eq.fpath) # equery file name
  msrap.dpath <- file.path(files.dir, msrap.dname) # msrap dir path
  
  # get collapsed gse ids
  x <- scan(eq.path, what="", sep="\n");dat1 <- unlist(strsplit(x, " "))
  gsmt <- data.frame(gsm = dat1[2:length(dat1)]);gsmt$gse <- dat1[1]
  for(i in x[2:length(x)]){
    dati <- unlist(strsplit(i, " "))
    gsmti <- data.frame(gsm = dati[2:length(dati)])
    gsmti$gse <- dati[1];gsmt <- rbind(gsmt, gsmti)}
  
  # collapse on repeat gsm/multiple gse id's
  gsev <- unlist(lapply(gsmt[,1], function(x){
    paste0(gsmt[gsmt[,1] == x, 2], collapse=";")}))
  gsmt$gseid <- gsev
  
  # filter vars, repeated values
  gmap <- gsmt[!duplicated(gsmt$gsm), c(1, 2)]
  
  # get msrap mapped terms
  # append collapsed mapped metadata
  gsmid.fnindex <- 3; gmap$msrap_flatjson <- "NA"
  lfl.filt <- list.files(dpath)
  for(f in lfl.filt){
    gsmi <- unlist(strsplit(f,"\\."))[gsmid.fnindex]
    if(gsmi %in% gmap[,1]){
      jfi <- jsonlite::fromJSON(txt = file.path(dpath, f), flatten=T)
      xi <- jfi[,!colnames(jfi) %in% c("mapped ontology terms", "real-value properties")]
      pi <- paste0("'",names(xi),"':'",as.character(xi[1,]),"'",collapse=";")
      if("real-value properties" %in% names(jfi)){
        if(length(unlist(jfi[,"real-value properties"]))>0){
          pi <- paste0(pi,";'real-value properties':",
                       as.character(unlist(jfi[,"real-value properties"])), 
                       collapse="")}}
      mot <- paste0(as.character(unlist(jfi[,"mapped ontology terms"])), collapse=";")
      jfi.mapped.flat <- paste0(pi, ";", mot)
      gmap[gmap[,1]==gsmi,]$msrap_flatjson <- jfi.mapped.flat}}
  # save data
  gmap.fn <- "mdmap-gsm"
  save(gmap, file=paste0(gmap.fn, ".rda"))
  write.csv(gmap, file=paste0(gmap.fn, ".csv"))
  return(NULL)
}
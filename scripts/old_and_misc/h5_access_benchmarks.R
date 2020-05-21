#--------------
# vignette code
#--------------
{
  dbn = "remethdb2.h5"
  
  mdp = data.mdpost(dbn, dsn = "mdpost")
  
  dsn = "redsignal"
  rs.gsm = rhdf5::h5read(dbn, paste(dsn, "rownames", sep = ".")) # redsignal rownames (gsm ids)
  rs.gsm = gsub("\\..*", "", rs.gsm)
  mdf = mdp[mdp$gsm %in% rs.gsm,]
  unique(unlist(strsplit(mdf$tissue, ";"))) # available tissue terms
  
  # get sample id query by tissue term
  termi = "blood"
  var.query = "tissue"
  which.index = 1:2
  which.gsm = which(grepl(paste0(".*", termi, ".*"), 
                          mdf[,var.query]))[which.index]
  gsmvi = mdf$gsm[which.gsm]
  
  anno.name = "IlluminaHumanMethylation450kanno.ilmn12.hg19"
  man = eval(parse(text = paste(anno.name, "Manifest", sep = "::")))
  loc = eval(parse(text = paste(anno.name, "Locations", sep = "::")))
  identical(rownames(loc), rownames(man))
  
  chrname = "chr9"
  cgfilt = rownames(loc[grepl(paste0("^", chrname, "$"), loc$chr),])
  cgvi = unique(c(man[cgfilt,]$AddressA, man[cgfilt,]$AddressB))
  
  se.rgi = getrg(gsmv = gsmvi, cgv = cgvi, dbn = dbn, data.type = "se", 
                 metadata = T, verbose = T, md.dsn = "mdpost")
  
}

#-----------------------------
# args and data for benchmarks
#-----------------------------
{
  dbn = "remethdb2.h5"
  mdp = data.mdpost(dbn, dsn = "mdpost")
  dsn = "redsignal"
  
  # get gsm ids with available signal
  rs.gsm = rhdf5::h5read(dbn, paste(dsn, "rownames", sep = ".")) # redsignal rownames (gsm ids)
  rs.gsm = gsub("\\..*", "", rs.gsm)
  
  # num GSMs to query
  ngsm = c(2, seq(5, 20, 5), seq(20, 80, 10))
}

#--------------------------
# benchmark random samples
#--------------------------
{
  # get timing data
  dtr = matrix(nrow = 0, ncol = 2)
  for(ngi in ngsm){
    gsmvi = rs.gsm[sample(length(rs.gsm), ngi)]
    tdi = system.time(getrg(gsmv = gsmvi, cgv = "all", dbn = dbn, data.type = "se", 
                            metadata = T, verbose = T, md.dsn = "mdpost"))
    dtr = rbind(dtr, matrix(c(length(gsmvi), tdi[3]), nrow = 1))
    message("finished ngi ", ngi)
  }
}

#------------------------------
# benchmark consecutive samples
#------------------------------
{
  # get timing data
  dtc = matrix(nrow = 0, ncol = 2)
  
  for(ngi in ngsm){
    gsmvi = rs.gsm[1:ngi]
    tdi = system.time(getrg(gsmv = gsmvi, cgv = "all", dbn = dbn, data.type = "se", 
                            metadata = T, verbose = T, md.dsn = "mdpost"))
    dtc = rbind(dt, matrix(c(length(gsmvi), tdi[3]), nrow = 1))
    message("finished ngi ", ngi)
  }
  
  pdf("bmquery_gsm-consecutive.pdf", 6, 5)
  {
    par(mfrow = c(2, 1), 
        oma = c(3, 3, 3, 0), 
        mar = c(3, 5, 1, 1))
    plot(dtc[,1], dtc[,2], type = "o",
         xlab = "Num. GSMs", ylab = "Time Elapsed (sec)")
    plot(dtc[,1], dtc[,2]/dt[,1], type = "o",
         xlab = "Num. GSMs", ylab = "Time Per GSM (sec)")
    
    mtext("Consecutive GSM Query", side = 3, outer = T)
    mtext("Num. GSMs", side = 1, outer = T)
  }
  dev.off()
}

#------------------------------
# overlay plots with error bars
#------------------------------
# utilities from https://davetang.org/muse/2014/06/25/plotting-error-bars-with-r/
{
  # function for error bars
  error.bar <- function(x, y, upper, lower=upper, length=0.1,...){
    if(length(x) != length(y) | length(y) !=length(lower) | length(lower) != length(upper))
      stop("vectors must be same length")
    arrows(x,y+upper, x, y-lower, angle=90, code=3, length=length, ...)
  }
  
  # function for standard error of the mean
  sem <- function(x){
    sd(x)/sqrt(length(x))
  }
}

# calculate means
my_mean <- apply(df, 2, mean)

# calculate sem
my_sem <- apply(df, 2, sem)

pdf("bmquery_gsm-consecutive.pdf", 6, 5)
{
  par(mfrow = c(2, 1), 
      oma = c(3, 3, 3, 0), 
      mar = c(3, 5, 1, 1))
  plot(dtc[,1], dtc[,2], type = "o", col = "blue",
       xlab = "Num. GSMs", ylab = "Time Elapsed (sec)")
  lines(dtc[,1], dtc[,2], col = "red", type = "o")
  
  plot(dtc[,1], dtc[,2]/dtc[,1], type = "o", col = "blue",
       xlab = "Num. GSMs", ylab = "Time Per GSM (sec)")
  lines(dtr[,1], dtr[,2]/dtr[,1], col = "red", type = "o")
  
  mtext("Consecutive GSM Query", side = 3, outer = T)
  mtext("Num. GSMs", side = 1, outer = T)
}
dev.off()
library(rhdf5)

#-----------------
# make the test db
#-----------------
fnl = c("newred.comp", "1553728068.rawgrn.compilation.mdat", "1553728595.noobbeta.compilation.mdat")
dsnl = c("redsignal", "greensignal", "noobbeta")
dbn = "remethdb_test.h5"
h5createFile(dbn)
rmax = 10; cmax = 10 # 622399
# delete old tables
for(d in dsnl){
  h5delete(dbn, d)
  h5delete(dbn, paste0(d, ".colnames"))
  h5delete(dbn, paste0(d, ".rownames"))
}

h5ds.add = function(fnl, dsnl, rmax, cmax, 
                    nr.inc = 10, dbn = "remethdb_test.h5"){
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
h5ds.add(fnl = fnl, dsnl = dsnl, rmax = 100, cmax = 1000, nr.inc = 20)

# test h5 write results
dsn = "redsignal"; whichcol = c(1, 10); whichrow = c(1, 50)
hi = h5read(dbn, dsn, index = list(whichrow[1]:whichrow[2], whichcol[1]:whichcol[2]))
dff = NULL
con <- file(fnread, "r")
cn = readLines(con, n = 1)
dati = unlist(strsplit(readLines(con, n = whichrow[2]), " "))
wdi = which(grepl(".*GSM.*", dati))
dff = matrix(nrow = 0, ncol = whichcol[2])
# make the new data slice
for(w in wdi){
  wadd = w+1
  dff = rbind(dff, matrix(dati[wadd:(wadd + whichcol[2] - 1)], nrow = 1))
}
class(dff) = "numeric"
close(con)
identical(dff, hi)

hi = h5read(dbn, dsn, index = list(1:10, 1:10))
di = dff

r1 = h5read(dbn, paste0("redsignal", ".rownames"), index = list(whichrow[1]:whichrow[2]))
r2 = h5read(dbn, paste0("greensignal", ".rownames"), index = list(whichrow[1]:whichrow[2]))
r3 = h5read(dbn, paste0("noobbeta", ".rownames"), index = list(whichrow[1]:whichrow[2]))
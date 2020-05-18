library(HDF5Array)

# run from dir containing full h5se-gr object

path <- "remethdb_h5se_gr_00-00-01_1583780004"

gr <- loadHDF5SummarizedExperiment(path)

# subset gsmids
gsmv <- c("GSM1038308", "GSM1038309")
gsm.patt <- paste0("^", gsmv, ".*", collapse = "|")
which.gsm <- grepl(gsm.patt, colnames(gr))
grsub <- gr[,which.gsm]

head(getBeta(grsub))

# subset probes

new.dn <- "remethdb_h5se-test_gr_00-00-01_1583780004"
saveHDF5SummarizedExperiment(, dir = new.dn)
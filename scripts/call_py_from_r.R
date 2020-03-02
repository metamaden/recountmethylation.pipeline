devtools::install_github("metamaden/rmpipeline")

# option 1
pname <- "rmpipeline"
sname <- "get_timestamp.py"
path <- paste(system.file(package=pname), sname, sep="/")
command <- paste("python", path, sname, "to", sep = " ")
try(suppressWarnings(response <- system(command, intern=T)), silent = T)

# option 2
library(reticulate)

path <- paste(system.file(package=pname), sname, sep="/")
py_run_file(path)
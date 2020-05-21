devtools::install_github("metamaden/rmpipeline")

# option 1
pname <- "rmpipeline"
sname <- "get_timestamp.py"
path <- paste(system.file(package=pname), sname, sep="/")
command <- paste("python", path, sname, "to", sep = " ")
ts <- system(command, intern=T)

# option 2
library(reticulate)
pname <- "rmpipeline"
sname <- "get_timestamp.py"
path <- paste(system.file(package=pname), sname, sep="/")
reticulate::py_run_file(path)
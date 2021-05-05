# recountmethylation.pipeline

[![DOI](https://zenodo.org/badge/241738988.svg)](https://zenodo.org/badge/latestdoi/241738988)

Functions and utilities to compile and harmonize public GEO DNA methylation array datasets. This includes 
generation of the assay compilation objects as HDF5 and HDF5-SummarizedExperiment objects, as well as 
metadata mappings and predictions for sample metadata from study SOFT files. To see implementation, see the 
[`recountmethylation_instance` resource](https://github.com/metamaden/recountmethylation_instance).

# Install

## Devtools install 

Install this package from an R session using the devtools R package with:

```
devtools::install_github("metamaden/rmpipeline")
``` 

## Manual install

Obtain this package from github with:

```
git clone https://www.github.com/metamaden/recountmethylation.pipeline
```

Install the cloned repo using:

```
R CMD INSTALL recountmethylation.pipeline
```
## Load

After installation, load this package into an active R session `library(rmpipeline)`. 

#!/usr/bin/env R

# Author: Sean Maden
# 
# Snakemake script to prepare DNAm assays data.

library(recountmethylation.pipeline)

dnam_pipeline <- function(){
  get_rg_dtables()
  get_h5db_rg()
  get_h5se_rg()
  get_h5db_gm()
  get_h5se_gm()
  get_h5db_gr()
  get_h5se_gr()
  return(NULL)
}

suppressMessages(dnam_pipeline())
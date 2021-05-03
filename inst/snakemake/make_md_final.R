#!/usr/bin/env R

# Author: Sean Maden
# 
# Snakemake script to produce all sample metadata for a recountmethylation 
# instance. This includes metadata mapping and harmonization, DNAm-based metadata,
# and MetaSRA-pipeline sample type predictions.

library(recountmethylation.pipeline)

suppressMessages()
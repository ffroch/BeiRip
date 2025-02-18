##########################################################################
# This script creates phyloseq objects from the raw data
# Part B: Fungal data
##########################################################################
# load libraries
library(phyloseq)
library(tidyr)
library(tidyverse)
library(ape)
##########################################################################
# PART B: relative abundance phylo object
###########################################################################
# calculate relative abundance
phy <- readRDS("./03-Routputfiles/D_Phyloobjects/D2A_fun_raw_allsamples.rds")
phyrel <- transform_sample_counts(phy, function(x) x / sum(x))
###########################################################################
# save phyloseq object
saveRDS(phyrel, "./03-Routputfiles/D_Phyloobjects/D2B_fun_rel_allsamples.rds")
# remove the mock and negative control samples
phyrel <- subset_samples(phyrel, piece != "mo" & piece != "NK")
# save phyloseq object
saveRDS(phyrel, "./03-Routputfiles/D_Phyloobjects/D2B_fun_rel_nocontrols.rds")
###########################################################################
###########################################################################
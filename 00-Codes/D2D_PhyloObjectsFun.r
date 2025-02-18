##########################################################################
# This script creates phyloseq objects from the raw data
# Part D: Fungal data
##########################################################################
# load libraries
library(phyloseq)
library(tidyr)
library(tidyverse)
##########################################################################
# PART B: create a phylo object on genus level
##########################################################################
phyabs1 <- readRDS("./03-Routputfiles/D_Phyloobjects/D2C_fun_abs_nocontrols.rds")
phygen <- tax_glom(phyabs1, taxrank = "Genus")
# save phyloseq object
saveRDS(phygen, "./03-Routputfiles/D_Phyloobjects/D2D_fun_abs_nc_genus.rds")
##########################################################################
##########################################################################
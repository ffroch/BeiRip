##########################################################################
# This script creates phyloseq objects from the raw data
# Part D: absolute abundance phylo object
        # contains the absolute abundance of the 16S rRNA gene sequencing
        # corrected by the spike in bacteria Imtechella
        # does NOT contain the mock community samples and the NC
        # 2 files with and without Imtechella and Allobacillus
##########################################################################
# load libraries
library(phyloseq)
library(tidyr)
library(tidyverse)
library(ape)
##########################################################################
##########################################################################
# PART D: absolute abundance phylo object
# Calculate absolute abundance of all bacteria
##########################################################################
phy_rel <- readRDS("./03-Routputfiles/D_Phyloobjects/D1B_bac_rel_nocontrols_spike.rds")
# get otu table
otuf <- otu_table(phy_rel)
# get meta table
metaf <- sample_data(phy_rel)
# check if rownames of otu table and colnames of meta table are identical
identical(colnames(otuf), rownames(metaf))
# calculate absolute gene copies for each ASV based on Imtechella
otufabs <- sweep(otuf, 2, metaf$gencopiesbasedonimt, "*")
# correct for predilution before DNA extraction
otufabs <- sweep(otufabs, 2, metaf$dilpreextr, "/")
# create a new phyloseq object with absolute gene copies
phyabs <- phy_rel
otu_table(phyabs) <- otu_table(otufabs, taxa_are_rows = TRUE)
# remove samples with mo in piece
phyabs <- subset_samples(phyabs, piece != "mo")
# remove samples with NK in piece
phyabs <- subset_samples(phyabs, piece != "NK")
# remove features with 0 abundance
phyabs <- prune_taxa(taxa_sums(phyabs) > 0, phyabs)
##########################################################################
# make clean taxonomy, there are a lot of underscores with additional
# numbers and letters in the names, that should be removed
# check taxonomy generally
tax <- data.frame(tax_table(phyabs))
tax[, 1:6] <- data.frame(lapply(tax[, 1:6], function(x) gsub("_.*", "", x)))
# replace tax table in phyloseq object
tax_table(phyabs) <- tax_table(as.matrix(tax))
saveRDS(phyabs, "./03-Routputfiles/D_Phyloobjects/D1D_bac_abs_nocontrols_spike.rds")
##########################################################################
# remove spike-in bacteria
phyabs <- subset_taxa(phyabs, Genus != "Imtechella" & Genus != "Allobacillus")
saveRDS(phyabs, "./03-Routputfiles/D_Phyloobjects/D1D_bac_abs_nocontrols_nospike.rds")
##########################################################################
##########################################################################

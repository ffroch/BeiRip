##########################################################################
# This script creates phyloseq objects from the raw data
# Part B: relative abundance phylo object + Spike-in bacteria abundance
        # contains the relative abundance of the 16S rRNA gene sequencing
        # phyloseq object is saved as phy_rel and phy_rel_tea
        # does contain also the spike in bacteria Imtechella and Allobacillus
        # 2 files, with and without mock community samples and negative Controls
##########################################################################
# load libraries
library(phyloseq)
library(tidyr)
library(tidyverse)
library(ape)
##########################################################################
# PART B: relative abundance phylo object + Spike-in bacteria abundance
# Calculate rel. abundance of Spike-in bacteria (Imtechella, Allobacillus)
##########################################################################
# import phyloseq object
phy <- readRDS("./03-Routputfiles/D_Phyloobjects/D1A_bac_raw_allsamples_spike.rds")
# calculate relative abundance
phy_rel <- transform_sample_counts(phy, function(x) x / sum(x))
##########################################################################
# Extract relative abundance of "Imtechella" and "Allobacillus" for each sample
# access taxonomy and abundance data
taxtable <- as.data.frame(tax_table(phy_rel))
otutable <- as.data.frame(otu_table(phy_rel))
# subset forthe genus "Imtechella"
taxtable[taxtable$Genus ==  "Imtechella", ]
imt <- which(taxtable$Genus ==  "Imtechella")
# calculate imtechella abundance per sample
imt_abundance <- data.frame(colSums(otutable[imt, ]))
colnames(imt_abundance) <- "imt_abundance"
imt_abundance$sampleid <- rownames(imt_abundance)
# subset forthe genus "Allobacillus"
taxtable[taxtable$Genus ==  "Allobacillus", ]
all <- which(taxtable$Genus ==  "Allobacillus")
# calculate allobacillus abundance per sample
all_abundance <- data.frame(colSums(otutable[all, ]))
colnames(all_abundance) <- "all_abundance"
all_abundance$sampleid <- rownames(all_abundance)
##########################################################################
# Calculate total 16S gene copies per sample
##########################################################################
metadf <- data.frame(sample_data(phy_rel))
# add imtechella and allobacillus abundance to metadata
metadf <- merge(metadf, imt_abundance, by = "sampleid", all.x = TRUE)
metadf <- merge(metadf, all_abundance, by = "sampleid", all.x = TRUE)
meta <- metadf
rownames(meta) <- meta$sampleid
# calculate the 16s gene copies
meta$imt16scopies <- 6e7 / 2e7 * meta$spikecells
meta$all16scopies <- 1.4e8 / 2e7 * meta$spikecells
# calculate absolute cell number per sample
meta$absbasedonimpt <- meta$imt16scopies * (1 / meta$imt_abundance - 1) / 3
meta$absbasedonall <- meta$all16scopies * (1 / meta$all_abundance - 1) / 7
meta$gencopiesbasedonimt <- meta$imt16scopies * (1 / meta$imt_abundance)
meta$gencopiesbasedonall <- meta$all16scopies * (1 / meta$all_abundance)
sample_data(phy_rel) <- meta
# reduce sample names
sample_names(phy_rel) <- substr(sample_names(phy_rel), 1, 8)
# reduce names in sample_data sampleid
sample_data(phy_rel)$sampleid <- substr(sample_data(phy_rel)$sampleid, 1, 8)
# change - to _ in sample names
sample_names(phy_rel) <- gsub("-", "_", sample_names(phy_rel))
# replace - with _ in sample_data sampleid
sample_data(phy_rel)$sampleid <- gsub("-", "_", sample_data(phy_rel)$sampleid)
# replace - with _ in sample_data internID
sample_data(phy_rel)$internID <- gsub("-", "_", sample_data(phy_rel)$internID)
# save phy_rel 
saveRDS(phy_rel, "./03-Routputfiles/D_Phyloobjects/D1B_bac_rel_allsamples_spike.rds")
# remove mock community samples and negative controls
phy_rel <- subset_samples(phy_rel, piece != "mo" & piece != "NK")
# save phy_rel filtered
saveRDS(phy_rel, "./03-Routputfiles/D_Phyloobjects/D1B_bac_rel_nocontrols_spike.rds")
##########################################################################
##########################################################################

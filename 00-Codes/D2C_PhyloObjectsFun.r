##########################################################################
# This script creates phyloseq objects from the raw data
# Part B: Fungal data
##########################################################################
# load libraries
library(phyloseq)
library(tidyr)
library(tidyverse)
library(ape)
library(zoo)
##########################################################################
# PART B: absolute abundance phylo object
# calculate absolute abundance of fungi based on qPCR data
###########################################################################
phyrel <- readRDS("./03-Routputfiles/D_Phyloobjects/D2B_fun_rel_nocontrols.rds")
# remove the mock samples
phyrel <- subset_samples(phyrel, !grepl("mock", sample_names(phyrel)))
# import primer matching asvs
targets0mm <- readRDS("./03-Routputfiles/Z_18SqPCR/matchingasvsall0mm.rds")
targets1mm <- readRDS("./03-Routputfiles/Z_18SqPCR/matchingasvsall1mm.rds")
# import qpcr results
qpcr <- readRDS("./03-Routputfiles/qpcrsampleresults.rds")
colnames(qpcr)[1] <- "internID"
# remove duplicates from qpcr based on internID
qpcr <- qpcr[!duplicated(qpcr$internID), ]
# add FR- before internID
qpcr$internID <- paste("FR_", qpcr$internID, sep = "")
# replace - with _ in internID
qpcr$internID <- gsub("-", "_", qpcr$internID)
###########################################################################
# extract otu table
otudf <- data.frame(otu_table(phyrel))
# keep only the rows that are in targets with 0 mismatches
otudffil0mm <- otudf[rownames(otudf) %in% targets0mm, ]
# calculate the sum of the columns
matchsums0mm <- colSums(otudffil0mm)
# make matchsums a data frame
matchsums0mm <- data.frame(internID = names(matchsums0mm), matchsums0mm)
# keep only the rows that are in targets with 1 mismatches
otudffil1mm <- otudf[rownames(otudf) %in% targets1mm, ]
# calculate the sum of the columns
matchsums1mm <- colSums(otudffil1mm)
# make matchsums a data frame
matchsums1mm <- data.frame(internID = names(matchsums1mm), matchsums1mm)
plot(matchsums1mm$matchsums1mm, matchsums0mm$matchsums0mm)

matchsums <- data.frame(matchsums0mm, matchsums1mm = matchsums1mm$matchsums1mm)
# merge matchsums with qpcr
qpcr <- plyr::join(qpcr, matchsums, by = "internID", type = "left")

qpcr$totalcopies0mm <- qpcr$MeanConc/qpcr$matchsums0mm
qpcr$totalcopies1mm <- qpcr$MeanConc/qpcr$matchsums1mm
# get meta data from phyloseq object
meta <- data.frame(sample_data(phyrel))
# keep only internID and days
meta1 <- meta[, c("internID", "days", "dilpreextr")]
# merge meta with qpcr
qpcr <- plyr::join(qpcr, meta1, by = "internID", type = "left")
# correct for predilution before DNA extraction and correct for elution volume
# elution volume was 35µl 5µl were used for qPCR = factor 7
# the 35µl corresponded to 1.5g meat sample to get the total copies per gram meat factor 1.5
qpcr$totalcopies0mm <- qpcr$totalcopies0mm/qpcr$dilpreextr*7/1.5
qpcr$totalcopies1mm <- qpcr$totalcopies1mm/qpcr$dilpreextr*7/1.5

##########################################################################
# add total copies to meta
meta <- plyr::join(meta, qpcr[, c("internID", "totalcopies0mm", "totalcopies1mm")], by = "internID", type = "left")
rownames(meta) <- meta$sampleid
# get otu table
# check if rownames of otu table and colnames of meta table are identical
identical(colnames(otudf), rownames(meta))
# calculate absolute gene copies for each ASV based qPCR data
otudfabs <- sweep(otudf, 2, meta$totalcopies0mm, "*")
# create a new phyloseq object with absolute gene copies
phyabs <- phyrel
# replace - with . in sample names
otu_table(phyabs) <- otu_table(otudfabs, taxa_are_rows = TRUE)
phyabsextra <- phyabs
phyabsextra <- sample_data(meta)
saveRDS(phyabsextra, "./03-Routputfiles/D_Phyloobjects/D2C_fun_abs_nocontrols_unfiltered.rds")
# remove samples from phyabs with Inf or NA in meta$totalcopies
phyabs1 <- subset_samples(phyabs, !is.na(meta$totalcopies0mm) & !is.infinite(meta$totalcopies0mm))
# remove features with 0 abundance
phyabs1 <- prune_taxa(taxa_sums(phyabs1) > 0, phyabs1)
# add meta data to phyabs1
sample_data(phyabs1) <- sample_data(meta)
# save phyabs1 as rds
saveRDS(phyabs1, "./03-Routputfiles/D_Phyloobjects/D2C_fun_abs_nocontrols.rds")
###########################################################################


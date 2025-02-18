##########################################################################
# This script creates phyloseq objects from the raw data
# Part C: copy number corrected version of the raw data phylo object
        # contains the copy number corrected version of the 16S rRNA gene sequencing
        # phyloseq object is saved as phycnc
        # does contain also the spike in bacteria Imtechella and Allobacillus
        # does contain the mock community samples and the negative Controls
##########################################################################
# load libraries
library(phyloseq)
library(tidyr)
library(tidyverse)
library(ape)
library(zoo)
##########################################################################
# functions
# write a function called find lowest taxonomic level
find_lowest_present_taxon <- function(row, copynr) {
  for (rank in rev(names(row))) {
    if (row[[rank]] %in% copynr$name) {
      return(row[[rank]])
    }
  }
  return(NA) # Return NA if no match is found
}
##########################################################################
##########################################################################
# PART C: copy number correction of the raw data phylo object
##########################################################################
phycnc <- readRDS("./03-Routputfiles/D_Phyloobjects/D1A_bac_raw_allsamples_spike.rds")
# get and transpose otu table
otufcnct <- data.frame(t(otu_table(phycnc)))
# colnames of data frame can not start with a number
# add an X to all colnames of otufcnct if there is not alread an X
colnames(otufcnct) <- ifelse(grepl("X", colnames(otufcnct)),
                             colnames(otufcnct),
                             paste0("X", colnames(otufcnct)))
# calculate column sums
colsums <- colSums(otufcnct)
# order the columns by column sums, to get the most abundant ASVs first
otufcnct <- otufcnct[, order(colsums, decreasing = TRUE)]
##########################################################################
# get meta data
metafcnc <- data.frame(sample_data(phycnc))
names <- rownames(metafcnc)
# add colsumsraw to metafcnc
# check if rownames of otu table and meta data are identical
identical(rownames(otufcnct), rownames(metafcnc))
# combine otu table and meta data
comb <- cbind(metafcnc, otufcnct)
# get the start index of the first asv
start <- ncol(metafcnc) + 1
##########################################################################
# get taxonomy
taxx <- data.frame(tax_table(phycnc))
# add an X to all rownames of taxonomy if there is not alread an X
rownames(taxx) <- ifelse(grepl("X", rownames(taxx)),
                         rownames(taxx), paste0("X", rownames(taxx)))
# reomve rows of taxonomy that are not in comb
taxx <- taxx[rownames(taxx) %in% colnames(comb)[start:ncol(comb)], ]
# sort rows of taxonomy by column order in comb
taxx <- taxx[match(colnames(comb)[start:ncol(comb)], rownames(taxx)), ]
rownamestaxx <- rownames(taxx)
# check if rownames of otu table and taxonomy are identical
identical(colnames(comb)[start:ncol(comb)], rownames(taxx))
##########################################################################
# get known copy numbers for each species from rrnDB
copynr <- read.table("./02-rawdata/rrnDB-5.8_pantaxa_stats_NCBI.tsv",
                     sep = "\t", header = TRUE)
# remove the column sum16slist from copynr
copynr <- copynr[, 1:10]
##########################################################################
# make taxonomy and copynr compatible
##########################################################################
# species is not cleaned yet and has to be cleaned
# replace everything between "_" and " " with " " in taxx
# like Leuconostoc_B inhae
taxx <- data.frame(lapply(taxx, function(x) sub("_.* ", " ", x)))
##########################################################################
# combine copynr and taxx with repeated search
##########################################################################
taxx1 <- taxx
##########################################################################
# make a column called ASV in taxx1
taxx1$ASV <- taxx1$sseqid
# find the lowest present taxon and add it to name
taxx1$name <- apply(taxx1[1:7], 1, function(row) find_lowest_present_taxon(row, copynr))
# join taxx1 and copynr by name
taxx1 <- dplyr::left_join(taxx1, copynr, by = "name")
# correct copy number for rank superkingdom to 5
taxx1$mean[taxx1$rank == "superkingdom"] <- 5
# correct name for Unassigned
taxx1$name[taxx1$Species == "Unassigned"] <- "Unassigned"
# correct copy number for Unassigned to 5
taxx1$mean[taxx1$name == "Unassigned"] <- 5
##########################################################################
# make a column called label in taxx2
taxx1$label <- paste(taxx1$ASV, "\n",
                     taxx1$Species, "\n",
                     taxx1$rank, "\ncorrected for copy numbers:",
                     taxx1$mean)
##########################################################################
# rename rownames of taxx2 by rownames of taxx
rownames(taxx1) <- taxx1$sseqid
##########################################################################
# recalculate the counts for each ASV based on the
# mean copy number of the species
comb2 <- comb
for (i in start:ncol(comb2)) {
  # get mean copy number
  gencop <- taxx1$mean[i - (start - 1)]
  # calculate counts
  comb2[, i] <- comb2[, i] / gencop
}
newotu <- t(comb2[, start:ncol(comb2)])
# remove X from rownames of newotu
rownames(newotu) <- gsub("X", "", rownames(newotu))
otu_table(phycnc) <- otu_table(newotu, taxa_are_rows = TRUE)
tax_table(phycnc) <- tax_table(as.matrix(taxx1))
# save phycorr
saveRDS(phycnc, "./03-Routputfiles/D_Phyloobjects/D1C_bac_rawcnc_allsamples_spike.rds")
# remove mock community samples and negative controls
phycnc <- subset_samples(phycnc, piece != "mo" & piece != "NK")
# save phycnc
saveRDS(phycnc, "./03-Routputfiles/D_Phyloobjects/D1C_bac_rawcnc_nocontrols_spike.rds")
saveRDS(taxx1, "./03-Routputfiles/D1C_CNCtaxonomy.rds")
##########################################################################
##########################################################################

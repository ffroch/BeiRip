##########################################################################
# This script creates phyloseq objects from the raw data
# Part E: copy number corrected absolute data phylo object
        # contains the absolute abundance of the 16S rRNA gene sequencing
        # corrected mean copy numbers based on rrnDB
        # phyloseq object is saved as phycorr
        # does NOT contain the mock community samples and the NC
        # does NOT contain the spike in bacteria Imtechella and Allobacillus
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
# PART E: copy number corrected absolute data phylo object
##########################################################################
phyabs <- readRDS("./03-Routputfiles/D_Phyloobjects/D1D_bac_abs_nocontrols_nospike.rds")
# get and transpose otu table
otufabst <- data.frame(t(otu_table(phyabs)))
# colnames of data frame can not start with a number
# add an X to all colnames of otufabst if there is not alread an X
colnames(otufabst) <- ifelse(grepl("X", colnames(otufabst)),
                             colnames(otufabst),
                             paste0("X", colnames(otufabst)))
# calculate column sums
colsums <- colSums(otufabst)
# order the columns by column sums, to get the most abundant ASVs first
otufabst <- otufabst[, order(colsums, decreasing = TRUE)]
##########################################################################
# get meta data
metafabs <- data.frame(sample_data(phyabs))
names <- rownames(metafabs)
# add colsumsraw to metafabs
# check if rownames of otu table and meta data are identical
identical(rownames(otufabst), rownames(metafabs))
# combine otu table and meta data
comb <- cbind(metafabs, otufabst)
# get the start index of the first asv
start <- ncol(metafabs) + 1
##########################################################################
# get taxonomy
taxx <- data.frame(tax_table(phyabs))
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
phycorr <- phyabs
otu_table(phycorr) <- otu_table(newotu, taxa_are_rows = TRUE)
tax_table(phycorr) <- tax_table(as.matrix(taxx1))
# save phycorr
saveRDS(taxx1, "./03-Routputfiles/D1E_CNCtaxonomy.rds")
saveRDS(phycorr, "./03-Routputfiles/D_Phyloobjects/D1E_bac_abscnc_nocontrols_nospike.rds")
##########################################################################
##########################################################################

##########################################################################
# This script loads the bacterial and fungal absolute abundance phyloseq objects
# transform it into relative abundance data
# log transform the data for normalization
# rescale the values using Z-score scaling
# combine the data into one phyloseq object
# and save the phyloseq object
##########################################################################
# load libraries
library(phyloseq)
library(tidyr)
library(tidyverse)
library(ape)
##########################################################################
# PART A: Combined data for CoNet and Cytoscape
##########################################################################
# load data
# for bacteria the copy number corrected absolute data was used
# for fungi the absolute data was used
phybac <- readRDS("./03-Routputfiles/D_Phyloobjects/D1E_bac_abscnc_nocontrols_nospike.rds")
phyfun <- readRDS("./03-Routputfiles/D_Phyloobjects/D2C_fun_abs_nocontrols.rds")
# there are some samples in the bacterial that are not in the fungal data
# remove samples that are not in both phyloseq objects
common <- intersect(sample_names(phybac), sample_names(phyfun))
phybac2 <- prune_samples(sample_names(phybac) %in% common, phybac)
phyfun2 <- prune_samples(sample_names(phyfun) %in% common, phyfun)
##########################################################################
# transform the data into relative abundance
phybacrel <- transform_sample_counts(phybac2, function(x) x/sum(x))
phyfunrel <- transform_sample_counts(phyfun2, function(x) x/sum(x))
##########################################################################
# log transform the data
phybaclog <- transform_sample_counts(phybacrel, function(x) log10(x + 1))
phyfunlog <- transform_sample_counts(phyfunrel, function(x) log10(x + 1))
##########################################################################
# rescale the data using Z-score scaling
bac_otu <- otu_table(phybaclog)
fun_otu <- otu_table(phyfunlog)
bac_otu_scaled <- apply(bac_otu, 2, scale)
fun_otu_scaled <- apply(fun_otu, 2, scale)
rownames(bac_otu_scaled) <- rownames(bac_otu)
rownames(fun_otu_scaled) <- rownames(fun_otu)
phybacz <- phybaclog
otu_table(phybacz) <- otu_table(bac_otu_scaled, taxa_are_rows = taxa_are_rows(phybaclog))
phyfunz <- phyfunlog
otu_table(phyfunz) <- otu_table(fun_otu_scaled, taxa_are_rows = taxa_are_rows(phyfunlog))
# check if the rescaling worked
identical(taxa_names(phybaclog), taxa_names(phybacz))
identical(taxa_names(phyfunlog), taxa_names(phyfunz))
##########################################################################
# remove the trees from the phyloseq objects
# otherwise the merge_phyloseq function will not work
phybacz <- phyloseq(otu_table(phybacz),
                    sample_data(phybacz),
                    tax_table(phybacz))
phyfunz <- phyloseq(otu_table(phyfunz),
                    sample_data(phyfunz),
                    tax_table(phyfunz))
phycomb <- merge_phyloseq(phybacz, phyfunz)
##########################################################################
# save the phyloseq object
saveRDS(phycomb, "./03-Routputfiles/D_Phyloobjects/D3A_comb_rel.rds")
##########################################################################
# filter the phyloseq object for ASVs that are at least present in 10% of the samples
otucomb <- otu_table(phycomb)
otu_present_count <- apply(otucomb, 1, function(x) sum(x > 0))
filtered_indices <- otu_present_count >= 0.1 * ncol(otucomb)
filtered_otucomb <- otucomb[filtered_indices, ]
phycomb_filtered <- prune_taxa(filtered_indices, phycomb)
##########################################################################
# save the filtered phyloseq object
saveRDS(phycomb_filtered, "./03-Routputfiles/D_Phyloobjects/D3A_comb_rel_filtered.rds")
##########################################################################
# generate otu tables for Cytoscape
otucomb <- data.frame(otu_table(phycomb_filtered))
otucomb <- data.frame("ASV" = rownames(otucomb), otucomb)
# save as txt
write.table(otucomb, file = "./03-Routputfiles/D_Phyloobjects/D3A_comb_rel_otu.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
# create tax table for conet
tax <- data.frame(tax_table(phycomb_filtered))
tax$genusnew <- ifelse(is.na(tax$isoGenus), tax$Genus, tax$isoGenus)
tax$ASV <- rownames(tax)
taxclean <- tax[, c("ASV", "genusnew", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "qseqid")]
# save as txt
write.table(tax, file = "./03-Routputfiles/D_Phyloobjects/D3A_comb_rel_tax.txt",
            sep = "\t", quote = FALSE)
# save as csv
write.csv(taxclean, file = "./03-Routputfiles/D_Phyloobjects/D3A_comb_rel_tax.csv",
          row.names = FALSE)

##########################################################################


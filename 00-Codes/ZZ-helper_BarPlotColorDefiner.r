##########################################################################
# load libraries
library(phyloseq)
library(ape)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(httpgd)
library(RColorBrewer)
library(cowplot)
##########################################################################
# load data
##########################################################################
# phylo object
phy <- readRDS("./03-Routputfiles/D_Phyloobjects/D1F_bac_abscnc_nc_ns_genus.rds")
phyrel <- transform_sample_counts(phy, function(x) x / sum(x))
##########################################################################
# get the asv table
otu <- data.frame(otu_table(phyrel))
# get the taxonomy table
tax <- data.frame(tax_table(phyrel))
# get the metadata
meta <- data.frame(sample_data(phyrel))
##########################################################################
# combine otu and taxonomy table
if (identical(rownames(otu), rownames(tax))) {
  comb <- cbind(tax, otu)
} else {
  stop("rownames of otu and tax are not identical")
}
##########################################################################
start <- ncol(tax) + 1
# make it long format
long <- data.frame(pivot_longer(comb, cols = c(start:ncol(comb)),
                                names_to = "sampleid", values_to = "relab"))
##########################################################################
# find the maximum relative abundance per genus
maxrelab <- long %>%
  group_by(Genus) %>%
  summarise(maxrelab = max(relab, na.rm = TRUE))
# sort the genera by maxrelab
maxrelab <- data.frame(maxrelab[order(maxrelab$maxrelab, decreasing = TRUE), ])
plot(maxrelab$maxrelab, type = "l", xlab = "Genera", ylab = "Max. relative abundance")

# get the number of genera with maxrelab > 0.05
sum(maxrelab$maxrelab > 0.05)

# save the genera with maxrelab > 0.05
owncolgen <- maxrelab$Genus[maxrelab$maxrelab > 0.05]

saveRDS(owncolgen, "./03-Routputfiles/Z_Helperfiles/owncolgenbac.rds")

##########################################################################
# FUNGAL DATA
##########################################################################

# load data
phy <- readRDS("./03-Routputfiles/D_Phyloobjects/D2D_fun_abs_nc_genus.rds")
phyrel <- transform_sample_counts(phy, function(x) x / sum(x))
##########################################################################
# get the asv table
otu <- data.frame(otu_table(phyrel))
# get the taxonomy table
tax <- data.frame(tax_table(phyrel))
# get the metadata
meta <- data.frame(sample_data(phyrel))
##########################################################################
# combine otu and taxonomy table
if (identical(rownames(otu), rownames(tax))) {
  comb <- cbind(tax, otu)
} else {
  stop("rownames of otu and tax are not identical")
}
##########################################################################
start <- ncol(tax) + 1
# make it long format
long <- data.frame(pivot_longer(comb, cols = all_of(start:ncol(comb)),
                                names_to = "sampleid", values_to = "relab"))
##########################################################################
# find the maximum relative abundance per genus
maxrelab <- long %>%
  group_by(Genus) %>%
  summarise(maxrelab = max(relab, na.rm = TRUE))
# sort the genera by maxrelab
maxrelab <- data.frame(maxrelab[order(maxrelab$maxrelab, decreasing = TRUE), ])
plot(maxrelab$maxrelab, type = "l", xlab = "Genera", ylab = "Max. relative abundance")

# get the number of genera with maxrelab > 0.05
sum(maxrelab$maxrelab > 0.1)

# save the genera with maxrelab > 0.05
owncolgen <- maxrelab$Genus[maxrelab$maxrelab > 0.1]

saveRDS(owncolgen, "./03-Routputfiles/Z_Helperfiles/owncolgenfun.rds")

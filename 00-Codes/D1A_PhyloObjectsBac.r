##########################################################################
# This script creates phyloseq objects from the raw data
# Part A: raw data phylo object
        # contains the raw data of the 16S rRNA gene sequencing, includes 
        # the otu table, taxonomy table, metadata and phylogenetic tree
        # phyloseq object is saved as phy_ferdl and phy_tea
        # does contain also the spike in bacteria Imtechella and Allobacillus
        # does contain the mock community samples and the negative Controls
##########################################################################
# load libraries
library(phyloseq)
library(tidyr)
library(tidyverse)
library(ape)
##########################################################################
# PART A: Create Raw Data Phyloseq Object
##########################################################################

##########################################################################
# META DATA ##############################################################
##########################################################################
# import meta data
metadata <- read.csv("./01-metadata/samplemeta.tsv", sep = "\t")
# add slice information to metadata
slicekey <- read.csv("./01-metadata/slice_day_key.txt", sep = "\t")
metadata <- plyr::join(metadata, slicekey, by = "slice", type = "left")
rownames(metadata) <- metadata[, 1]
##########################################################################
# add plate count information to metadata (based on PCA plates)
counts <- read.table("./02-rawdata/Plate_counts_txt.txt",
                     header = TRUE, sep = "\t")
#select only data from PCA
counts <- counts %>%
  filter(medium == "PCA")
counts <- plyr::join(counts, slicekey, type = "left")
#calculate the mean log_counts for each sampleID, medium and day
bmc <- counts %>%
  group_by(sampleID, medium, slice) %>%
  summarise(samplemean_logcounts = mean(log_counts, na.rm = TRUE)) %>%
  data.frame()
colnames(bmc)[1] <- "piece"
rownames(metadata) <- metadata[, 1]
metadata <- plyr::join(metadata, bmc, by = c("piece", "slice"))
# replace Inf with NA
metadata[metadata ==  Inf] <- NA
# unlog samplemean_logcounts
metadata$samplemean_counts <- 10^metadata$samplemean_logcounts
##########################################################################
# keep only samples with 16S in target
metadata <- metadata[metadata$target ==  "16S", ]
meta <- sample_data(metadata)
rownames(meta) <- meta$sampleid
##########################################################################

##########################################################################
# ASV TABLE ##############################################################
##########################################################################
owd <- getwd()
setwd("./05-qiime/HTS_16S")
##########################################################################
# import asv table
asv <- read.table("./04-asvtable/exported-table/feature-table.txt",
                  sep = "\t", check.names = FALSE, header = TRUE, row.names = 1)
otu <- otu_table(asv, taxa_are_rows = TRUE)
# show me the samples in meta that are not in otu
meta[!(rownames(meta) %in% sample_names(otu)), ]
# remove samples in meta that are not in otu
meta <- meta[rownames(meta) %in% sample_names(otu), ]
##########################################################################
# TAXONOMY TABLE
##########################################################################
taxonomy <- read.table("./06-taxonomy/01-exported-taxonomy-raw/taxonomy.tsv",
                       sep = "\t", check.names = FALSE,
                       header = TRUE, row.names = 1)
taxlist <- strsplit(as.character(taxonomy$Taxon), ";",
                    fixed = TRUE)
max_length <- max(sapply(taxlist, length))
taxlist <- lapply(taxlist, function(x) {
  length(x) <- max_length
    x
})
taxonomy_df <- as.data.frame(do.call(rbind, taxlist))
colnames(taxonomy_df) <- c("Kingdom", "Phylum", "Class", "Order",
                           "Family", "Genus", "Species")
# remove everything before "__" including "__" in the taxonomy_df
taxonomy_df <- data.frame(lapply(taxonomy_df, function(x) sub(".*__", "", x)))
# add rownames
rownames(taxonomy_df) <- rownames(taxonomy)
##########################################################################
# import blast results
matches <- read.table("./08-blastdb/comparison.blastout", sep = "\t",
                      header = FALSE)
colnames(matches) <- c("qseqid", "sseqid", "pident", "qcovs",
                       "length", "mismatch", "gapopen", "qstart",
                       "qend", "sstart", "send", "evalue", "bitscore")
# keep only matches with 100% identity
matches <- matches[matches$pident == 100, ]
# keep only matches with length > 300
matches <- matches[matches$length > 260, ]
# order by qseqid
matches <- matches[order(matches$qseqid), ]
# keep only unique matches
matches <- matches[!duplicated(matches$sseqid), ]
# add matching isolates to taxonomy_df
taxonomy_df$sseqid <- rownames(taxonomy)
taxonomy_df <- plyr::join(taxonomy_df, matches, by = "sseqid", type = "left")
rownames(taxonomy_df) <- taxonomy_df$sseqid
# replace NA by the last classified taxonomy level
taxonomy_df[, 1:7] <- t(apply(taxonomy_df[,1:7], 1,
                              function(x) zoo::na.locf(x, na.rm = FALSE)))
tax <- tax_table(as.matrix(taxonomy_df))
##########################################################################

##########################################################################
# TREE ###################################################################
##########################################################################
tree <- read.tree("./08-phylogenetictreemafft/exported-rooted-tree/tree.nwk")
##########################################################################
setwd(owd)
saveRDS(matches, "./03-Routputfiles/D1A_isolatematchesbac.rds")
##########################################################################

##########################################################################
# CREATE PHYLOSEQ OBJECT #################################################
##########################################################################
phy <- phyloseq(otu, tax, meta, tree)
# add reads per sample to the meta data
colsumsraw <- colSums(otu_table(phy))
meta <- sample_data(phy)
names <- rownames(meta)
# add colsumsraw to metafabs
# check if rownames of metafabs and names of colsumsraw are identical
identical(names, names(colsumsraw))
meta$colsumsraw <- colsumsraw
# replace meta data in phyloseq object
sample_data(phy) <- sample_data(meta)
# save phy separately for Ferdl and Tea
phy_ferdl <- subset_samples(phy, owner == "Ferdl")
# reduce sample names
sample_names(phy_ferdl) <- substr(sample_names(phy_ferdl), 1, 8)
# reduce names in sample_data sampleid
sample_data(phy_ferdl)$sampleid <- substr(sample_data(phy_ferdl)$sampleid, 1, 8)
# change - to _ in sample names
sample_names(phy_ferdl) <- gsub("-", "_", sample_names(phy_ferdl))
# replace - with _ in sample_data sampleid
sample_data(phy_ferdl)$sampleid <- gsub("-", "_", sample_data(phy_ferdl)$sampleid)
# replace - with _ in sample_data internID
sample_data(phy_ferdl)$internID <- gsub("-", "_", sample_data(phy_ferdl)$internID)
saveRDS(phy_ferdl, "./03-Routputfiles/D_Phyloobjects/D1A_bac_raw_allsamples_spike.rds")
# subset and save samples for Tea
#phy_tea <- subset_samples(phy, owner == "Tea")
#saveRDS(phy_tea, "./03-Routputfiles/Dbac_phyraw_tea.rds")
##########################################################################
##########################################################################
phy_ferdl <- subset_samples(phy_ferdl, piece != "mo" & piece != "NK")
# save phy_ferdl filtered
saveRDS(phy_ferdl, "./03-Routputfiles/D_Phyloobjects/D1A_bac_raw_nocontrols_spike.rds")


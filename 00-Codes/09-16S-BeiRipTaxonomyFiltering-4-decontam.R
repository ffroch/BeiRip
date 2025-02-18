# Extracted from decontam pipeline:
# https://benjjneb.github.io/decontam/vignettes/decontam_intro.html
###############################################################################
# load libraries
###############################################################################
library(decontam)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(phyloseq)
library(tidyverse)
###############################################################################
# load data
###############################################################################
# Metadata
metadata <- read.csv("./01-metadata/samplemeta.tsv", sep  =  "\t")
# remove samples with 18S in target
metadata <- metadata[metadata$target !=  "18S", ]
# make internID as rownames
row.names(metadata) <- metadata[, 1]
# just select relvant columns
relmeta <- metadata %>%
  select("samplecontrol", "libconc")
################################################################################
owd <- getwd()
setwd("./05-qiime/HTS_16S/")
################################################################################
# ASV-table
# Beware the "#" included in the header of the fwature table (#OTU ID).
# If you followed my Amplicon_sequencing_data_analysis_pipeline.sh
# it was already removed.
asvtable <- read.csv("./04-asvtable/02-taxfilt/exported-table-taxfilt-wobh/feature-table.txt",
                     sep = "\t", check.names  =  FALSE)
numsamples <- length(asvtable)
otumat <- as.matrix(asvtable[, 2:numsamples])
rownames(otumat) <- as.data.frame(asvtable)[, 1]
otu  <-  otu_table(otumat, taxa_are_rows  =  TRUE)
samplesdata <- sample_data(metadata)
#check if rownames of OTU and samplesdata are the same
identical(sample_names(otu), rownames(samplesdata))
phy <- phyloseq(otu, samplesdata)
################################################################################
# taxonomy
taxonomy <- read.table("./06-taxonomy/01-exported-taxonomy-raw/taxonomy.tsv",
                       sep = "\t", check.names  =  FALSE,
                       header = TRUE, row.names = 1)
################################################################################
setwd("./04-asvtable/decontam/")
################################################################################
# Inspect library sizes #
################################################################################
df <- as.data.frame(sample_data(phy))
df$LibrarySize <- sample_sums(phy)
df <- df[order(df$LibrarySize), ]
df$Index <- seq_len(nrow(df))
p <- ggplot(data = df, aes(x = Index, y = LibrarySize, color = samplecontrol)) +
  geom_point()
#save plot as tiff
ggsave("./plot_librarysize.tiff", plot = p, width = 10,
       height = 5, units = "cm", dpi = 300)
################################################################################
# Identify contaminants - prevalence # selected method
################################################################################
sample_data(phy)$is.neg <- sample_data(phy)$samplecontrol  ==  "control"
# The threshold value will identify as contaminants all sequences/ASVs
# there are more prevalent (than that value) in negative controls than
# in positive (real) samples.
# Function to use with high biomass samples: isContaminant
contamdf_prev <- isContaminant(phy, method = "prevalence",
                               neg = "is.neg", threshold = 0.5)
table(contamdf_prev$contaminant)
head(which(contamdf_prev$contaminant))
# LIST CONTAMINANTS
cont_prev <- subset(contamdf_prev, contamdf_prev$contaminant  ==  TRUE)
cont_tax <- merge(cont_prev, taxonomy, by = 0, all.x = TRUE)
write.table(cont_tax, "./contaminantsprev.tax.txt",
            col.names  =  NA, quote  =  FALSE, row.names  =  TRUE, sep  =  "\t")
# PLOT CONTAMINANTS
ps_pa <- transform_sample_counts(phy, function(abund) 1 * (abund > 0))
ps_pa_neg <- prune_samples(sample_data(ps_pa)$samplecontrol  ==  "control",
                           ps_pa)
ps_pa_pos <- prune_samples(sample_data(ps_pa)$samplecontrol  ==  "sample",
                           ps_pa)
da_pa <- data.frame(pa.pos = taxa_sums(ps_pa_pos),
                    pa.neg = taxa_sums(ps_pa_neg),
                    contaminant = contamdf_prev$contaminant)
p <- ggplot(data = da_pa, aes(x = pa.neg, y = pa.pos, color = contaminant)) +
  geom_point() +
  xlab("Prevalence (Negative Controls)") +
  ylab("Prevalence (True Samples)")
#save plot as tiff
ggsave("./plot_prevalence.tiff", plot = p, width = 10, height = 5,
       units = "cm", dpi = 300)
# REMOVE CONTAMINANTS FROM PHYLOSEQ OBJECT
phy_noncontam_prev <- prune_taxa(!contamdf_prev$contaminant, phy)
otu_nocontam_prev <- otu_table(phy_noncontam_prev)
# EXPORT FILTERED ASV TABLE
write.table(otu_nocontam_prev, "./feature-table.decontam.txt", col.names  =  NA,
            quote  =  FALSE, row.names  =  TRUE, sep  =  "\t")
###############################################################################
# print library versions used and save as txt file
session_info <- capture.output(sessionInfo())
sink("./09-16S-4-SessionInfo.txt")
cat(session_info, sep = "\n")
sink()
################################################################################
setwd(owd)
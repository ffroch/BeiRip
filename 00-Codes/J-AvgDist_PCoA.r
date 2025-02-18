library(phyloseq)
library(ape)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(vegan)

#############################################################################
# DATA IMPORT AND PREPARATION
#############################################################################
# import phyloseq object
phy <- readRDS("./03-Routputfiles/D_Phyloobjects/D1A_bac_raw_nocontrols_spike.rds")
# remove spike-in asvs
phy <- subset_taxa(phy, Genus != "Allobacillus" & Genus != "Imtechella")
meta <- data.frame(sample_data(phy))
taxbac <- data.frame(tax_table(phy))
#############################################################################
# import fungal data
phyfun <- readRDS("./03-Routputfiles/D_Phyloobjects/D2A_fun_raw_nocontrols.rds")
# remove samples with less than 20 reads
otufun <- otu_table(phyfun)
seq_depth <- colSums(otufun)
keep_asvs <- names(seq_depth[seq_depth >= 20])
phyfun <- prune_samples(keep_asvs, phyfun)
otufun <- otu_table(phyfun)
taxfun <- data.frame(tax_table(phyfun))
#############################################################################
# filter 16s samples from 18s samples
phy <- prune_samples(sample_names(phy) %in% sample_names(phyfun), phy)
otubac <- otu_table(phy)
#############################################################################
mindepth <- min(colSums(otubac))
matrix <- as.matrix(data.frame(t(otubac)))
if (file.exists("./03-Routputfiles/J_AvgDist_PCoA/bray.rds")) {
  bray <- readRDS("./03-Routputfiles/J_AvgDist_PCoA/bray.rds")
} else {
  set.seed(123)
  bray <- avgdist(matrix, dmethod = "bray", iterations = 1000, sample = mindepth)
  saveRDS(bray, "./03-Routputfiles/J_AvgDist_PCoA/bray.rds")
}
#############################################################################
mindepth <- min(colSums(otufun))
matrix <- as.matrix(data.frame(t(otufun)))
if (file.exists("./03-Routputfiles/J_AvgDist_PCoA/brayfun.rds")) {
  brayfun <- readRDS("./03-Routputfiles/J_AvgDist_PCoA/brayfun.rds")
} else {
  set.seed(123)
  brayfun <- avgdist(matrix, dmethod = "bray", iterations = 1000, sample = mindepth)
  saveRDS(brayfun, "./03-Routputfiles/J_AvgDist_PCoA/brayfun.rds")
}

# agglomerate data to species level considering isolate
#taxbac$Species <- paste0(taxbac$Species, "_", taxbac$qseqid)
#tax_table(phy) <- as.matrix(taxbac)
phyreds <- tax_glom(phy, "Species")
otubacreds <- otu_table(phyreds)
#taxfun$Species <- paste0(taxfun$Species, "_", taxfun$qseqid)
phyfunreds <- tax_glom(phyfun, "Species")
otufunreds <- otu_table(phyfunreds)
# aggregate data to genus level
phyredg <- tax_glom(phy, "Genus")
otubacredg <- otu_table(phyredg)
phyfunredg <- tax_glom(phyfun, "Genus")
otufunredg <- otu_table(phyfunredg)
#############################################################################
#PCoA
pcoa16s <- cmdscale(bray, eig = TRUE, k = 2)

# fit data using envfit
if (file.exists("./03-Routputfiles/J_AvgDist_PCoA/envit16SASV.rds")) {
  fit16sa <- readRDS("./03-Routputfiles/J_AvgDist_PCoA/envit16SASV.rds")
} else {
  set.seed(123)
  fit16sa <- envfit(pcoa16s, t(otubac), perm = 999)
  saveRDS(fit16sa, "./03-Routputfiles/J_AvgDist_PCoA/envit16SASV.rds")
}
vectors_16sa <- as.data.frame(fit16sa$vectors$arrows)
vectors_16sa$pvals <- fit16sa$vectors$pvals
vectors_16sa$r <- fit16sa$vectors$r
vectors_16sa_sig <- vectors_16sa[vectors_16sa$pvals < 0.05 & vectors_16sa$r > 0.3,]
vectors_16sa_sig[,1:2] <- vectors_16sa_sig[,1:2] * 0.5
print(vectors_16sa_sig)
vectors_16sa_sig$sseqid <- rownames(vectors_16sa_sig)
vectors_16sa_sig <- plyr::join(vectors_16sa_sig, taxbac, by = "sseqid", type = "left")

if (file.exists("./03-Routputfiles/J_AvgDist_PCoA/envit16Sspecies.rds")) {
  fit16ss <- readRDS("./03-Routputfiles/J_AvgDist_PCoA/envit16Sspecies.rds")
} else {
  set.seed(123)
  fit16ss <- envfit(pcoa16s, t(otubacreds), perm = 999)
  saveRDS(fit16ss, "./03-Routputfiles/J_AvgDist_PCoA/envit16Sspecies.rds")
}
vectors_16ss <- as.data.frame(fit16ss$vectors$arrows)
vectors_16ss$pvals <- fit16ss$vectors$pvals
vectors_16ss$r <- fit16ss$vectors$r
vectors_16ss_sig <- vectors_16ss[vectors_16ss$pvals < 0.05 & vectors_16ss$r > 0.3,]
vectors_16ss_sig[,1:2] <- vectors_16ss_sig[,1:2] * 0.6
print(vectors_16ss_sig)
vectors_16ss_sig$sseqid <- rownames(vectors_16ss_sig)
vectors_16ss_sig <- plyr::join(vectors_16ss_sig, taxbac, by = "sseqid", type = "left")
#remove vectors with Bacteria in Species
vectors_16ss_sig <- vectors_16ss_sig[!grepl("Bacteria", vectors_16ss_sig$Species),]

if (file.exists("./03-Routputfiles/J_AvgDist_PCoA/envit16Sgenus.rds")) {
  fit16sg <- readRDS("./03-Routputfiles/J_AvgDist_PCoA/envit16Sgenus.rds")
} else {
  set.seed(123)
  fit16sg <- envfit(pcoa16s, t(otubacredg), perm = 999)
  saveRDS(fit16sg, "./03-Routputfiles/J_AvgDist_PCoA/envit16Sgenus.rds")
}
vectors_16sg <- as.data.frame(fit16sg$vectors$arrows)
vectors_16sg$pvals <- fit16sg$vectors$pvals
vectors_16sg$r <- fit16sg$vectors$r
vectors_16sg_sig <- vectors_16sg[vectors_16sg$pvals < 0.05 & vectors_16sg$r > 0.3,]
vectors_16sg_sig[,1:2] <- vectors_16sg_sig[,1:2] * 0.3
print(vectors_16sg_sig)
vectors_16sg_sig$sseqid <- rownames(vectors_16sg_sig)
vectors_16sg_sig <- plyr::join(vectors_16sg_sig, taxbac, by = "sseqid", type = "left")

if (file.exists("./03-Routputfiles/J_AvgDist_PCoA/envit18SASV.rds")) {
  fit18sa <- readRDS("./03-Routputfiles/J_AvgDist_PCoA/envit18SASV.rds")
} else {
  set.seed(123)
  fit18sa <- envfit(pcoa16s, t(otufun), perm = 999)
  saveRDS(fit18sa, "./03-Routputfiles/J_AvgDist_PCoA/envit18SASV.rds")
}
vectors_18sa <- as.data.frame(fit18sa$vectors$arrows)
vectors_18sa$pvals <- fit18sa$vectors$pvals
vectors_18sa$r <- fit18sa$vectors$r
vectors_18sa_sig <- vectors_18sa[vectors_18sa$pvals < 0.05 & vectors_18sa$r > 0.1,]
vectors_18sa_sig[,1:2] <- vectors_18sa_sig[,1:2] * 0.5
print(vectors_18sa_sig)
vectors_18sa_sig$sseqid <- rownames(vectors_18sa_sig)
vectors_18sa_sig <- plyr::join(vectors_18sa_sig, taxfun, by = "sseqid", type = "left")

if (file.exists("./03-Routputfiles/J_AvgDist_PCoA/envit18Sspecies.rds")) {
  fit18ss <- readRDS("./03-Routputfiles/J_AvgDist_PCoA/envit18Sspecies.rds")
} else {
  set.seed(123)
  fit18ss <- envfit(pcoa16s, t(otufunreds), perm = 999)
  saveRDS(fit18ss, "./03-Routputfiles/J_AvgDist_PCoA/envit18Sspecies.rds")
}
vectors_18ss <- as.data.frame(fit18ss$vectors$arrows)
vectors_18ss$pvals <- fit18ss$vectors$pvals
vectors_18ss$r <- fit18ss$vectors$r
vectors_18ss_sig <- vectors_18ss[vectors_18ss$pvals < 0.05 & vectors_18ss$r > 0.1,]
vectors_18ss_sig[,1:2] <- vectors_18ss_sig[,1:2] * 0.6
print(vectors_18ss_sig)
vectors_18ss_sig$sseqid <- rownames(vectors_18ss_sig)
vectors_18ss_sig <- plyr::join(vectors_18ss_sig, taxfun, by = "sseqid", type = "left")

if (file.exists("./03-Routputfiles/J_AvgDist_PCoA/envit18Sgenus.rds")) {
  fit18sg <- readRDS("./03-Routputfiles/J_AvgDist_PCoA/envit18Sgenus.rds")
} else {
  set.seed(123)
  fit18sg <- envfit(pcoa16s, t(otufunredg), perm = 999)
  saveRDS(fit18sg, "./03-Routputfiles/J_AvgDist_PCoA/envit18Sgenus.rds")
}
vectors_18sg <- as.data.frame(fit18sg$vectors$arrows)
vectors_18sg$pvals <- fit18sg$vectors$pvals
vectors_18sg$r <- fit18sg$vectors$r
vectors_18sg_sig <- vectors_18sg[vectors_18sg$pvals < 0.05 & vectors_18sg$r > 0.1,]
vectors_18sg_sig[,1:2] <- vectors_18sg_sig[,1:2] * 0.3
print(vectors_18sg_sig)
vectors_18sg_sig$sseqid <- rownames(vectors_18sg_sig)
vectors_18sg_sig <- plyr::join(vectors_18sg_sig, taxfun, by = "sseqid", type = "left")

# quick and dirty correlations
cor16s <- as.data.frame(cor(t(otubac), pcoa16s$points[,1:2]))
cor18s <- as.data.frame(cor(t(otufun), pcoa16s$points[,1:2]))
names(cor16s) <- c("PC1", "PC2")
cor16s$sseqid <- rownames(cor16s)
names(cor18s) <- c("PC1", "PC2")
cor18s$sseqid <- rownames(cor18s)
# filter for correlation values > 0.5
threshold <- 0.2
cor16sfil <- cor16s[abs(cor16s$PC1) > threshold*2.5 | abs(cor16s$PC2) > threshold*2.5,]
# remove NA values
cor16sfil <- cor16sfil[!is.na(cor16sfil$PC1),]
cor18sfil <- cor18s[abs(cor18s$PC1) > threshold | abs(cor18s$PC2) > threshold,]
# remove NA values
cor18sfil <- cor18sfil[!is.na(cor18sfil$PC1),]
cor16sfil <- plyr::join(cor16sfil, taxbac, by = "sseqid", type = "left")
cor18sfil <- plyr::join(cor18sfil, taxfun, by = "sseqid", type = "left")

pcoa_points <- as.data.frame(pcoa16s$points)
pcoa_points$sampleid <- rownames(pcoa_points)
pcoa_points <- plyr::join(pcoa_points, meta, by = "sampleid", type = "left")

pcoa_plot <- ggplot(pcoa_points, aes(x = V1, y = V2)) +
  geom_point(aes(color = factor(batch), size = days)) +
  theme_bw() +
  geom_segment(data = cor16sfil, aes(x = 0, y = 0, xend = PC1, yend = PC2), color = "red") +
  geom_label(data = cor16sfil, aes(x = PC1, y = PC2, label = Genus), color = "red", nudge_x = 0.02, nudge_y = 0.02, size = 2.5, fill = NA, label.size = NA) +
  geom_segment(data = vectors_16ss_sig, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), color = "orange") +
  geom_label(data = vectors_16ss_sig, aes(x = Dim1, y = Dim2, label = Genus), color = "orange", nudge_x = 0.02, nudge_y = 0.02, size = 2.5, fill = NA, label.size = NA) +
  geom_segment(data = cor18sfil, aes(x = 0, y = 0, xend = PC1, yend = PC2), color = "blue") +
  geom_label(data = cor18sfil, aes(x = PC1, y = PC2, label = Genus), color = "blue", nudge_x = 0.02, nudge_y = 0.02, size = 2.5, fill = NA, label.size = NA) +
  geom_segment(data = vectors_18ss_sig, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), color = "green") +
  geom_label(data = vectors_18ss_sig, aes(x = Dim1, y = Dim2, label = Genus), color = "green", nudge_x = 0.02, nudge_y = 0.02, size = 2.5, fill = NA, label.size = NA)

pcoa_plot

library(ggrepel)
set.seed(1234)
pcoa_plot <- ggplot(pcoa_points, aes(x = V1, y = V2)) +
  geom_point(aes(color = factor(batch), size = days), alpha = 0.7) +
  theme_bw() +
  geom_segment(data = vectors_16ss_sig, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), color = "darkred") +
  geom_label_repel(data = vectors_16ss_sig, aes(x = Dim1, y = Dim2, label = Species), color = "darkred", nudge_x = 0.0, nudge_y = 0.0, size = 2.5, fill = NA, label.size = NA, box.padding = 0) +
  geom_segment(data = vectors_18ss_sig, aes(x = 0, y = 0, xend = Dim1, yend = Dim2), color = "#585858") +
  geom_label_repel(data = vectors_18ss_sig, aes(x = Dim1, y = Dim2, label = Species), color = "#585858", nudge_x = 0.02, nudge_y = 0.02, size = 2.5, fill = NA, label.size = NA) +
  scale_color_manual(values = c("#9E001B", "#002D41", "#388E3C"), name = "batch")+
  xlab("PCoA1") +
  ylab("PCoA2")

pcoa_plot

ggsave("./03-Routputfiles/J_AvgDist_PCoA/J_PCoA_plot.svg", pcoa_plot, width = 15, height = 15, units = "cm", dpi = 300)
ggsave("./03-Routputfiles/J_AvgDist_PCoA/J_PCoA_plot.png", pcoa_plot, width = 15, height = 15, units = "cm", dpi = 300)
saveRDS(pcoa_plot, "./03-Routputfiles/J_AvgDist_PCoA/J_PCoA_plot.rds")

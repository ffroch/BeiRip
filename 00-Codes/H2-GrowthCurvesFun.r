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
##########################################################################
# PART 1: ASV level
##########################################################################
##########################################################################
# load data
##########################################################################
# phylo object
phyabs <- readRDS("./03-Routputfiles/D_Phyloobjects/D2C_fun_abs_nocontrols.rds")
##########################################################################
tt <- data.frame(tax_table(phyabs))
colnames(tt)[which(colnames(tt) == "sseqid")] <- "ASV"
tt$label <- paste0(tt$ASV, "\n", tt$Genus, "\n", tt$qseqid)


# log10 transform the abundance data in phyabs
otu <- otu_table(phyabs)
otu <- log10(otu + 1)
totu <- t(otu)
# get metadata
meta <- as.data.frame(sample_data(phyabs))
# join metadata with otu
if(identical(rownames(meta), rownames(totu))) {
  comb <- cbind(meta, as.data.frame(totu))
} else {
  stop("rownames of metadata and otu are not identical")
}
##########################################################################
# Part 1.1: Taxaplots
##########################################################################
# make a sequence starting with at start with a step of 6 until 100+start
start <- ncol(meta) + 1
# sort the columns start:ncol(comb) by the sum of the columns
comb1 <- comb[, c(1:(start - 1), order(colSums(comb[, start:ncol(comb)]),
              decreasing = TRUE) + start - 1)]
sequence <- seq(start, start + 24, 12)
# prepare color set
colsel <- c("#9E001B", "#3DB2C8", "#03546C")
# plot
for (i in sequence) {
  # generate asubset dataset
  sub <- comb1[, c(1:(start - 1), i:(i + 11))]
  # make it long
  sublong <- data.frame(pivot_longer(sub,
                                     cols =  start:(start + 11),
                                     names_to = "ASV",
                                     values_to = "counts"))
  sublong <- left_join(sublong, tt, by = c("ASV"))
  # replace 0 in sublong$counts with -Inf
  sublong$counts[sublong$counts == 0] <- -Inf
  # plot
  p <- ggplot(sublong) +
    geom_point(aes(sublong$days, sublong$counts, col = factor(sublong$batch))) +
    stat_smooth(method = "loess", se = FALSE, aes(sublong$days, sublong$counts, group = factor(sublong$batch), col = factor(sublong$batch)), alpha = 0.2, linewidth = 0.1, linetype = 2) +
    stat_smooth(method = "loess", se = TRUE, alpha = 0.2, group = 1, aes(sublong$days, sublong$counts), col = "black") +
    theme_bw() +
    theme(legend.position = "bottom", 
          strip.background = element_rect(fill = "white")) +
    guides(color = guide_legend(nrow = 1, title = "batch")) +
    ylab("log10(gene copies)") +
    ylim(0, 5) +
    scale_color_manual(values = colsel) +
    scale_fill_manual(values = colsel) +
    facet_wrap(~label) +
    xlab("days")
  print(p)
  ggsave(paste("./03-Routputfiles/H_GrowthCurves/H2_GrowthCurvesASVFun",
               i - (start - 1), "-", i - (start - 12), ".svg", sep = ""),
         width = 25, height = 25, units = "cm", dpi = 300)
  ggsave(paste("./03-Routputfiles/H_GrowthCurves/H2_GrowthCurvesASVFun",
               i - (start - 1), "-", i - (start - 12), ".png", sep = ""),
         width = 25, height = 25, units = "cm", dpi = 300)
}
##########################################################################
##########################################################################
# Part 1.2: Heatmap
##########################################################################
# normalized absolute abundance plots
sequence <- seq(start, start + 24, 1)
complete <- data.frame()
for (i in sequence) {
  sub <- comb1[, c(1:(start - 1), i)]
  sublong <- data.frame(pivot_longer(sub,
                                     cols =  start,
                                     names_to = "ASV",
                                     values_to = "counts"))
  # calculate mean per day but take 0 as NA
  sublong$counts[sublong$counts == 0] <- NA
  sublongmean <- data.frame(aggregate(sublong$counts, by = list(sublong$days),
                                      FUN = function(x) mean(x, na.rm = TRUE)))

  minx <- min(sublongmean$x, na.rm = TRUE)
  sublongmean$norm <- (sublongmean$x-minx)/(max(sublongmean$x, na.rm = TRUE)-minx)
  sublongmean$ASV <- colnames(comb1)[i]
  # get column index of ASV in tt
  index <- which(colnames(tt) == "ASV")
  sublongmean <- plyr::join(sublongmean, tt[, c(1:7, index)], by = "ASV")
  complete <- rbind(complete, sublongmean)
}
colnames(complete) <- c("days", "mean", "norm", "ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
complete$labeln <- paste0(substr(complete$ASV,2,6), "-", complete$Order, "-", complete$Family,
                        "-", complete$Species)
# keep only norm of 1
compfil <- complete[!is.na(complete$norm), ]
compfil <- compfil[compfil$norm == 1, ]
# order comfil by days decreasing
compfil <- compfil[order(-compfil$days), ]
# add a column called index 
compfil$index <- 1:nrow(compfil)
compfil$fil <- ifelse(compfil$mean < 2, 1, 0)
# keep only ASV and index
compfil <- compfil[, c("ASV", "index", "fil")]
# join compfil with complete
completen <- plyr::join(compfil, complete, by = "ASV", type = "left")
compred <- completen[completen$fil == 0, ]
################################################################################
# additional ordering
compredn <- compred %>%
  select(labeln, days, norm) %>%
  spread(days, norm)
# give me the first column index of each row that is above 0.75
compredn$first <- apply(compredn[, 2:ncol(compredn)], 1, function(x) which(x > 0.75)[1])
compredn <- compredn %>%
  select(labeln, first)

compredadd <- plyr::join(compred, compredn, by = "labeln")
compredadd <- compredadd[order(-compredadd$first, compredadd$index), ]
compredadd$index2 <- 1:nrow(compredadd)
ggplot(compredadd, aes(x = factor(days), y = reorder(labeln, index2), fill = norm)) +
  geom_tile() +
  scale_fill_gradient2(low = "#1d4877", mid = "#fbb021", high = "#ee3e32", 
                       midpoint = 0.5, na.value = "white", name = "") +
  theme_bw() +
  labs(x = "days", y = "ASV")


# save plot as tiff
ggsave("./03-Routputfiles/H_GrowthCurves/H2_TopNormHeatMapASVFun.svg", width = 25,
       height = 25, units = "cm", dpi = 300)
ggsave("./03-Routputfiles/H_GrowthCurves/H2_TopNormHeatMapASVFun.png", width = 25,
       height = 25, units = "cm", dpi = 300)
##########################################################################
##########################################################################
# PART 2: Genus level
##########################################################################
##########################################################################
# load data
##########################################################################
# phylo object
phygen <- readRDS("./03-Routputfiles/D_Phyloobjects/D2D_fun_abs_nc_genus.rds")
##########################################################################
# log10 transform the abundance data in phygen
otu <- otu_table(phygen)
otu <- log10(otu + 1)
totu <- t(otu)
# get metadata
meta <- as.data.frame(sample_data(phygen))
# calculate the detection limit per sample
meta$detlim <- log10((1 - (1 - 0.95)^(1/meta$colsumsraw))*meta$samplemean_counts)
# join metadata with otu
if(identical(rownames(meta), rownames(totu))) {
  comb <- cbind(meta, as.data.frame(totu))
} else {
  stop("rownames of metadata and otu are not identical")
}
# get taxonomy
tt <- data.frame(tax_table(phyabs))
colnames(tt)[which(colnames(tt) == "sseqid")] <- "ASV"
tt$label <- paste0(tt$ASV, "\n", tt$Genus, "\n", tt$qseqid)
##########################################################################
# Part 2.1: Taxaplots
##########################################################################
# make a sequence starting with at start with a step of 6 until 100+start
start <- ncol(meta) + 1
# sort the columns start:ncol(comb) by the sum of the columns
comb1 <- comb[, c(1:(start - 1), order(colSums(comb[, start:ncol(comb)]),
              decreasing = TRUE) + start - 1)]
sequence <- seq(start, start + 18, 12)
# prepare color set
colsel <- c("#D32F2F", "#303F9F", "#00796B")
# plot
for (i in sequence) {
  # generate asubset dataset
  sub <- comb1[, c(1:(start - 1), i:(i + 11))]
  # make it long
  sublong <- data.frame(pivot_longer(sub,
                                     cols =  start:(start + 11),
                                     names_to = "ASV",
                                     values_to = "counts"))
  sublong <- left_join(sublong, tt, by = c("ASV"))
  # replace 0 in sublong$counts with -Inf
  sublong$counts[sublong$counts == 0] <- -Inf
  # plot
  p <- ggplot(sublong) +
    geom_point(aes(sublong$days, sublong$counts, col = factor(sublong$batch))) +
    stat_smooth(method = "loess", se = FALSE, aes(sublong$days, sublong$counts, group = factor(sublong$batch), col = factor(sublong$batch)), alpha = 0.2, linewidth = 0.1, linetype = 2) +
    stat_smooth(method = "loess", se = TRUE, alpha = 0.2, group = 1, aes(sublong$days, sublong$counts), col = "black") +
    theme_bw() +
    theme(legend.position = "bottom", 
          strip.background = element_rect(fill = "white")) +
    guides(color = guide_legend(nrow = 1, title = "batch")) +
    ylab("log10(gene copies)") +
    ylim(0, 5) +
    scale_color_manual(values = colsel) +
    scale_fill_manual(values = colsel) +
    facet_wrap(~label) +
    xlab("days")
  print(p)
  ggsave(paste("./03-Routputfiles/H_GrowthCurves/H2_GrowthCurvesGenusFun_",
               i - (start - 1), "-", i - (start - 12), ".svg", sep = ""),
         width = 25, height = 25, units = "cm", dpi = 300)
  ggsave(paste("./03-Routputfiles/H_GrowthCurves/H2_GrowthCurvesGenusFun_",
               i - (start - 1), "-", i - (start - 12), ".png", sep = ""),
         width = 25, height = 25, units = "cm", dpi = 300)
}
##########################################################################
# Part 2.2: Heatmap
##########################################################################
# normalized absolute abundance plots
sequence <- seq(start, start + 18, 1)
complete <- data.frame()
for (i in sequence) {
  sub <- comb1[, c(1:(start - 1), i)]
  sublong <- data.frame(pivot_longer(sub,
                                     cols =  start,
                                     names_to = "ASV",
                                     values_to = "counts"))
  # calculate mean per day but take 0 as NA
  sublong$counts[sublong$counts == 0] <- NA
  sublongmean <- data.frame(aggregate(sublong$counts, by = list(sublong$days),
                                      FUN = function(x) mean(x, na.rm = TRUE)))

  minx <- min(sublongmean$x, na.rm = TRUE)
  sublongmean$norm <- (sublongmean$x-minx)/(max(sublongmean$x, na.rm = TRUE)-minx)
  sublongmean$ASV <- colnames(comb1)[i]
  # get column index of ASV in tt
  index <- which(colnames(tt) == "ASV")
  sublongmean <- plyr::join(sublongmean, tt[, c(1:7, index)], by = "ASV")
  complete <- rbind(complete, sublongmean)
}
colnames(complete) <- c("days", "mean", "norm", "ASV", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
complete$labeln <- paste0(substr(complete$ASV,2,6), "-", complete$Order, "-", complete$Family,
                        "-", complete$Genus)
# keep only norm of 1
compfil <- complete[!is.na(complete$norm), ]
compfil <- compfil[compfil$norm == 1, ]
# order comfil by days decreasing
compfil <- compfil[order(-compfil$days), ]
# add a column called index 
compfil$index <- 1:nrow(compfil)
compfil$fil <- ifelse(compfil$mean < 2, 1, 0)
# keep only ASV and index
compfil <- compfil[, c("ASV", "index", "fil")]
# join compfil with complete
completen <- plyr::join(compfil, complete, by = "ASV", type = "left")
compred <- completen[completen$fil == 0, ]
# replace NA with 0
#compred$norm[is.na(compred$norm)] <- 0
# adapt label
compred$labeln <- substring(compred$labeln, 7, 100)
################################################################################
# additional ordering
compredn <- compred %>%
  select(labeln, days, norm) %>%
  spread(days, norm)
# give me the first column index of each row that is above 0.75
compredn$first <- apply(compredn[, 2:ncol(compredn)], 1, function(x) which(x > 0.75)[1])
compredn <- compredn %>%
  select(labeln, first)

compredadd <- plyr::join(compred, compredn, by = "labeln")
compredadd <- compredadd[order(-compredadd$first, compredadd$index), ]
compredadd$index2 <- 1:nrow(compredadd)
ggplot(compredadd, aes(x = factor(days), y = reorder(labeln, index2), fill = norm)) +
  geom_tile() +
  scale_fill_gradient2(low = "#1d4877", mid = "#fbb021", high = "#ee3e32", 
                       midpoint = 0.5, na.value = "white", name = "") +
  theme_bw() +
  labs(x = "days", y = "Taxon")

ggsave("./03-Routputfiles/H_GrowthCurves/H2_TopNormHeatMapGenFun.svg", width = 25,
       height = 25, units = "cm", dpi = 300)
ggsave("./03-Routputfiles/H_GrowthCurves/H2_TopNormHeatMapGenFun.png", width = 25,
       height = 25, units = "cm", dpi = 300)

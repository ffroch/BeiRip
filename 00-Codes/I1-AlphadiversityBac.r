# I will come directly to the facit right now. There is a problem with the
# absolute abundance data and the iNEXT function. But since the absolute
# abundance data has the same proportions as the raw data anyway, we will 
# use the copy number corrected raw data. Since the function need integers
# the numbers are rounded to the next integer. 
# still i want to have a comparison between the copy number corrected data,
# the raw data and their relative abundance data.
##########################################################################
library(phyloseq)
library(ape)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(httpgd)
library(iNEXT)
library(cowplot)
##########################################################################
# define functions from iNEXT package
Chat.Ind <- function(x, m){
  x <- x[x>0]
  n <- sum(x)
  f1 <- sum(x == 1)
  f2 <- sum(x == 2)
  f0.hat <- ifelse(f2 == 0, (n - 1) / n * f1 * (f1 - 1) / 2, (n - 1) / n * f1 ^ 2/ 2 / f2)  #estimation of unseen species via Chao1
  A <- ifelse(f1>0, n*f0.hat/(n*f0.hat+f1), 1)
  Sub <- function(m){
    #if(m < n) out <- 1-sum(x / n * exp(lchoose(n - x, m)-lchoose(n - 1, m)))
    if(m < n) {
      xx <- x[(n-x)>=m]
      out <- 1-sum(xx / n * exp(lgamma(n-xx+1)-lgamma(n-xx-m+1)-lgamma(n)+lgamma(n-m)))
    }
    if(m == n) out <- 1-f1/n*A
    if(m > n) out <- 1-f1/n*A^(m-n+1)
    out
  }
  sapply(m, Sub)		
}
##########################################################################
# DATA IMPORT AND PREPARATION
phyraw <- readRDS("./03-Routputfiles/D_Phyloobjects/D1A_bac_raw_nocontrols_spike.rds")
phycnc <- readRDS("./03-Routputfiles/D_Phyloobjects/D1C_bac_rawcnc_nocontrols_spike.rds")
# remove spike in cells from both phyloseq objects
phyraw <- subset_taxa(phyraw, Genus != "Imtechella" & Genus != "Allobacillus")
phycnc <- subset_taxa(phycnc, Genus != "Imtechella" & Genus != "Allobacillus")
# get the otu tables
oturaw <- data.frame(otu_table(phyraw))
oturaw_t <- t(oturaw)
otucnc <- data.frame(otu_table(phycnc))
otucnc_t <- t(otucnc)
# correct otucnc to the same total number of reads as oturaw
otucnc <- sweep(otucnc, 2, colSums(oturaw)/colSums(otucnc), "*")
# the numbers in otucnc will be rounded to the next integer, so it is possible 
# that asvs with values <0.5 will be lost
# so i will use ceiling function to round up instead of round
otucncc <- ceiling(otucnc)
# make a new phyloseq object with the copy number corrected data
otu_table(phycnc) <- otu_table(as.matrix(otucncc), taxa_are_rows = TRUE)
# calculate coverage
cov_per_sampleraw <- apply(oturaw_t, 1, function(x) Chat.Ind(x, sum(x)))
cov_dfraw <- data.frame(sample = rownames(oturaw_t), coverage = cov_per_sampleraw)
cov_dfraw
cov_per_samplecnc <- apply(otucncc, 1, function(x) Chat.Ind(x, sum(x)))
cov_dfcnc <- data.frame(sample = rownames(otucncc), coverage = cov_per_samplecnc)
cov_dfcnc
# estimateD
set.seed(123)
bycoverageraw <- estimateD(oturaw, base = "coverage", level = 0.995, nboot = 50)
set.seed(123)
bycoveragecnc <- estimateD(otucnc, base = "coverage", level = 0.995, nboot = 50)
# combine the data
bycoverageraw$data <- "raw"
bycoveragecnc$data <- "cnc"
hilldiv <- rbind(bycoverageraw, bycoveragecnc)
# get meta data
meta <- data.frame(sample_data(phycnc))
# add meta
hilldiv_m <- dplyr::left_join(hilldiv, meta, by = c("Assemblage" = "sampleid"))
##########################################################################
hilldiv_m$Order.q <- factor(hilldiv_m$Order.q,
                            levels = c(0, 1, 2),
                            labels = c("species richness", "Hill-Shannon",
                                       "Hill-Simpson"))
hilldiv_m$data <- factor(hilldiv_m$data, levels = c("raw", "cnc"),
                         labels = c("raw reads", "copy number corrected reads"))
# plotting
# plot alpha diversity over time
#hilldiv_m <- hilldiv_m[hilldiv_m$Order.q != "species richness", ]
#hilldiv_m <- droplevels(hilldiv_m)
p1 <- ggplot(hilldiv_m, aes(x = factor(days), y = qD)) +
  geom_boxplot(outlier.shape = NA, aes(color = factor(batch), fill = factor(batch))) +
  geom_boxplot(outlier.shape = NA,  color = alpha("black", 1), alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.5, aes(color = factor(batch))) +
  facet_grid(Order.q ~ data, scales = "free_y") +
  scale_color_manual(values = c("#9E001B", "#002D41", "#388E3C"), name = "batch") +
  scale_fill_manual(values = c("#9E001B", "#002D41", "#388E3C"), name = "batch") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "Days", y = "qD")
print(p1)
# save as tiff 25x25 cm
ggsave("./03-Routputfiles/I_Alphadiv/I1-HilldiversitySliceBac_boxbox.svg", p1, width = 25,
       height = 25, units = "cm", dpi = 300)
ggsave("./03-Routputfiles/I_Alphadiv/I1-HilldiversitySliceBac_boxbox.png", p1, width = 25,
       height = 25, units = "cm", dpi = 300)
# alternative plotting
hilldiv_m <- hilldiv_m %>%
  mutate(
    position = case_when(
      batch == 1 ~ as.numeric(as.factor(days)) - 0.25,
      batch == 2 ~ as.numeric(as.factor(days)),
      batch == 3 ~ as.numeric(as.factor(days)) + 0.25
    )
  )

p2 <- ggplot(hilldiv_m, aes(x = factor(days), y = qD)) +
  geom_boxplot(outlier.shape = NA,  color = alpha("black", 1), alpha = 0.1) +
  geom_point(alpha = 0.7, aes(color = factor(batch), x = position)) +
  facet_grid(Order.q ~ data, scales = "free_y") +
  scale_color_manual(values = c("#9E001B", "#002D41", "#388E3C"), name = "batch") +
  theme_bw() +
  xlab("days") +
  theme(strip.background = element_rect(fill = "white"))
print(p2)
ggsave("./03-Routputfiles/I_Alphadiv/I1-HilldiversitySliceBac_boxjit.svg", p2, width = 25,
       height = 25, units = "cm", dpi = 300)#
ggsave("./03-Routputfiles/I_Alphadiv/I1-HilldiversitySliceBac_boxjit.png", p2, width = 25,
       height = 25, units = "cm", dpi = 300)

##########################################################################
# Hill Shannon only
hilldiv_f <- hilldiv_m[hilldiv_m$Order.q == "Hill-Shannon", ]
hilldiv_f <- hilldiv_f[hilldiv_f$data == "copy number corrected reads", ]
hilldiv_f <- hilldiv_f %>%
  mutate(
    position = case_when(
      batch == 1 ~ as.numeric(as.factor(days)) - 0.25,
      batch == 2 ~ as.numeric(as.factor(days)),
      batch == 3 ~ as.numeric(as.factor(days)) + 0.25
    )
  )

# plotting
# plot alpha diversity over time
p2 <- ggplot(hilldiv_f, aes(x = factor(days), y = qD)) +
  geom_boxplot(outlier.shape = NA,  color = alpha("black", 1), alpha = 0.1) +
  geom_point(alpha = 0.7, aes(color = factor(batch), x = position)) +
  scale_color_manual(values = c("#9E001B", "#002D41", "#388E3C"), name = "batch") +
  theme_bw() +
  xlab("days") +
  ylab("Hill-Shannon diversity")
print(p2)
ggsave("./03-Routputfiles/I_Alphadiv/I1-HilldiversitySliceBac_boxjit_shannononly.svg", p2, width = 25,
       height = 25, units = "cm", dpi = 300)#
ggsave("./03-Routputfiles/I_Alphadiv/I1-HilldiversitySliceBac_boxjit_shannononly.png", p2, width = 25,
       height = 25, units = "cm", dpi = 300)
saveRDS(p2, "./03-Routputfiles/I_Alphadiv/I1-HilldiversitySliceBac_boxjit_shannononly.rds")
################################################################################
# comparison alpha diversity over batch
# create new group for comparison
kwtest <- list()
dunntest <- list()
plotlist <- list()
for (j in unique(hilldiv_m$data)) {
  hilldiv_mfil <- hilldiv_m[hilldiv_m$data == j, ]
  for (i in unique(hilldiv_mfil$Order.q)) {
  key <- paste(j, i, sep = "_")
  kwtest[[key]] <- kruskal.test(qD ~ batch, data = hilldiv_mfil[hilldiv_mfil$Order.q == i, ])
  dunntest[[key]] <- dunn.test::dunn.test(hilldiv_mfil$qD[hilldiv_mfil$Order.q == i],
                                       hilldiv_mfil$batch[hilldiv_mfil$Order.q == i],
                                       method = "bonferroni")
  annotations <- data.frame(comparison = dunntest[[key]]$comparisons,
                            p.value = dunntest[[key]]$P.adjusted,
                            label = ifelse(dunntest[[key]]$P.adjusted < 0.025, "*", ""),
                            pkw = rep(kwtest[[key]]$p.value, length(dunntest[[key]]$P.adjusted))) %>%
                            separate(comparison, into = c("BatchA", "BatchB"), sep = " - ")
  annotations$show <- ifelse(annotations$p.value < 0.025 & annotations$pkw < 0.05, TRUE, FALSE)
  ymax <- max(hilldiv_mfil[hilldiv_mfil$Order.q == i, ]$qD)
  annotations$y_position <- c(ymax * 1.05, ymax * 1.1, ymax * 1.15)
  annotations$formatted_pvalue <- sprintf("%.3f", annotations$p.value)
  plotlist[[key]] <- ggplot(hilldiv_mfil[hilldiv_mfil$Order.q == i, ], aes(x = factor(batch), y = qD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, alpha = 0.5, aes(color = days)) +
  facet_grid(Order.q ~ data, scales = "free_y") +
  scale_color_continuous(low = "#3DB2C8", high = "#9E001B") + 
  theme_bw() +
  geom_segment(data = annotations,
               aes(x = BatchA, xend = BatchB,
                   y = y_position,
                   yend = y_position),
                   color = ifelse(annotations$show == TRUE, "black", "red"),
                   alpha = ifelse(annotations$show == TRUE, 1, 0)) +
  annotate("text", 
           x = (as.numeric(annotations$BatchA) + 
                as.numeric(annotations$BatchB)) / 2,
           y = annotations$y_position, 
           label = annotations$formatted_pvalue, 
           vjust = -0.1, color = ifelse(annotations$show == TRUE, "black", "red"),
                   alpha = ifelse(annotations$show == TRUE, 1, 0), size = 5) +
  labs(x = "Batch", y = "qD")
  }
}
#plots <- plot_grid(plotlist[[1]], plotlist[[3]],  plotlist[[2]], plotlist[[4]], ncol = 2)
plots <- plot_grid(plotlist[[1]], plotlist[[4]], plotlist[[2]], plotlist[[5]], plotlist[[3]], plotlist[[6]], ncol = 2)
plots
ggsave("./03-Routputfiles/I_Alphadiv/I1-HilldiversityBatchBac.svg", plots, width = 25,
       height = 25, units = "cm", dpi = 300)
################################################################################


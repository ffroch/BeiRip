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
# BACTERIAL PART    
##########################################################################
# phylo object
phyabs <- readRDS("./03-Routputfiles/D_Phyloobjects/D1E_bac_abscnc_nocontrols_nospike.rds")
metadf <- data.frame(sample_data(phyabs))
##########################################################################
# Prepare data for plotting
##########################################################################
# make absbasedonimpt and absbasedonall to one column with pivot_longer
metadflong <- data.frame(pivot_longer(metadf,
                                      cols = c("absbasedonimpt",
                                               "absbasedonall",
                                               "gencopiesbasedonimt",
                                               "gencopiesbasedonall",
                                               "samplemean_counts"),
                                      names_to = "indicator",
                                      values_to = "cellnumber"))
# keep only samplemean_counts, gencopiesbasedonimt and gencopiesbasedonall
metadflongred <- metadflong %>%
  filter(indicator != "absbasedonimpt", indicator != "absbasedonall")
# replace inf with NA
metadflongred[metadflongred == Inf] <- NA
# replace samplemean_counts in indicator with "CFU on PCA"
metadflongred$indicator[metadflongred$indicator == "samplemean_counts"] <- "CFU on PCA"
##########################################################################
# FUNGAL PART
##########################################################################
phyabsf <- readRDS("./03-Routputfiles/D_Phyloobjects/D2C_fun_abs_nocontrols_unfiltered.rds")
metadff <- data.frame(sample_data(phyabsf))
##########################################################################
# Prepare data for plotting
##########################################################################
# make absbasedonimpt and absbasedonall to one column with pivot_longer
metadflongf <- data.frame(pivot_longer(metadff,
                                      cols = c("totalcopies0mm",
                                                "totalcopies1mm",
                                               "samplemean_counts"),
                                      names_to = "indicator",
                                      values_to = "cellnumber"))

# replace inf with NA
metadflongf[metadflongf == Inf] <- NA
# replace samplemean_counts in indicator with "CFU on YGC"
metadflongf$indicator[metadflongf$indicator == "samplemean_counts"] <- "CFU on YGC"
##########################################################################
# combine data
mdfred <- metadflongred %>%
  select("days", "batch", "cellnumber", "indicator", "piece", "slice")
mdffred <- metadflongf %>%
  select("days","batch", "cellnumber", "indicator", "piece", "slice")
alldata <- rbind(mdfred, mdffred)
##########################################################################
# add index numbers to reorder the indicator
alldata$index <- ifelse(alldata$indicator == "gencopiesbasedonall", 1,
                        ifelse(alldata$indicator == "gencopiesbasedonimt", 2,
                               ifelse(alldata$indicator == "CFU on PCA", 3,
                                      ifelse(alldata$indicator == "totalcopies0mm", 4,
                                            ifelse(alldata$indicator == "totalcopies1mm", 5,
                                             ifelse(alldata$indicator == "CFU on YGC", 6, NA))))))
alldata <- alldata[!is.na(alldata$cellnumber), ]
##########################################################################
# plot
##########################################################################
colors <- c("#670A1E", "#9E001B", "#FD5A52", "#f58287", "#C2185B",
            "#7B1FA2", "#512DA8", "#002D41", "#03546C", "#2D829D",
            "#3DB2C8", "#71949B", "#00796B", "#388E3C", "#689F38",
            "#d3b928", "#FFA000", "#F57C00", "#E64A19", "#5D4037",
            "#C9B8A7", "#616161", "#b4b4b4")
batchcol <- c("#9E001B", "#002D41","#388E3C")

# remove CFU on YGC from alldata

p <- ggplot(alldata, aes(days, log10(cellnumber), color = reorder(indicator, index))) +
  geom_point(alpha = 0.7) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.5,
              linetype = "longdash") +
  facet_wrap(~piece) +
  scale_color_manual(values = c("#670A1E", "#9E001B", "#FD5A52", "#002D41", "#2D829D", "#71949B"),
                     labels = c(bquote("genomic equivalents (" * italic(Allobacillus) * ")"), bquote("genomic equivalents (" * italic(Imtechella) * ")"),
                                "CFU on PCA", "18S gene copies (qPCR 0 mismatches)", "18S gene copies (qPCR 1 mismatch)", "CFU on YGC"),
                      name = "measurement type") +
  theme_bw() +
  ylab("log10(measurement type)") +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 8))
print(p)

# save plot
ggsave("./03-Routputfiles/E_CountsVsGenecopies/E_CountsVsGenecopies.svg",
       width = 21, height = 29, units = "cm")
ggsave("./03-Routputfiles/E_CountsVsGenecopies/E_CountsVsGenecopies.png",
       width = 21, height = 29, units = "cm")


##########################################################################
# show the correlations of the methods
##########################################################################
alldata <- alldata[, 1:6]
alldatawide <- data.frame(pivot_wider(alldata, names_from = indicator, values_from = cellnumber))

corr1 <- cor.test(log10(alldatawide$gencopiesbasedonall), log10(alldatawide$CFU.on.PCA))
corr1
x_pos <- min(log10(alldatawide$gencopiesbasedonall), na.rm = TRUE) + 0.995 * (max(log10(alldatawide$gencopiesbasedonall), na.rm = TRUE) - min(log10(alldatawide$gencopiesbasedonall), na.rm = TRUE))
y_pos <- min(log10(alldatawide$CFU.on.PCA), na.rm = TRUE) + 0.1 * (max(log10(alldatawide$CFU.on.PCA), na.rm = TRUE) - min(log10(alldatawide$CFU.on.PCA), na.rm = TRUE))

p1 <- ggplot(alldatawide, aes(x = log10(gencopiesbasedonall), y = log10(CFU.on.PCA))) +
  geom_point(aes(color = factor(batch)),size = 4, alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = batchcol, name = "batch") +
  theme_bw() +
  theme(legend.position.inside = c(0,1),
  legend.justification = c(0,1)) +
  xlab(bquote("log10(genomic equivalents " * italic(Allobacillus) * ")")) +
  ylab("log10(CFU on PCA)")
  
if (corr1$p.value > 0.001) {
  p1 <- p1 + geom_text(aes(label = paste("r = ", round(corr1$estimate, 2),
                                        "\n95% CI [", round(corr1$conf.int[1], 2), ", ",
                                        round(corr1$conf.int[2], 2), "]",
                                        "\np = ", round(corr1$p.value, 3))),
                       x = x_pos, y = y_pos, size = 3, hjust = 1, color = "black")
} else {
  p1 <- p1 + geom_text(aes(label = paste("r = ", round(corr1$estimate, 2),
                                        "\n95% CI [", round(corr1$conf.int[1], 2), ", ",
                                        round(corr1$conf.int[2], 2), "]",
                                        "\np = < 0.001")),
                       x = x_pos, y = y_pos, size = 3, hjust = 1, color = "black")
}

corr2 <- cor.test(log10(alldatawide$gencopiesbasedonimt), log10(alldatawide$CFU.on.PCA))
corr2
x_pos <- min(log10(alldatawide$gencopiesbasedonimt), na.rm = TRUE) + 0.995 * (max(log10(alldatawide$gencopiesbasedonimt), na.rm = TRUE) - min(log10(alldatawide$gencopiesbasedonimt), na.rm = TRUE))
y_pos <- min(log10(alldatawide$CFU.on.PCA), na.rm = TRUE) + 0.1 * (max(log10(alldatawide$CFU.on.PCA), na.rm = TRUE) - min(log10(alldatawide$CFU.on.PCA), na.rm = TRUE))
p2 <- ggplot(alldatawide, aes(x = log10(gencopiesbasedonimt), y = log10(CFU.on.PCA))) +
  geom_point(aes(color = factor(batch)),size = 4, alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = batchcol, name = "batch") +
  theme_bw() +
  theme(legend.position = c(0.1,0.8)) +
  xlab(bquote("log10(genomic equivalents " * italic(Imtechella) * ")")) +
  ylab("log10(CFU on PCA)")
if (corr2$p.value > 0.001) {
  p2 <- p2 + geom_text(aes(label = paste("r = ", round(corr2$estimate, 2),
                                        "\n95% CI [", round(corr2$conf.int[1], 2), ", ",
                                        round(corr2$conf.int[2], 2), "]",
                                        "\np = ", round(corr2$p.value, 3))),
                       x = x_pos, y = y_pos, size = 4, hjust = 1 , color = "black")
} else {
  p2 <- p2 + geom_text(aes(label = paste("r = ", round(corr2$estimate, 2),
                                        "\n95% CI [", round(corr2$conf.int[1], 2), ", ",
                                        round(corr2$conf.int[2], 2), "]",
                                        "\np = < 0.001")),
                       x = x_pos, y = y_pos, size = 4, hjust = 1, color = "black")
}
p2
corr3 <- cor.test(log10(alldatawide$totalcopies0mm), log10(alldatawide$CFU.on.YGC))
corr3
x_pos <- min(log10(alldatawide$totalcopies0mm), na.rm = TRUE) + 0.995 * (max(log10(alldatawide$totalcopies0mm), na.rm = TRUE) - min(log10(alldatawide$totalcopies0mm), na.rm = TRUE))
y_pos <- min(log10(alldatawide$CFU.on.YGC), na.rm = TRUE) + 0.1 * (max(log10(alldatawide$CFU.on.YGC), na.rm = TRUE) - min(log10(alldatawide$CFU.on.YGC), na.rm = TRUE))
p3 <- ggplot(alldatawide, aes(x = log10(totalcopies0mm), y = log10(CFU.on.YGC))) +
  geom_point(aes(color = factor(batch)),size = 4, alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = batchcol, name = "batch") +
  theme_bw() +
  xlab("log10(18S gene copies (qPCR))") +
  ylab("log10(CFU on YGC)")
if (corr3$p.value > 0.001) {
  p3 <- p3 + geom_text(aes(label = paste("r = ", round(corr3$estimate, 2),
                                        "\n95% CI [", round(corr3$conf.int[1], 2), ", ",
                                        round(corr3$conf.int[2], 2), "]",
                                        "\np = ", round(corr3$p.value, 3))),
                       x = x_pos, y = y_pos, size = 3, hjust = 1)
} else {
  p3 <- p3 + geom_text(aes(label = paste("r = ", round(corr3$estimate, 2),
                                        "\n95% CI [", round(corr3$conf.int[1], 2), ", ",
                                        round(corr3$conf.int[2], 2), "]",
                                        "\np = < 0.001")),
                       x = x_pos, y = y_pos, size = 3, hjust = 1)
}
#p3 <- p3 + geom_text(aes(x = log10(totalcopies0mm), y = log10(CFU.on.YGC), label = paste0(piece,"_",days)), size = 2, hjust = 1)

# remove alldatawide with B7, B8 and B9 in piece
alldatawide <- alldatawide[alldatawide$piece != "B7" & alldatawide$piece != "B8" & alldatawide$piece != "B9", ]

corr4 <- cor.test(log10(alldatawide$totalcopies0mm), log10(alldatawide$CFU.on.YGC))
corr4
x_pos <- min(log10(alldatawide$totalcopies0mm), na.rm = TRUE) + 0.995 * (max(log10(alldatawide$totalcopies0mm), na.rm = TRUE) - min(log10(alldatawide$totalcopies0mm), na.rm = TRUE))
y_pos <- min(log10(alldatawide$CFU.on.YGC), na.rm = TRUE) + 0.1 * (max(log10(alldatawide$CFU.on.YGC), na.rm = TRUE) - min(log10(alldatawide$CFU.on.YGC), na.rm = TRUE))
p4 <- ggplot(alldatawide, aes(x = log10(totalcopies0mm), y = log10(CFU.on.YGC))) +
  geom_point(aes(color = factor(batch)),size = 4, alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = batchcol, name = "batch") +
  theme_bw() +
  xlab("log10(18S gene copies (qPCR))") +
  ylab("log10(CFU on YGC)")
if (corr4$p.value > 0.001) {
  p4 <- p4 + geom_text(aes(label = paste("r = ", round(corr4$estimate, 2),
                                        "\n95% CI [", round(corr4$conf.int[1], 2), ", ",
                                        round(corr4$conf.int[2], 2), "]",
                                        "\np = ", round(corr4$p.value, 3))),
                       x = x_pos, y = y_pos, size = 3, hjust = 1)
} else {
  p4 <- p4 + geom_text(aes(label = paste("r = ", round(corr4$estimate, 2),
                                        "\n95% CI [", round(corr4$conf.int[1], 2), ", ",
                                        round(corr4$conf.int[2], 2), "]",
                                        "\np = < 0.001")),
                       x = x_pos, y = y_pos, size = 3, hjust = 1)
}
#p4 <- p4 + geom_text(aes(x = log10(totalcopies0mm), y = log10(CFU.on.YGC), label = paste0(piece,"_",days)), size = 2, hjust = 1)




# identify samples with more than 10^6 CFU of Serratia in the 16S data
phy <- readRDS("./03-Routputfiles/D_Phyloobjects/D1F_bac_abscnc_nc_ns_genus.rds")
otu <- otu_table(phy)
tax <- tax_table(phy)
genus_index <- which(tax[,  "Genus"] == "Serratia")

serratia_counts <- colSums(otu[ genus_index, , drop = FALSE])
high_serratia_samples <- names(serratia_counts[serratia_counts > 10^6])
high_serratia_samples

genus_index <- which(tax[,  "Genus"] == "Pseudomonas")
pseudomonas_counts <- colSums(otu[ genus_index, , drop = FALSE])
high_pseudomonas_samples <- names(pseudomonas_counts[pseudomonas_counts > 10^6])
high_pseudomonas_samples

alldata$selector <- paste0("FR_", alldata$piece, "_0", alldata$slice)
alldata$selector <- gsub("_010", "_10", alldata$selector)

alldatasel <- alldata[!(alldata$selector %in% high_serratia_samples), ]
alldatasel <- alldatasel[!(alldatasel$selector %in% high_pseudomonas_samples), ]
alldatasel <- alldatasel[, 1:6]
alldataselwide <- data.frame(pivot_wider(alldatasel, names_from = indicator, values_from = cellnumber))

corr5 <- cor.test(log10(alldataselwide$totalcopies0mm), log10(alldataselwide$CFU.on.YGC))
corr5
x_pos <- min(log10(alldataselwide$totalcopies0mm), na.rm = TRUE) + 0.995 * (max(log10(alldataselwide$totalcopies0mm), na.rm = TRUE) - min(log10(alldataselwide$totalcopies0mm), na.rm = TRUE))
y_pos <- min(log10(alldataselwide$CFU.on.YGC), na.rm = TRUE) + 0.1 * (max(log10(alldataselwide$CFU.on.YGC), na.rm = TRUE) - min(log10(alldataselwide$CFU.on.YGC), na.rm = TRUE))
p5 <- ggplot(alldataselwide, aes(x = log10(totalcopies0mm), y = log10(CFU.on.YGC))) +
  geom_point(aes(color = factor(batch)),size = 4, alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
  scale_color_manual(values = batchcol, name = "batch") +
  theme_bw() +
  xlab("log10(18S gene copies (qPCR))") +
  ylab("log10(CFU on YGC)")
if (corr5$p.value > 0.001) {
  p5 <- p5 + geom_text(aes(label = paste("r = ", round(corr5$estimate, 2),
                                        "\n95% CI [", round(corr5$conf.int[1], 2), ", ",
                                        round(corr5$conf.int[2], 2), "]",
                                        "\np = ", round(corr5$p.value, 3))),
                       x = x_pos, y = y_pos, size = 3, hjust = 1)
} else {
  p5 <- p5 + geom_text(aes(label = paste("r = ", round(corr5$estimate, 2),
                                        "\n95% CI [", round(corr5$conf.int[1], 2), ", ",
                                        round(corr5$conf.int[2], 2), "]",
                                        "\np = < 0.001")),
                       x = x_pos, y = y_pos, size = 3, hjust = 1)
}
#p5 <- p5 + geom_text(aes(x = log10(totalcopies0mm), y = log10(CFU.on.YGC), label = paste0(piece,"_",days)), size = 2, hjust = 1)

p0 <- plot_grid(p1, p2, p3, p4, p5, ncol = 2, labels = letters[1:5])
p0

ggsave("./03-Routputfiles/E_CountsVsGenecopies/E_CountsVsGenecopies_corr.svg",
       width = 21, height = 29, units = "cm")
ggsave("./03-Routputfiles/E_CountsVsGenecopies/E_CountsVsGenecopies_corr.png",
       width = 21, height = 29, units = "cm")

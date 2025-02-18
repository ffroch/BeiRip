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
phyraw <- readRDS("./03-Routputfiles/D_Phyloobjects/D2A_fun_raw_nocontrols.rds")
# get the otu tables
oturaw <- data.frame(otu_table(phyraw))
oturaw_t <- t(oturaw)
# calculate coverage
cov_per_sample <- apply(oturaw_t, 1, function(x) Chat.Ind(x, sum(x)))
cov_df <- data.frame(sample = rownames(oturaw_t), coverage = cov_per_sample)
cov_df
# estimateD
set.seed(123)
bycoverageraw <- estimateD(oturaw, base = "coverage", level = 0.995, nboot = 50)
hilldiv <- bycoverageraw
# get meta data
meta <- data.frame(sample_data(phyraw))
# add meta
hilldiv_m <- dplyr::left_join(hilldiv, meta, by = c("Assemblage" = "sampleid"))
##########################################################################
hilldiv_m$Order.q <- factor(hilldiv_m$Order.q,
                            levels = c(0, 1, 2),
                            labels = c("species richness", "Hill-Shannon",
                                       "Hill-Simpson"))

# plotting
# plot alpha diversity over time
#hilldiv_m <- hilldiv_m[hilldiv_m$Order.q != "species richness", ]
#hilldiv_m <- droplevels(hilldiv_m)
p1 <- ggplot(hilldiv_m, aes(x = factor(days), y = qD)) +
  geom_boxplot(outlier.shape = NA, aes(color = factor(batch), fill = factor(batch))) +
  geom_boxplot(outlier.shape = NA,  color = alpha("black", 1), alpha = 0.5) +
  #geom_jitter(width = 0.3, alpha = 0.5, aes(color = factor(batch))) +
  facet_grid(Order.q ~ ., scales = "free_y") +
  scale_color_manual(values = c("#9E001B", "#002D41", "#388E3C"), name = "batch") +
  scale_fill_manual(values = c("#9E001B", "#002D41", "#388E3C"), name = "batch") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white")) +
  labs(x = "days", y = "qD")
print(p1)
# save as tiff 25x25 cm
ggsave("./03-Routputfiles/I_Alphadiv/I2-HilldiversitySliceFun_boxbox.svg", p1, width = 25,
       height = 25, units = "cm", dpi = 300)
ggsave("./03-Routputfiles/I_Alphadiv/I2-HilldiversitySliceFun_boxbox.png", p1, width = 25,
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
  facet_grid(Order.q ~ ., scales = "free_y") +
  scale_color_manual(values = c("#9E001B", "#002D41", "#388E3C"), name = "batch") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"))
print(p2)
ggsave("./03-Routputfiles/I_Alphadiv/I2-HilldiversitySliceFun_boxjit.svg", p2, width = 25,
       height = 25, units = "cm", dpi = 300)
ggsave("./03-Routputfiles/I_Alphadiv/I2-HilldiversitySliceFun_boxjit.png", p2, width = 25,
       height = 25, units = "cm", dpi = 300)
       

##########################################################################
# shannon only
hilldiv_f <- hilldiv_m[hilldiv_m$Order.q == "Hill-Shannon", ]

# plotting
# alternative plotting
hilldiv_f <- hilldiv_f %>%
  mutate(
    position = case_when(
      batch == 1 ~ as.numeric(as.factor(days)) - 0.25,
      batch == 2 ~ as.numeric(as.factor(days)),
      batch == 3 ~ as.numeric(as.factor(days)) + 0.25
    )
  )

p2 <- ggplot(hilldiv_f, aes(x = factor(days), y = qD)) +
  geom_boxplot(outlier.shape = NA,  color = alpha("black", 1), alpha = 0.1) +
  geom_point(alpha = 0.7, aes(color = factor(batch), x = position)) +
  scale_color_manual(values = c("#9E001B", "#002D41", "#388E3C"), name = "batch") +
  theme_bw() +
  xlab("days") +
  ylab("Hill-Shannon diversity")
print(p2)
ggsave("./03-Routputfiles/I_Alphadiv/I2-HilldiversitySliceFun_boxjit_shannononly.svg", p2, width = 25,
       height = 25, units = "cm", dpi = 300)
ggsave("./03-Routputfiles/I_Alphadiv/I2-HilldiversitySliceFun_boxjit_shannononly.png", p2, width = 25,
       height = 25, units = "cm", dpi = 300)
saveRDS(p2, "./03-Routputfiles/I_Alphadiv/I2-HilldiversitySliceFun_boxjit_shannononly.rds")
################################################################################
# comparison alpha diversity over batch
# create new group for comparison
kwtest <- list()
dunntest <- list()
plotlist <- list()
for (i in unique(hilldiv_m$Order.q)) {
  kwtest[[i]] <- kruskal.test(qD ~ batch, data = hilldiv_m[hilldiv_m$Order.q == i, ])
  dunntest[[i]] <- dunn.test::dunn.test(hilldiv_m$qD[hilldiv_m$Order.q == i],
                                       hilldiv_m$batch[hilldiv_m$Order.q == i],
                                       method = "bonferroni")
  annotations <- data.frame(comparison = dunntest[[i]]$comparisons,
                            p.value = dunntest[[i]]$P.adjusted,
                            label = ifelse(dunntest[[i]]$P.adjusted < 0.025, "*", ""),
                            pkw = rep(kwtest[[i]]$p.value, length(dunntest[[i]]$P.adjusted))) %>%
                            separate(comparison, into = c("BatchA", "BatchB"), sep = " - ")
  annotations$show <- ifelse(annotations$p.value < 0.025 & annotations$pkw < 0.05, TRUE, FALSE)
  ymax <- max(hilldiv_m[hilldiv_m$Order.q == i, ]$qD)
  annotations$y_position <- c(ymax * 1.05, ymax * 1.1, ymax * 1.15)
  annotations$formatted_pvalue <- sprintf("%.3f", annotations$p.value)
  plotlist[[i]] <- ggplot(hilldiv_m[hilldiv_m$Order.q == i, ], aes(x = factor(batch), y = qD)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.3, alpha = 0.5, aes(color = days)) +
  facet_grid(Order.q ~ ., scales = "free_y") +
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

#plots <- plot_grid(plotlist[[1]],  plotlist[[2]],  ncol = 1)
plots <- plot_grid(plotlist[[1]],  plotlist[[2]],  plotlist[[3]],  ncol = 1)
plots
ggsave("./03-Routputfiles/I_Alphadiv/I2-HilldiversityBatchFun.svg", plots, width = 25,
       height = 25, units = "cm", dpi = 300)
################################################################################
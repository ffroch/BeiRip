library(tidyverse)
data <- read.table("./02-rawdata/Co_Culture_qPCR_finalResults.txt", header = T, sep = "\t")

data$ratiobac <- log2(data$bac_fun_16S/data$bac_16S)
data$ratiofun <- log2(data$bac_fun_18S/data$fun_18S)
data$ratiobaccfu <- log2(data$bac_cfu_co/data$bac_cfu)
data$ratiofuncfu <- log2(data$fun_cfu_co/data$fun_cfu)

ggplot(data, aes(x = log10(bac_fun_16S), y = log10(bac_cfu_co))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
ggplot(data, aes(x = log10(bac_fun_18S), y = log10(fun_cfu_co))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
ggplot(data, aes(x = log10(bac_16S), y = log10(bac_cfu))) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1)
ggplot(data, aes(x = log10(fun_18S), y = log10(fun_cfu))) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1)

ggplot(data, aes(x = ratiobac, y = ratiobaccfu)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1)
ggplot(data, aes(x = ratiofun, y = ratiofuncfu)) +
    geom_point() +
    geom_abline(intercept = 0, slope = 1)

datalong <- pivot_longer(data, cols = c(ratiobac, ratiofun, ratiobaccfu, ratiofuncfu), names_to = "ratio", values_to = "value")

meta <- read.table("./01-metadata/isolateinfo.txt", header = T, sep = "\t")
metashort <- meta[, c("shortIDreplacement", "GenusWGS", "Group")]
colnames(metashort) <- c("bac", "bacGenus", "Group")

datalong <- data.frame(datalong)
datalong <- left_join(datalong, metashort, by = c("bac" = "bac"))
metaus <- meta[, c("shortIDreplacement", "GenusWGS")]
colnames(metaus) <- c("fun", "funGenus")
datalong <- left_join(datalong, metaus, by = c("fun" = "fun"))

ggplot(datalong, aes(x = fun, y = value, color = bac)) +
  geom_boxplot() +
  facet_wrap(~ratio)

ggplot(datalong, aes(x = fun, y = value, color = Group)) +
  geom_boxplot() +
  facet_wrap(~ratio)

ggplot(datalong, aes(x = ratio, y = value, color = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(aes(group = Group), position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75))

datalong$group2 <- ifelse(grepl("bac", datalong$ratio), "bacteria", "fungi")
datalong$highlights <- ifelse(datalong$bac == "ID-1609", "ID-1609", "other strains")
datalong$highlights <- ifelse(datalong$bac == "ID-1626", "ID-1626", datalong$highlights)
datalong$highlights <- ifelse(datalong$bac == "ID-1456", "ID-1456", datalong$highlights)
datalong$labs <- ifelse(datalong$ratio == "ratiobac", "qPCR 16S", 
                         ifelse(datalong$ratio == "ratiofun", "qPCR 18S", 
                                ifelse(datalong$ratio == "ratiobaccfu", "CFU", "CFU")))

ggplot(datalong, aes(x = labs, y = value, fill = Group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.3) +
  geom_jitter(aes(group = Group, color = highlights, shape = fun), alpha = 0.8, position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.75)) +
    facet_wrap(~group2, scales = "free_x") +
    theme_bw() +
    scale_fill_manual(values = c("white", "grey"), name = "group") +
    scale_color_manual(values = c("#D45500", "#004455", "#00738E", "#626567"), name = "bacterial partner") +
    scale_shape_manual(values = c(15, 16, 17, 3), name = "fungal partner") +
    ylab("log2 ratio") +
    xlab("") +
    theme(strip.background = element_rect(fill = "transparent", color = NA),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          legend.key = element_rect(fill = "transparent", color = NA),
          legend.background = element_rect(fill = "transparent", color = NA))

ggsave("./03-Routputfiles/S_CoCulture/S_qPCR.svg", width = 16, height = 16, units = "cm")
ggsave("./03-Routputfiles/S_CoCulture/S_qPCR.png", width = 16, height = 16, units = "cm", dpi = 300)

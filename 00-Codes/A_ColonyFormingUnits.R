# A_ColonyFormingUnits.R
# Plot the colony forming units (CFU) of the bacteria and fungi
# in the different media over time.
###############################################################################
library(ggplot2)
library(dplyr)#for "%>%"
library(cowplot)#for plot_grid
###############################################################################
# import data
data <- read.table("./02-rawdata/Plate_counts_txt.txt",
                   header = TRUE, sep = "\t")
key <- read.table("./01-metadata/slice_day_key.txt", header = TRUE, sep = "\t")
key$days <- key$days
###############################################################################
# filter the data - RCM is not necessary to focus on
# remove manually excluded counts (wrong plating dilution)
data <- data %>%
  filter(ex3 == 0, medium != "RCM")
data <- plyr::join(data, key, type = "left")
###############################################################################
# plot all counts of different dilutions including under or overestimation
ggplot(data, aes(factor(days), undest_logcounts,
                 color = factor(underestimation))) +
  geom_point() +
  geom_boxplot() +
  facet_grid(. ~ medium)
###############################################################################
# split data to fungal and bacterial counts
#select only data from medium BR, YGC25 and YGC50
fungi <- data %>%
  filter(medium == "BR" | medium == "YGC25" | medium == "YGC4")
#calculate the mean log_counts for each sampleID, medium and day
fmc <- fungi %>%
  group_by(sampleID, medium, days) %>%
  summarise(samplemean_logcounts = mean(log_counts, na.rm = TRUE))
fmc <- data.frame(fmc)
# save fmc as rds
saveRDS(fmc, "./03-Routputfiles/A_PlateCounts/A_fungalplatecounts.rds")
###############################################################################
#select only data from medium MRS and PCA
bacteria <- data %>%
  filter(medium == "MRS" | medium == "PCA")
#calculate the mean log_counts for each sampleID, medium and day
bmc <- bacteria %>%
  group_by(sampleID, medium, days) %>%
  summarise(samplemean_logcounts = mean(log_counts, na.rm = TRUE))
bmc <- data.frame(bmc)
# save bmc as rds
saveRDS(bmc, "./03-Routputfiles/A_PlateCounts/A_bacterialplatecounts.rds")
###############################################################################
# plot everything together
#add a column to data2 called kingdom with "fungi" for medium BR, YGC25 and
# YGC50 and "bacteria" for medium MRS and PCA
data$kingdom <- ifelse(data$medium == "BR" | data$medium == "YGC25" |
                        data$medium == "YGC4", "fungi", "bacteria")
p1 <- data %>%
filter(kingdom == "fungi") %>%
ggplot(aes(days, log_counts)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(se = TRUE, linewidth = 0.5, color = "darkred") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
        strip.text = element_text(size = 12)) +
  ylab("log10 cfu") +
  facet_grid(. ~ medium)
p1
#remove data from day 9 of medium VRBD, there went something wrong
data2 <- data %>%
  filter(!(medium == "VRBD" & slice == 9))
p2 <- data2 %>%
  filter(kingdom == "bacteria") %>%
ggplot(aes(days, log_counts)) +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(se = TRUE, linewidth = 0.5, color = "darkred") +
  theme_bw() +
  theme(strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12)) +
  ylab("log10 cfu") +
  facet_grid(. ~ medium)
p2

#combine plots
plot_grid(p2, p1, align = "v", nrow = 2)
#save plot
ggsave("./03-Routputfiles/A_PlateCounts/A_growthcurves.svg",
       width = 20, height = 20, units = "cm", dpi = 300)
ggsave("./03-Routputfiles/A_PlateCounts/A_growthcurves.png",
       width = 20, height = 20, units = "cm", dpi = 300)
###############################################################################
#make a summary table with the mean and standard deviation of
# log_counts for each medium and day
summary <- data %>%
  group_by(medium, days) %>%
  summarise(mean_logcounts = mean(log_counts, na.rm = TRUE),
  sd_logcounts = sd(log_counts, na.rm = TRUE))
#save summary table
write.table(summary,
            "./03-Routputfiles/A_PlateCounts/A_growthcurves_summary.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
##############################################################################
#make a new summary table only with days 1, 16 and 56
summary2 <- summary %>%
  filter(days == 0 | days == 15 | days == 55)
summary2 <- data.frame(summary2)
#paste the medium and days columns together
summary2$md <- paste0(summary2$medium, summary2$days)
#round the mean_logcounts and sd_logcounts columns to 2 decimal places
summary2$mean_logcounts <- round(summary2$mean_logcounts, 2)
summary2$sd_logcounts <- round(summary2$sd_logcounts, 2)
#save summary2 as table
saveRDS(summary2, "./03-Routputfiles/A_PlateCounts/A_growthcurves_summary2.rds")
###############################################################################
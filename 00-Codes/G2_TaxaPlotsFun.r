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
phyabscorr <- readRDS("./03-Routputfiles/D_Phyloobjects/D2C_fun_abs_nocontrols.rds")
##########################################################################
# Prepare data
##########################################################################
# transpose otufabs
otufabst <- data.frame(t(otu_table(phyabscorr)))
otufabst$sampleid <- rownames(otufabst)
# add metadata
metafabs <- data.frame(sample_data(phyabscorr))
names <- rownames(metafabs)
# combine otu table and meta data
df <- plyr::join(metafabs, otufabst, by = c("sampleid"))
# get taxonomy
taxa <- data.frame(tax_table(phyabscorr))
taxa$sseqid <- rownames(taxa)
# sort taxa by Kingdom, Phylum, Class, Order, Family, Genus = necessary for plotting
taxa <- taxa %>% arrange(Kingdom, Phylum, Class, Order, Family, Genus)
taxa$newid <- seq.int(nrow(taxa))
taxa$newid <- paste0("id", taxa$newid)
##########################################################################
# calculate the mean per day (col 21 = days, col 25 to end = ASVs)
##########################################################################
# remove columns that are not needed
di <- which(colnames(df) == "days")
start <- ncol(metafabs) + 1
dfred <- df[, c(di, start:ncol(df))]
# calculate the mean per days
dfmean <- dfred %>%
  group_by(days) %>%
  summarise(across(everything(), \(x) mean(x, na.rm = TRUE)))
# make it a data frame
dfmean <- data.frame(dfmean)
##########################################################################
# Transform data for plotting
##########################################################################
# transform dfmean to long format
dfmeanlong <- data.frame(pivot_longer(dfmean, cols = c(2:ncol(dfmean)),
                                      names_to = "sseqid", values_to = "meancounts"))
# remove X from sseqid if it starts with X
dfmeanlong$sseqid <- ifelse(grepl("X", dfmeanlong$sseqid),
                           sub("X", "", dfmeanlong$sseqid),
                           dfmeanlong$sseqid)
# add taxonomy
dfmeanlong <- plyr::join(dfmeanlong, taxa, by = "sseqid")
# calculate the total sum of all ASVs per day
cumsum <- dfmeanlong %>%
  group_by(days) %>%
  summarise(total = log10(sum(meancounts, na.rm = TRUE)))
cumsum <- data.frame(cumsum)
# calculate the proportion of each ASV per day
dfmeanlong <- dfmeanlong %>%
  group_by(days) %>%
  mutate(total_value = sum(meancounts, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(proportion = meancounts / total_value)
# merge cumsum with dfmeanlong
dfmeanlong <- merge(dfmeanlong, cumsum, by = "days")
##########################################################################
# calculate the maximimum proportion of each sseqid 
maxproportion <- dfmeanlong %>%
  group_by(Genus) %>%
  summarise(maxproportion = max(proportion, na.rm = TRUE))
# merge maxproportion with dfmeanlong
dfmeanlong <- plyr::join(dfmeanlong, maxproportion, by = "Genus")
# create a new column and name low abundant genera as other
dfmeanlong$newgenus <- ifelse(dfmeanlong$maxproportion > 0.01,
                              dfmeanlong$Genus, "Others")
##########################################################################
# import data if isolates
frequency <- read.table("./03-Routputfiles/F_blastn_fun_seqvsiso_frequencies.txt", sep = "\t", header = TRUE)
colnames(frequency)[1] <- "sseqid"
dfmeanlong <- plyr::join(dfmeanlong, frequency, by = "sseqid")
# create a column in dfmeanlong that contains a 1 if the frequency is not na. else a 0.5
dfmeanlong$iso <- ifelse(is.na(dfmeanlong$Freq), "noiso", "iso")
# save dfmeanlong for later use
saveRDS(dfmeanlong, "./03-Routputfiles/G_Taxaplots/G2-dfmeanlongFun.rds")
##########################################################################
owncolgen <- readRDS("./03-Routputfiles/Z_Helperfiles/owncolgenfun.rds")
dfmeanlong$generataxaplot <- ifelse(dfmeanlong$Genus %in% owncolgen, dfmeanlong$Genus, "Others")
# further corrections
dfmeanlong$generataxaplot <- ifelse(dfmeanlong$generataxaplot == "Debaryomycetaceae", "uncl. Debaryomycetaceae", dfmeanlong$generataxaplot)
dfmeanlong$generataxaplot <- ifelse(dfmeanlong$generataxaplot == "Trichosporonaceae", "uncl. Trichosporonaceae", dfmeanlong$generataxaplot)
dfmeanlong$generataxaplot <- ifelse(dfmeanlong$generataxaplot == "Dothideomycetes", "uncl. Dothideomycetes", dfmeanlong$generataxaplot)
dfmeanlong$generataxaplot <- ifelse(dfmeanlong$generataxaplot == "Cystobasidiales", "uncl. Cystobasidiales", dfmeanlong$generataxaplot)
dfmeanlong$generataxaplot <- ifelse(dfmeanlong$generataxaplot == "Leotiomycetes", "uncl. Dothideomycetes", dfmeanlong$generataxaplot)
dfmeanlong$generataxaplot <- ifelse(dfmeanlong$generataxaplot == "Plectosphaerellaceae", "uncl. Plectosphaerellaceae", dfmeanlong$generataxaplot)
dfmeanlong$generataxaplot <- ifelse(dfmeanlong$generataxaplot == "Pleosporales", "uncl. Pleosporales", dfmeanlong$generataxaplot)
##########################################################################
# sort dfmeanlong by taxonomy
dfmeanlong <- dfmeanlong %>% arrange(Kingdom, Phylum, Class, Order, Family, generataxaplot)
dfmeanlong$index <- seq.int(nrow(dfmeanlong))
dfmeanlong$index <- ifelse(dfmeanlong$generataxaplot == "Others", 100000, dfmeanlong$index)
##########################################################################
# create a colorset
colors <- c("#670A1E", "#9E001B", "#FD5A52", "#f58287", "#C2185B",
            "#7B1FA2", "#512DA8", "#002D41", "#03546C", "#2D829D",
            "#3DB2C8", "#71949B", "#00796B", "#388E3C", "#689F38",
            "#d3b928", "#FFA000", "#F57C00", "#E64A19", "#5D4037",
            "#C9B8A7", "#616161", "#b4b4b4")
##########################################################################
# plot taxaplots
##########################################################################
# adapt the labels for the plot
# Function to create custom labels
legendlabel<-unique(dfmeanlong$generataxaplot)
legendlabel<-legendlabel[c(2:(length(legendlabel)),1)]
llsplit1 <- ifelse(grepl("uncl.", legendlabel), "uncl.", "")
llsplit2 <- legendlabel
llsplit2 <-gsub("uncl. ", "", llsplit2)
legendlabel<-ape::mixedFontLabel(llsplit1, llsplit2,sep=" ",italic=2,always.upright = c("Others", "Dothideomycetes", "Cystobasidiales", "Leotiomycetes", "Pleosporales"))
##########################################################################
# additional the same as bar plot
p5 <- ggplot(dfmeanlong, aes(factor(days), proportion * total, fill = reorder(generataxaplot, index)),
             col = "lightgrey") +
  geom_bar(stat = "identity") +
  theme_bw() +
  scale_fill_manual(values = colors, labels = legendlabel) +
  #scale_alpha_manual(values = c(1, 0.5)) +
  #scale_color_manual(values = c("red","black")) +
  geom_col(color = ifelse(dfmeanlong$iso == "iso", "black", "lightgrey"),
           linewidth = ifelse(dfmeanlong$iso == "iso", 0.5, 0.0001),
           alpha = ifelse(dfmeanlong$iso == "iso", 1, 1)) +
  geom_col(color = ifelse(dfmeanlong$iso == "iso", "black", "lightgrey"),
           linewidth = ifelse(dfmeanlong$iso == "iso", 0.5, 0.0001),
           linetype = ifelse(dfmeanlong$iso == "iso", "solid", "blank"),
           alpha = ifelse(dfmeanlong$iso == "iso", 1, 0)) +
  guides(fill = guide_legend(title = "Genera", ncol = 1)) +
  xlab("days") +
  ylab("log10(estimated 18S rRNA gene copies)")
print(p5)
# save plot as tiff 25x25cm
ggsave("./03-Routputfiles/G_Taxaplots/G2_fun_mean_time.svg", p5,
       width = 20, height = 20, units = "cm")
ggsave("./03-Routputfiles/G_Taxaplots/G2_fun_mean_time.png", p5,
       width = 20, height = 20, units = "cm")
saveRDS(p5, "./03-Routputfiles/G_Taxaplots/G2_fun_mean_time.rds")
###############################################################################
# Individual Plotting
###############################################################################
# transform dfmean to long format
dflong <- data.frame(pivot_longer(df, cols = c(start:ncol(df)),
                                      names_to = "sseqid", values_to = "mean"))
# remove X from sseqid if it starts with X
dflong$sseqid <- ifelse(grepl("X", dflong$sseqid),
                           sub("X", "", dflong$sseqid),
                           dflong$sseqid)
# add taxonomy
dflong <- plyr::join(dflong, taxa[,1:8], by = "sseqid")
##########################################################################
# split up the data on slice
plot_list <- list()
for (i in 1:10) {
  owncolgen <- readRDS("./03-Routputfiles/Z_Helperfiles/owncolgenfun.rds")
dflong$generataxaplot <- ifelse(dflong$Genus %in% owncolgen, dflong$Genus, "Others")
# further corrections
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Debaryomycetaceae", "uncl. Debaryomycetaceae", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Trichosporonaceae", "uncl. Trichosporonaceae", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Dothideomycetes", "uncl. Dothideomycetes", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Cystobasidiales", "uncl. Cystobasidiales", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Leotiomycetes", "uncl. Dothideomycetes", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Plectosphaerellaceae", "uncl. Plectosphaerellaceae", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Pleosporales", "uncl. Pleosporales", dflong$generataxaplot)
# sort dflong by taxonomy
dflong <- dflong %>% arrange(Kingdom, Phylum, Class, Order, Family, generataxaplot)
dflong$index <- seq.int(nrow(dflong))
dflong$index <- ifelse(dflong$generataxaplot == "Others", 10000000000, dflong$index)

  dflongsub <- dflong %>%
    filter(slice == i)
# calculate the total sum of all ASVs per day and piece
cumsum <- dflongsub %>%
  group_by(piece) %>%
  summarise(total = log10(sum(mean, na.rm = TRUE)))
cumsum <- data.frame(cumsum)
# calculate the proportion of each ASV per day and piece
dflongsub <- dflongsub %>%
  group_by(piece) %>%
  mutate(total_value = sum(mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(proportion = mean / total_value)
# merge cumsum with dflong
dflongsub <- merge(dflongsub, cumsum, by = c("piece"))

##########################################################################
# import data if isolates
frequency <- read.table("./03-Routputfiles/F_blastn_fun_seqvsiso_frequencies.txt", sep = "\t", header = TRUE)
colnames(frequency)[1] <- "sseqid"
dflongsub <- plyr::join(dflongsub, frequency, by = "sseqid")
# create a column in dflongsub that contains a 1 if the frequency is not na. else a 0.5
dflongsub$iso <- ifelse(is.na(dflongsub$Freq), "noiso", "iso")
##########################################################################
##########################################################################
# plot taxaplots
##########################################################################
# adapt the labels for the plot
# Function to create custom labels
legendlabel<-unique(dflong$generataxaplot)
legendlabel<-legendlabel[c(2:(length(legendlabel)),1)]
llsplit1 <- ifelse(grepl("uncl.", legendlabel), "uncl.", "")
llsplit2 <- legendlabel
llsplit2 <-gsub("uncl. ", "", llsplit2)
legendlabel<-ape::mixedFontLabel(llsplit1, llsplit2,sep=" ",italic=2,always.upright = c("Others"))

  p <- ggplot(dflongsub, aes(factor(piece), proportion * total, fill = reorder(generataxaplot,index)),
               col = "lightgrey") +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = colors) +
    geom_col(color = ifelse(dflongsub$iso == "iso", "black", "lightgrey"),
             linewidth = ifelse(dflongsub$iso == "iso", 0.5, 0.0001),
             alpha = ifelse(dflongsub$iso == "iso", 1, 1)) +
    geom_col(color = ifelse(dflongsub$iso == "iso", "black", "lightgrey"),
             linewidth = ifelse(dflongsub$iso == "iso", 0.5, 0.0001),
             linetype = ifelse(dflongsub$iso == "iso", "solid", "blank"),
             alpha = ifelse(dflongsub$iso == "iso", 1, 0))+
    guides(fill = guide_legend(title = "Genera", ncol = 1)) +
    xlab("") +
    ylab("") +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) 
  plot_list[[i]] <- p
}
 p <- ggplot(dflongsub, aes(factor(piece), proportion * total, fill = reorder(generataxaplot,index)),
               col = "lightgrey") +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = colors, labels = legendlabel) +
    geom_col(color = ifelse(dflongsub$iso == "iso", "black", "lightgrey"),
             linewidth = ifelse(dflongsub$iso == "iso", 0.5, 0.0001),
             alpha = ifelse(dflongsub$iso == "iso", 1, 1)) +
    geom_col(color = ifelse(dflongsub$iso == "iso", "black", "lightgrey"),
             linewidth = ifelse(dflongsub$iso == "iso", 0.5, 0.0001),
             linetype = ifelse(dflongsub$iso == "iso", "solid", "blank"),
             alpha = ifelse(dflongsub$iso == "iso", 1, 0)) +
    guides(fill = guide_legend(title = "Genera", ncol = 3)) +
    xlab("") +
    ylab("")
legend <- get_legend(p)
p1 <- plot_grid(legend)
plot_list[[i+1]] <- p1

combined_plot1 <- plot_grid(plotlist = plot_list[1:9], ncol = 3, nrow = 3, labels = letters[1:9], label_size = 12, hjust = -1, vjust = 1.5)
# We need to create a nested grid for the last two plots
bottom_left_plot <- plot_grid(plot_list[[10]], align = 'v')
bottom_right_plot <- plot_grid(plot_list[[11]], align = 'v')

# Combine the last two plots into one row with different widths
combined_plot2 <- plot_grid(bottom_left_plot, bottom_right_plot, ncol = 2, rel_widths = c(1, 2), labels = c("J", ""), label_size = 12, hjust = -1, vjust = 1.5)

combined_plot <- plot_grid(combined_plot1, combined_plot2, nrow = 2, rel_heights = c(3, 1))

# Create an empty drawing layer with the combined plot
# Add margins around the grid to ensure space for labels
final_grid_with_margins <- insert_xaxis_grob(combined_plot, grob = grid::nullGrob(), position = "top", height = grid::unit(0.05, "null"))
final_grid_with_margins <- insert_yaxis_grob(final_grid_with_margins, grob = grid::nullGrob(), position = "left", width = grid::unit(0.05, "null"))

# Now add overarching labels using ggdraw()
final_plot <- ggdraw(final_grid_with_margins) +
  draw_label("meat piece", x = 0.5, y = 0.23, hjust = 0.5, vjust = 0) +
  draw_label("log10(estimated 18S rRNA gene copies)", x = 0.03, y = 0.5, angle = 90, hjust = 0.5, vjust = 1)

# save plot as svg
ggsave("./03-Routputfiles/G_Taxaplots/G2_fun_individual_piece.svg", final_plot,
       width = 27.72, height = 38.94, units = "cm")
ggsave("./03-Routputfiles/G_Taxaplots/G2_fun_individual_piece.png", final_plot,
       width = 27.72, height = 38.94, units = "cm")

##########################################################################
plot_list <- list()
for (i in c("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9")) {
    owncolgen <- readRDS("./03-Routputfiles/Z_Helperfiles/owncolgenfun.rds")
dflong$generataxaplot <- ifelse(dflong$Genus %in% owncolgen, dflong$Genus, "Others")
# further corrections
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Debaryomycetaceae", "uncl. Debaryomycetaceae", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Trichosporonaceae", "uncl. Trichosporonaceae", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Dothideomycetes", "uncl. Dothideomycetes", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Cystobasidiales", "uncl. Cystobasidiales", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Leotiomycetes", "uncl. Dothideomycetes", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Plectosphaerellaceae", "uncl. Plectosphaerellaceae", dflong$generataxaplot)
dflong$generataxaplot <- ifelse(dflong$generataxaplot == "Pleosporales", "uncl. Pleosporales", dflong$generataxaplot)
# sort dflong by taxonomy
dflong <- dflong %>% arrange(Kingdom, Phylum, Class, Order, Family, generataxaplot)
dflong$index <- seq.int(nrow(dflong))
dflong$index <- ifelse(dflong$generataxaplot == "Others", 10000000000, dflong$index)

  dflongsub <- dflong %>%
    filter(piece == i)
# calculate the total sum of all ASVs per day and piece
cumsum <- dflongsub %>%
  group_by(slice) %>%
  summarise(total = log10(sum(mean, na.rm = TRUE)))
cumsum <- data.frame(cumsum)
# calculate the proportion of each ASV per day and piece
dflongsub <- dflongsub %>%
  group_by(slice) %>%
  mutate(total_value = sum(mean, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(proportion = mean / total_value)
# merge cumsum with dflong
dflongsub <- merge(dflongsub, cumsum, by = c("slice"))

##########################################################################
# import data if isolates
frequency <- read.table("./03-Routputfiles/F_blastn_fun_seqvsiso_frequencies.txt", sep = "\t", header = TRUE)
colnames(frequency)[1] <- "sseqid"
dflongsub <- plyr::join(dflongsub, frequency, by = "sseqid")
# create a column in dflongsub that contains a 1 if the frequency is not na. else a 0.5
dflongsub$iso <- ifelse(is.na(dflongsub$Freq), "noiso", "iso")
##########################################################################
# plot taxaplots
##########################################################################
# Function to create custom labels
legendlabel<-unique(dflong$generataxaplot)
legendlabel<-legendlabel[c(2:(length(legendlabel)),1)]
llsplit1 <- ifelse(grepl("uncl.", legendlabel), "uncl.", "")
llsplit2 <- legendlabel
llsplit2 <-gsub("uncl. ", "", llsplit2)
legendlabel<-ape::mixedFontLabel(llsplit1, llsplit2,sep=" ",italic=2,always.upright = c("Others"))

  p <- ggplot(dflongsub, aes(factor(days), proportion * total, fill = reorder(generataxaplot,index)),
               col = "lightgrey") +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = colors) +
    geom_col(color = ifelse(dflongsub$iso == "iso", "black", "lightgrey"),
             linewidth = ifelse(dflongsub$iso == "iso", 0.5, 0.0001),
             alpha = ifelse(dflongsub$iso == "iso", 1, 1)) +
    geom_col(color = ifelse(dflongsub$iso == "iso", "black", "lightgrey"),
             linewidth = ifelse(dflongsub$iso == "iso", 0.5, 0.0001),
             linetype = ifelse(dflongsub$iso == "iso", "solid", "blank"),
             alpha = ifelse(dflongsub$iso == "iso", 1, 0)) +
    guides(fill = guide_legend(title = "Genera", ncol = 1)) +
    xlab("") +
    ylab("") +
    theme(legend.position = "none") +
    scale_y_continuous(labels = scales::label_number(accuracy = 0.1)) 
  plot_list[[i]] <- p
}
p <- ggplot(dflongsub, aes(factor(days), proportion * total, fill = reorder(generataxaplot,index)),
               col = "lightgrey") +
    geom_bar(stat = "identity") +
    theme_bw() +
    scale_fill_manual(values = colors, labels = legendlabel) +
    geom_col(color = ifelse(dflongsub$iso == "iso", "black", "lightgrey"),
             linewidth = ifelse(dflongsub$iso == "iso", 0.5, 0.0001),
             alpha = ifelse(dflongsub$iso == "iso", 1, 1)) +
    geom_col(color = ifelse(dflongsub$iso == "iso", "black", "lightgrey"),
             linewidth = ifelse(dflongsub$iso == "iso", 0.5, 0.0001),
             linetype = ifelse(dflongsub$iso == "iso", "solid", "blank"),
             alpha = ifelse(dflongsub$iso == "iso", 1, 0)) +
    guides(fill = guide_legend(title = "Genera", ncol = 4)) +
    xlab("") +
    ylab("")
legend <- get_legend(p)
p1 <- plot_grid(legend)
plot_list[[10]] <- p1

combined_plot1 <- plot_grid(plotlist = plot_list[1:9], ncol = 3, nrow = 3, labels = letters[1:9], label_size = 12, hjust = -1, vjust = 1.5)
# We need to create a nested grid for the last two plots
bottom_left_plot <- plot_grid(plot_list[[10]], align = 'v')

# Combine the last two plots into one row with different widths
combined_plot2 <- plot_grid(bottom_left_plot, ncol = 1)

combined_plot <- plot_grid(combined_plot1, combined_plot2, nrow = 2, rel_heights = c(3, 1))

# Create an empty drawing layer with the combined plot
# Add margins around the grid to ensure space for labels
final_grid_with_margins <- insert_xaxis_grob(combined_plot, grob = grid::nullGrob(), position = "top", height = grid::unit(0.05, "null"))
final_grid_with_margins <- insert_yaxis_grob(final_grid_with_margins, grob = grid::nullGrob(), position = "left", width = grid::unit(0.05, "null"))

# Now add overarching labels using ggdraw()
final_plot <- ggdraw(final_grid_with_margins) +
  draw_label("days", x = 0.5, y = 0.23, hjust = 0.5, vjust = 0) +
  draw_label("log10(estimated 18S rRNA gene copies)", x = 0.03, y = 0.65, angle = 90, hjust = 0.5, vjust = 1)


#print(final_plot)
# save plot as svg
ggsave("./03-Routputfiles/G_Taxaplots/G2_fun_individual_days.svg", final_plot,
       width = 21, height = 29.5, units = "cm")

ggsave("./03-Routputfiles/G_Taxaplots/G2_fun_individual_days.png", final_plot,
       width = 21, height = 29.5, units = "cm")

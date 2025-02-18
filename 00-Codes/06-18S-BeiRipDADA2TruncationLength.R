##############################################################################
# load libraries
library(stringr)
library(tidyr)
library(ggplot2)
##############################################################################
prodir <- getwd()
# set working directory
setwd("./05-qiime/HTS_18S/05-trimmingoptions")
##############################################################################
# print library versions used and save as txt file
session_info <- capture.output(sessionInfo())
sink("./06-18S-SessionInfo.txt")
cat(session_info, sep = "\n")
sink()
##############################################################################
# generate a file list
filelist <- dir(".", pattern = "\\.tsv$", full.names = TRUE, recursive = TRUE)
# create an empty data frame
sumstats <- data.frame("mean_input_passed" = NA, "median_input_passed" = NA,
                       "mean_input_merged" = NA, "median_input_merged" = NA,
                       "mean_non_chimeric" = NA, "median_non_chimeric" = NA,
                       "trunc" = NA)
# read stats from the files in the file list and put them into the data frame
k <- 1
for (file in filelist){
  data <- read.table(file, skip = 2, sep = "\t")
  names(data) <- read.table(file, nrow = 1, sep = "\t", comment.char = "")[1, ]
  stats <- data.frame("mean_input_passed" = NA, "median_input_passed" = NA,
                      "mean_input_merged" = NA, "median_input_merged" = NA,
                      "mean_non_chimeric" = NA, "median_non_chimeric" = NA)
  j <- 1
  for (i in c(4, 7, 9)){
    stats[1, j] <- mean(data[, i])
    stats[1, j + 1] <- median(data[, i])
    j <- j + 2
  }
  name <- str_remove(file, "./stats-dada2_")
  name <- str_remove(name, "/metadata.tsv")
  sumstats[k, ] <- stats[1, ]
  sumstats[k, 7] <- name
  k <- k + 1
}
# plot the data
datalong <- pivot_longer(sumstats, 1:6, names_to = "feature",
                         values_to = "percentage")
p <- ggplot(datalong, aes(trunc, percentage)) +
  geom_boxplot() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1,
                                   vjust = 0.5, size = 10))
#save plot as tiff
ggsave("./dada2-trimming-stats-comparison.tiff", plot = p,
       width = 10, height = 10, units = "cm", dpi = 300)
##############################################################################
setwd(prodir)
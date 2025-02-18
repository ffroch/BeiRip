library(tidyverse)
library(ggplot2)
# Import data from qPCR machine
plate1 <- read.csv("./02-rawdata/BeiRipSamples-Plate1.csv", header = TRUE, skip = 17, sep = ";", dec = ",")
plate2 <- read.csv("./02-rawdata/BeiRipSamples-Plate2.csv", header = TRUE, skip = 17, sep = ";", dec = ",")
plate3 <- read.csv("./02-rawdata/BeiRipSamples-Plate3.csv", header = TRUE, skip = 17, sep = ";", dec = ",")

# combine the plates
plate <- rbind(plate1, plate2, plate3)

# keep only the columns that are needed
platered <- plate %>%
  select("Sample.name", "Sample.type", "Ct", "Mean.Ct", "Conc..Std.", "Mean.Conc.",
         "Std.Dev..Ct", "Std.Dev..Mean.Conc.")
# rename the columns
colnames(platered) <- c("SampleName", "SampleType", "Ct", "MeanCt", "ConcStd", "MeanConc",
                        "StdDevCt", "StdDevMeanConc")

# keep only samples
platered <- platered %>%
  filter(SampleType == "Unknown")
# remove samples with NK or Mock in sample name
plateredsam <- platered %>%
  filter(!grepl("NK", SampleName)) %>%
  filter(!grepl("Mock", SampleName)) %>%
  filter(SampleName != "")

platredcon <- platered %>%
  filter(!grepl("B", SampleName))

# add a column named piece with the first 2 letters of the sample name
plateredsam$piece <- substr(plateredsam$SampleName, 1, 2)
# add a column named slice with the last two letters of the sample name
plateredsam$slice <- substr(plateredsam$SampleName, 4, 5)
# make slice numeric
plateredsam$slice <- as.numeric(plateredsam$slice)
# save plateredsam as rds
saveRDS(plateredsam, "./03-Routputfiles/qpcrsampleresults.rds")

ggplot(plateredsam, aes(x = slice, y = log10(MeanConc), group = piece)) +
  geom_point(aes(color = piece)) +
  geom_smooth(method = "loess", na.rm = TRUE, se = FALSE) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Slice", y = "Mean copy numbers") +
  facet_wrap(~piece, ncol = 3)

ggplot(plateredsam, aes(x = slice, y = log10(MeanConc), group = 1)) +
  geom_point() +
  geom_smooth(method = "loess", na.rm = TRUE, se = FALSE) +
  theme_bw() +
  theme(legend.position = "none") +
  labs(x = "Slice", y = "Mean copy numbers") 

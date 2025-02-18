################################################################################
# load libraries
library(Biostrings)
library(DECIPHER)
library(seqinr)
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
##########################################################################
# PART 2: Primer design based on isolates
##########################################################################
# This is an alternative approach to the NCBI blast primer approach.
# I used the sequences of the 18S sequencing data, with main target 
# df9111b84ad2c6dd3e644f52434bbbd7
# with a later trimming 3d4847a1a754b6437afe9c164583ef25
##########################################################################
##########################################################################
# data import and preparation
##########################################################################
# import 18S ASVs as DNA stringset (raw and filtered asvs)
asvs <- readDNAStringSet(paste0("./05-qiime/HTS_18S/04-asvtable/",
                                "03-final/exportedrepseqs.final/",
                                "dna-sequences.fasta"))
asvsraw <- readDNAStringSet(paste0("./05-qiime/HTS_18S/04-asvtable/",
                                   "01-raw/exported-repseqs/",
                                   "dna-sequences.fasta"))
# get the sequence of the ASV named df9111b84ad2c6dd3e644f52434bbbd7; with a later trimming 3d4847a1a754b6437afe9c164583ef25
target <- asvs[names(asvs) == "3d4847a1a754b6437afe9c164583ef25"]
targetc <- as.character(target)
# split target into 20 bp long strings and match it with all sequences in asvs
empty <- data.frame()
for (i in (1: nchar(targetc))) {
    target1 <- substr(targetc, i, i + 19)
matches <- vmatchPattern(target1, asvsraw, fixed = TRUE)
matching_asvs <- asvsraw[unlist(lapply(matches, length)) > 0]
empty[i, 1] <- paste0(i)
empty[i, 2] <- length(matching_asvs)
}
empty$V1 <- as.numeric(empty$V1)
# plot the number of ASVs matching the target sequence
ggplot(empty, aes(V1, log10(V2))) +
    geom_point() +
    geom_line() +
    theme_bw() +
    xlab("position in query sequence") +
    ylab("number of matching ASVs")
ggsave("./03-Routputfiles/Z_18SqPCR/matchingASVsOn20mers.png", width = 20, height = 15)
###############################################################################
# add information of the primers
primers <- c("AACCGAGCCTTTCCTTCTGG", "GCAAAGGCCTGCTTTGAACA", #PP1
             "AGTTGAACCTTGGGCTTGGT", "AGAAGGAAAGGCTCGGTTGG", #PP2
             "GGTCCGCTTTATGGCGAGTA", "GTTCGCCAAACACCACAAGG", #PP3
             "TGTTTGGCGAACCAGGACTT", "GAATACTGATGCCCCCGACC", #PP4
             "CCAACCGAGCCTTTCCTTCT", "AGCAAAGGCCTGCTTTGAAC", #PP5
             "GGTCGGGGGCATCAGTATTC", "TTGGCAAATGCTTTCGCAGT", #PP6
             "AGCAGGCCTTTGCTCGAATA", "TGAATACTGATGCCCCCGAC", #PP7
             "TGAACCTTGGGCTTGGTTGG", "CCAGAAGGAAAGGCTCGGTT", #PP8
             "CTGGCTAACCATTCGCCCTT", "AGGCCTGCTTTGAACACTCT", #PP10
             "TTCTGGCTAACCATTCGCCC", "ACTGAATACTGATGCCCCCG") #PP9

primersfwd <- primers[seq(1, length(primers), 2)]
primersrev <- primers[seq(2, length(primers), 2)]
primersrev <- reverseComplement(DNAStringSet(primersrev)) %>% as.character()
# get position for each primer in primers on target
primerpos <- data.frame(V1 = c("PP1", "PP2", "PP3", "PP4", "PP5", "PP6", "PP7", "PP8", "PP10", "PP9"))
for (i in seq_along(primersfwd)) {
    primer <- primersfwd[i]
    matchesfwd <- vmatchPattern(primer, target, fixed = TRUE) %>% data.frame()
    primerpos[i, 2] <- matchesfwd$start
    primerpos[i, 3] <- matchesfwd$end
}
for (i in seq_along(primersrev)) {
    primer <- primersrev[i]
    matchesrev <- vmatchPattern(primer, target, fixed = TRUE) %>% data.frame()
    primerpos[i, 4] <- matchesrev$start
    primerpos[i, 5] <- matchesrev$end
}
primerpos$V6 <- seq(3.48, 2.58, -0.1)
primerpos$V7 <- seq(3.42, 2.52, -0.1)
primerpos$V8 <- seq(3.45, 2.55, -0.1)

ggplot(empty, aes(V1, log10(V2))) +
    geom_point() +
    geom_line() +
    theme_bw() +
    xlab("position in query sequence") +
    ylab("log10(number of matching ASVs)") +
    geom_rect(data = primerpos[1:9, ], aes(xmin = V2, xmax = V3, ymin = V7, ymax = V6), fill = "#154360", alpha = 0.5, inherit.aes = FALSE) +
    geom_rect(data = primerpos[1:9, ], aes(xmin = V4, xmax = V5, ymin = V7, ymax = V6), fill = "#154360", alpha = 0.5, inherit.aes = FALSE) +
    geom_segment(data = primerpos[1:9, ], aes(x = V3, xend = V4, y = V8, yend = V8), color = "#154360", inherit.aes = FALSE) +
    geom_text(data  = primerpos[1:9, ], aes(x = (V3 + V4)/2, y = V8, label = V1), vjust = -0.5, hjust = 0.5, size = 5, color = "black", inherit.aes = FALSE) +
    geom_rect(data = primerpos[10, ], aes(xmin = V2, xmax = V3, ymin = 0, ymax = V6), fill = "darkred", alpha = 0.5, inherit.aes = FALSE) +
    geom_rect(data = primerpos[10, ], aes(xmin = V4, xmax = V5, ymin = 0, ymax = V6), fill = "darkred", alpha = 0.5, inherit.aes = FALSE) +
    geom_segment(data = primerpos[10, ], aes(x = V3, xend = V4, y = V8, yend = V8), color = "darkred", inherit.aes = FALSE) +
    geom_text(data  = primerpos[10, ], aes(x = (V3 + V4)/2, y = V8, label = V1), vjust = -0.5, hjust = 0.5, size = 5, color = "black", inherit.aes = FALSE)
ggsave("./03-Routputfiles/Z_18SqPCR/matchingASVsOn20mers_plusprimers.png", width = 20, height = 15)
################################################################################
# filter empty for the minimum number of matches in V2
empty2 <- empty[empty$V2 == min(empty$V2), ]
# get select the positions in the target sequence listed in empty2
target2 <- substr(target, empty2$V1[1], empty2$V1[nrow(empty2)] + 19)
# list asvs that match target2
matches2 <- vmatchPattern(target2, asvsraw, fixed = TRUE)
matching_asvs2 <- asvsraw[unlist(lapply(matches2, length)) > 0]
matching_asvs2
################################################################################
# split target into 20bp long strings and match it with all sequences in matching_asvs2
empty3 <- data.frame()
for (i in (1: nchar(target))) {
    target1 <- substr(target, i, i + 19)
matches3 <- vmatchPattern(target1, matching_asvs2, fixed = TRUE)
matching_asvs3 <- matching_asvs2[unlist(lapply(matches3, length)) > 0]
empty3[i, 1] <- paste0(i)
empty3[i, 2] <- length(matching_asvs3)
}
empty3$V1 <- as.numeric(empty3$V1)
# plot the number of ASVs matching the target sequence
ggplot(empty3, aes(V1, V2)) +
    geom_point() +
    geom_line() +
    theme_bw() +
    xlab("position in query sequence") +
    ylab("number of matching ASVs")
################################################################################
# filter empty3 for the minimum number of matches in V2
empty4 <- empty3[empty3$V2 == min(empty3$V2), ]
# get select the positions in the target sequence listed in empty4
target3 <- substr(target, empty4$V1[1], empty4$V1[nrow(empty4)] + 19)
# list asvs that match target3
matches4 <- vmatchPattern(target3, matching_asvs3, fixed = TRUE)
matching_asvs4 <- matching_asvs3[unlist(lapply(matches4, length)) > 0]
matching_asvs4

# selected primer pair = PP9
# check the primer pair on the target sequences

primerpair9 <- c("TTCTGGCTAACCATTCGCCC", "ACTGAATACTGATGCCCCCG")
primer9fwd <- primerpair9[1]
primer9rev <- reverseComplement(DNAStringSet(primerpair9[2])) %>% as.character()

# 0 mismatches
allmatchesfwd <- vmatchPattern(primer9fwd, asvsraw, fixed = TRUE)
matchingasvsfwd <- asvsraw[unlist(lapply(allmatchesfwd, length)) > 0]
allmatches <- vmatchPattern(primer9rev, matchingasvsfwd, fixed = TRUE)
matchingasvsall0mm <- matchingasvsfwd[unlist(lapply(allmatches, length)) > 0]

# 1 mismatch per primer
allmatchesfwd <- vmatchPattern(primer9fwd, asvsraw, fixed = TRUE, max.mismatch = 1)
matchingasvsfwd <- asvsraw[unlist(lapply(allmatchesfwd, length)) > 0]
allmatches <- vmatchPattern(primer9rev, matchingasvsfwd, fixed = TRUE, max.mismatch = 1)
matchingasvsall1mm <- matchingasvsfwd[unlist(lapply(allmatches, length)) > 0]

saveRDS(names(matchingasvsall0mm), "./03-Routputfiles/Z_18SqPCR/matchingasvsall0mm.rds")
saveRDS(names(matchingasvsall1mm), "./03-Routputfiles/Z_18SqPCR/matchingasvsall1mm.rds")

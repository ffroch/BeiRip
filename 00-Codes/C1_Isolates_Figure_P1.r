###############################################################################
### FIGURE 1 ###
#Figure 1 shows the following:
#1. A tree with one representative of each Genus based on the 16S/18S rRNA
# gene sequences, representative is isolate with the longest sequence in the
# most prevalent Cluster within a Genus
#2. Additional information about the isolates and their diversity within a Genus
##2a. column for each batch with points
##2b. points are sized by the grouped number of 16S/18S rRNA gene genotypes
##2c. points are coloured by the last slice the isolate was found in
#only sequences >1000bp with taxonomy >=  2 matching databases were used
###############################################################################
### THIS IS PART 1 - DATA PREPARATION ###
###############################################################################
library(ggplot2)
library(DECIPHER)
library(dplyr)
library(Biostrings)
library(stringr)
###############################################################################
#### FIGURE 1A - BACTERIA ####
### Load data ###
#load stringset of bacteria
stringsetbac <- readDNAStringSet("./03-Routputfiles/B_SangerSeq/B1_seqstrimclean_35_15_100bp_bac.fa", "fasta")
#load bacterial classification data
lodibac <- readRDS("./03-Routputfiles/B_SangerSeq/B3_ListOfDoneIsolatesRevised.rds")
#check if ssbac1000 names and lodibac ids are in the same order
identical(names(stringsetbac), lodibac$ids)
## keep only isolates with taxonomy hits >=  2 databases
#remove all isolates in lodibac, that does not have a string in newGenus
lodibacred <- lodibac[!is.na(lodibac$newGenus), ]
#remove all isolates in stringsetbac, that are not present in lodibacred$ids
stringsetbacred <- stringsetbac[names(stringsetbac) %in% lodibacred$ids]
## keep only isolates with >=  1000bp sequence length
#remove all isolates in stringsetbac that are shorter than 1000bp and save it as ssbac1000
ssbac1000 <- stringsetbacred[width(stringsetbacred) >=  1000]
#count seqlength
lodibacred$seqlength <- str_count(lodibacred$row.names)
#remove all isolates in lodibac, that are shorter than 1000bp
lodibacred <- lodibacred[lodibacred$seqlength >=  1000, ]

### Clusterize data based on sequence similarity ###
## clusterize sequences
clusters <- Clusterize(ssbac1000, cutoff = 0.001, penalizeGapLetterMatches = TRUE)
### add cluster information to lodibacred
#check if clusters names and lodibacred ids are in the same order
identical(row.names(clusters), lodibacred$ids)
#add cluster information to lodibacred
lodibacred$cluster <- clusters$cluster
## count isolates based on newGenus and Cluster
clugen <- lodibacred %>%
  group_by(newGenus, cluster) %>%
  summarise(count = n())
## count clusters per Genus
div <- clugen %>%
  ungroup() %>%
  group_by(newGenus) %>%
  summarise(diversity = n())
## add clugen$count to lodibacred based on cluster
lodibacredc <- plyr::join(lodibacred, clugen, by = c("newGenus", "cluster"), type = "left")
## add diversity to lodibacredc based on newGenus
lodibacredc <- plyr::join(lodibacredc, div, by = "newGenus", type = "left")
## keep only one representative per Genus, from the most prevalent cluster with the longest sequence
# order lodibacredc by newGenus, count and seqlength
lodibacredc <- lodibacredc[order(lodibacredc$newGenus,
                                 lodibacredc$count,
                                 lodibacredc$seqlength,
                                 decreasing = TRUE), ]
# save lodibacredc as RDS for Part3
saveRDS(lodibacredc, "./03-Routputfiles/C_IsolateTree/C1_IsolatesUsedForFig1bac.RDS")
# keep only the first entry for each newGenus
representatives <- lodibacredc[!duplicated(lodibacredc$newGenus), ]
#remove all isolates in stringsetbac, that are not present in representative$ids
ssbrep <- stringsetbac[names(stringsetbac) %in% representatives$ids]
## add metadata from representatives to ssbrep
# make ids to rownames
row.names(representatives) <- representatives$ids
# sort representatives by ids
representatives <- representatives[order(representatives$ids), ]
#check if ssbrep names and representatives ids are in the same order
identical(names(ssbrep), row.names(representatives))
# add metadata from representatives to ssbrep
metadata(ssbrep) <- representatives

# save ssbrep as RDS
saveRDS(ssbrep, "./03-Routputfiles/C_IsolateTree/C1_RepresentativeIsolates_bac.RDS")
# write ssbrep as stringset to file



writeXStringSet(ssbrep, "./03-Routputfiles/C_IsolateTree/C1_RepresentativeIsolates_bac.fa")
###############################################################################
#### FIGURE 1A - FUNGI ####
# only sequences >1000bp with taxonomy >=  2 matching databases were used
# for fungi, isolates with lower taxonomy level than Genus were used if Genus was not available
### Load data ###
#load stringset of fungi
stringsetfun <- readDNAStringSet("./03-Routputfiles/B_SangerSeq/B2_seqstrimclean_35_15_100bp_fun.fa", "fasta")
#load fungal classification data with classification on Genus level at minBoot 100 using SILVA DB
lodifun <- read.table("./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_silva_genus_CI100.txt",
                      sep = ";", row.names = NULL)
lodifun <- data.frame(lodifun)
#get IDs from stringsetfun
idsfun <- names(stringsetfun)
#get seqlength from stringsetfun
lengthfun <- width(stringsetfun)
#make a dataframe with IDs, seqlength and classification on Genus level at minBoot 100
# using SILVA DB
lodifun <- cbind(idsfun, lengthfun, lodifun)
#remove isolates with NA in kingdom
lodifun <- lodifun[!is.na(lodifun$Kingdom), ]
#add new column named newtax and fill it with NA
lodifun$newtax <- NA
#find the column with the last NA entry counting from left
lodifun$firstNA <- apply(lodifun, 1, function(x) min(which(is.na(x))))
#give back a list of the columns that are indext in dffun2lowres$lastNA
lodifun$lastNAcol <- apply(lodifun, 1, function(x) names(x)[min(which(is.na(x)) - 1)])
#add column to dffun2lowres named newtax that contains "unclassified "
# and the content of the column that is indexed in lastNAcol
lodifun$newtax <- paste0("unclassified ",
                         apply(lodifun, 1, function(x) x[min(which(is.na(x)) - 1)]))
#replace newtax entry by Genus if lastNAcol is Genus
lodifun$newtax[lodifun$lastNAcol == "Genus"] <- lodifun$Genus[lodifun$lastNAcol == "Genus"]

### Filter data ###
#remove all isolates in lodifun with seqlength <900bp
lodifunred <- lodifun[lodifun$lengthfun >=  900, ]
#remove all isolates in stringsetfun that are shorter than 900bp and save it as ssfun900
ssfun900 <- stringsetfun[width(stringsetfun) >=  900]

### Clusterize data based on sequence similarity ###
## clusterize sequences
clusters <- Clusterize(ssfun900, cutoff = 0.001, penalizeGapLetterMatches = TRUE)
### add cluster information to lodifunred
#check if clusters names and lodifunred ids are in the same order
identical(row.names(clusters), lodifunred$idsfun)
#add cluster information to lodifunred
lodifunred$cluster <- clusters$cluster
## count isolates based on newtax and Cluster
clugenf <- lodifunred %>%
  group_by(newtax, cluster) %>%
  summarise(count = n())
## count clusters per Genus
divf <- clugenf %>%
  ungroup() %>%
  group_by(newtax) %>%
  summarise(diversity = n())
## add clugenf$count to lodifunred based on cluster
lodifunredc <- plyr::join(lodifunred, clugenf, by = c("newtax", "cluster"), type = "left")
## add diversity to lodifunredc based on newtax
lodifunredc <- plyr::join(lodifunredc, divf, by = "newtax", type = "left")
# save lodifunredc as RDS for Part3
saveRDS(lodifunredc, "./03-Routputfiles/C_IsolateTree/C1_IsolatesUsedForFig1fun.RDS")
## keep only one representative per Genus, from the most prevalent cluster with the longest sequence
# order lodifunredc by newtax, count and lengthfun
lodifunredc <- lodifunredc[order(lodifunredc$newtax,
                                 lodifunredc$count,
                                 lodifunredc$lengthfun,
                                 decreasing = TRUE), ]
# keep only the first entry for each newtax
representativesf <- lodifunredc[!duplicated(lodifunredc$newtax), ]
#remove all isolates in stringsetfun, that are not present in representativesf$idsfun
ssfrep <- stringsetfun[names(stringsetfun) %in% representativesf$idsfun]
## add metadata from representativesf to ssfrep
# make ids to rownames
row.names(representativesf) <- representativesf$idsfun
# sort representativesf by ids
representativesf <- representativesf[order(representativesf$idsfun), ]
#check if ssfrep names and representativesf ids are in the same order
identical(names(ssfrep), row.names(representativesf))
# add metadata from representativesf to ssfrep
metadata(ssfrep) <- representativesf

# save ssfrep as RDS
saveRDS(ssfrep, "./03-Routputfiles/C_IsolateTree/C1_RepresentativeIsolates_fun.RDS")

################################################################################
# Title: B3_GenerateListOfDoneIsolates.R
# Description: This script generates the following tables for the manuscript:
# B3_ListOfDoneIsolates.txt
# B3_classificationstats_bac.tsv (how many isolates with matching and
# mismatching classifications at minboot 100)
# B3_ListOfDoneIsolates_mismatches.txt (list of isolates, with mismatching
# classification at minboot 100)
# B3_mismatchedGenera100CI_bac.rds (table with isolate counts per mismatch
# "group" at minboot 90)
# B3_ListOfDoneIsolates_mismatches_shortseqs (EXKURS - checking for larger
# number of mismatches in short sequences)
# B3_ListOfDoneIsolatesRevised.rds (revised list with Consensus taxonomy as
# newGenus or new Genusweak
# B3_Isolates100CI.rds (Table SM1 - List of all bacterial isolates with their
# classification on Genus level at minBoot 100 with matches in all three
# databases or in two databases + one unclassified
# B3_mismatchedGenera90CI.txt (list of isolates, with mismatching classification
# at minboot 90)
# B3_mismatched Genera90CI_bac.rds (table with isolate counts per mismatch
# "group" at minboot 90)
# B3_Isolates90CI.rds (Table SM2 - List of all bacterial isolates with their
# classification on Genus level at minBoot 90 with matches in all three
# databases or in two databases + one unclassified
# B3_mismatchedGenera.rds (Table SM3 - List of all bacterial isolates with
#their classification on Genus level and the classification on Genus level
# from GTDB and RDP and the classification on Species level from GTDB and RDP
# B3_IsolatesFungi.rds (Table SM4 - List of all fungal isolates, with CI100
# or CI90)
# B3_IsolatesSpecies.rds (Table SM5 - List of all bacterial isolates with
# CI100 on species level
################################################################################
library(Biostrings)
library(dplyr)
###############################################################################
#load stringset that contains all bacterial isolates incl. their
# IDs and seqlength
stringsetbac <-
  readDNAStringSet("./03-Routputfiles/B_SangerSeq/B1_seqstrimclean_35_15_100bp_bac.fa",
                   "fasta")
#load bacterial classification data with classification on Genus level at
# minBoot 100 using SILVA DB
taxtabbac <-
  read.table("./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_silva_genus_CI100.txt",
             sep = ";", row.names  =  NULL)
taxtabbac1 <- data.frame(taxtabbac)
#get IDs from stringsetbac
idsbac <- names(stringsetbac)
#get seqlength from stringsetbac
lengthbac <- width(stringsetbac)
#make a dataframe with IDs, seqlength and classification on Genus level at
# minBoot 100 using SILVA DB
dfbac <- cbind(idsbac, lengthbac, taxtabbac1)
#load the other classification data
silva90 <- read.table("./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_silva_genus_CI90.txt",
                      sep = ";", row.names  =  NULL)
gtdb100 <- read.table("./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_GTDB_species_CI100.txt",
                      sep = ";", row.names  =  NULL)
gtdb90 <- read.table("./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_GTDB_species_CI90.txt",
                     sep = ";", row.names  =  NULL)
rdp100 <- read.table("./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_RDP_species_CI100.txt",
                     sep = ";", row.names  =  NULL)
rdp90 <- read.table("./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_RDP_species_CI90.txt",
                    sep = ";", row.names  =  NULL)
#add the other classification data to dfbac
addon <- data.frame("silva90" = silva90$Genus,
                    "gtdb100genus" = gtdb100$Genus,
                    "gtdb100species" = gtdb100$Species,
                    "gtdb90genus" = gtdb90$Genus,
                    "gtdb90species" = gtdb90$Species,
                    "rdp100genus" = rdp100$Genus,
                    "rdp100species" = rdp100$Species,
                    "rdp90genus" = rdp90$Genus,
                    "rdp90species" = rdp90$Species)
dfbac <- cbind(dfbac, addon)
###############################################################################
#load stringset that contains all fungal isolates incl. their IDs and seqlength
stringsetfun <- readDNAStringSet("./03-Routputfiles/B_SangerSeq/B2_seqstrimclean_35_15_100bp_fun.fa", "fasta")
#load fungal classification data with classification on Genus level at minBoot 100 using SILVA DB
taxtabfun <- read.table("./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_silva_genus_CI100.txt",
                        sep = ";", row.names = NULL)
taxtabfun1 <- data.frame(taxtabfun)
#get IDs from stringsetfun
idsfun <- names(stringsetfun)
#get seqlength from stringsetfun
lengthfun <- width(stringsetfun)
#make a dataframe with IDs, seqlength and classification on Genus level at
# minBoot 100 using SILVA DB
dffun <- cbind(idsfun, lengthfun, taxtabfun1)
#load the other classification data
silva90fun <- read.table("./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_silva_genus_CI90.txt",
                         sep = ";", row.names  =  NULL)
#add the other classification data to dffun
addonfun <- data.frame("silva90" = silva90fun$Genus)
dffun <- cbind(dffun, addonfun)
###############################################################################
#make a dataframe with all isolates and save it as B3_ListOfDoneIsolates.txt
#add empty columns for gtdb100genus, gtdb100species, gtdb90genus,
# gtdb90species, rdp100genus, rdp100species, rdp90genus, rdp90species
dffun$gtdb100genus <- NA
dffun$gtdb100species <- NA
dffun$gtdb90genus <- NA
dffun$gtdb90species <- NA
dffun$rdp100genus <- NA
dffun$rdp100species <- NA
dffun$rdp90genus <- NA
dffun$rdp90species <- NA
#if Phlyum in dffun is NA, fill with Phylum from silva90fun
dffun$Phylum[is.na(dffun$Phylum)] <- dffun$silva90[is.na(dffun$Phylum)]
#if Class in dffun is NA, fill with Class from silva90fun
dffun$Class[is.na(dffun$Class)] <- dffun$silva90[is.na(dffun$Class)]
#if Order in dffun is NA, fill with Order from silva90fun
dffun$Order[is.na(dffun$Order)] <- dffun$silva90[is.na(dffun$Order)]
#if Family in dffun is NA, fill with Family from silva90fun
dffun$Family[is.na(dffun$Family)] <- dffun$silva90[is.na(dffun$Family)]

colnames(dfbac)[1:2] <- c("ids", "seqlength")
colnames(dffun)[1:2] <- c("ids", "seqlength")
dfall <- rbind(dfbac, dffun)
write.table(dfall, "./03-Routputfiles/B_SangerSeq/B3_ListOfDoneIsolates.txt", sep = ";")
################################################################################
#statistic for dfbac to check the classification metrics
#gtdb has addons to the genusname like _A or _E
dfbac1 <- dfbac
#remove addons
dfbac1$gtdb100genus <- gsub("_.*", "", dfbac$gtdb100genus)
#silva (Genus) has addons to genusname in form of numbers
#remove addons
dfbac1$Genus <- gsub("[0-9]", "", dfbac$Genus)
#make an empty dataframe with colnames "parameter" and "value"
sdfbactax <- data.frame("parameter" = 0, "value" = 0)
#count number of isolates
sdfbactax[1, ] <- c("number of isolates", nrow(dfbac1))
#count number of isolates per genus that are not NA
sdfbactax[2, ] <- c("number of isolates in silva100 that are not NA",
                    nrow(dfbac1[!is.na(dfbac1$Genus), ]))
#count number of isolates per gtdb100genus that are not NA
sdfbactax[3, ] <- c("number of isolates in gtdb100genus that are not NA",
                    nrow(dfbac1[!is.na(dfbac1$gtdb100genus), ]))
#count number of isolates per gtdb100species that are not NA
sdfbactax[4, ] <- c("number of isolates in gtdb100species that are not NA",
                    nrow(dfbac1[!is.na(dfbac1$gtdb100species), ]))
#count number of isolates per rdp100genus that are not NA
sdfbactax[5, ] <- c("number of isolates in rdp100genus that are not NA",
                    nrow(dfbac1[!is.na(dfbac1$rdp100genus), ]))
#count number of isolates per rdp100species that are not NA
sdfbactax[6, ] <- c("number of isolates in rdp100species that are not NA",
                    nrow(dfbac1[!is.na(dfbac1$rdp100species), ]))
print(sdfbactax)
################################################################################
#make a dataframe as basis for Tables S1 - S3
#Part 1: minboot 100
################################################################################
#prep Part1:
#make an empty dataframe with colnames "parameter" and "value"
df2 <- data.frame("parameter" = 0, "value" = 0,
                  "nNA_matching" = 0, "nNA_not_matching" = 0)
#count NA per row and add to dfbac1
dfbac1$NAcount <- rowSums(is.na(dfbac1[, c(9, 11, 15)]))
#filter isolates per sum of NA in Genus, gtdb100genus and rdp100genus
dfbac10 <- dfbac1[dfbac1$NAcount == 0, ]
dfbac11 <- dfbac1[dfbac1$NAcount == 1, ]
dfbac12 <- dfbac1[dfbac1$NAcount == 2, ]
dfbac13 <- dfbac1[dfbac1$NAcount == 3, ]
################################################################################
#stats Part1:
#counts of 0 NA isolates
df2[1, ] <- c("number of isolates with 0 NA", nrow(dfbac10),
              nrow(dfbac10[dfbac10$Genus == dfbac10$gtdb100genus & dfbac10$Genus == dfbac10$rdp100genus, ]),
              nrow(dfbac10[dfbac10$Genus != dfbac10$gtdb100genus | dfbac10$Genus != dfbac10$rdp100genus, ]))
#counts of 1 NA isolates
df2[2, ] <- c("number of isolates with 1 NA", nrow(dfbac11),
              nrow(dfbac11[dfbac11$Genus == dfbac11$gtdb100genus & dfbac11$Genus == dfbac11$rdp100genus, ]),
              nrow(dfbac11) - nrow(dfbac11[dfbac11$Genus == dfbac11$gtdb100genus & dfbac11$Genus == dfbac11$rdp100genus, ]))
#counts of 2 NA isolates
df2[3, ] <- c("number of isolates with 2 NA", nrow(dfbac12), nrow(dfbac12), 0)
#counts of 3 NA isolates
df2[4, ] <- c("number of isolates with 3 NA", nrow(dfbac13), nrow(dfbac13), 0)
print(df2)
#save df2 as txt
saveRDS(df2, "./03-Routputfiles/B_SangerSeq/B3_classificationstats_bac.rds")
################################################################################
#mismatches Part1:
#show the isolates where the strings in Genus, gtdb100genus and rdp100genus are not the same
mismatch <- dfbac1[dfbac1$Genus != dfbac1$gtdb100genus | dfbac1$Genus != dfbac1$rdp100genus | dfbac1$gtdb100genus != dfbac1$rdp100genus, ]
#reduce mismatch to the columns ids, Family, Genus, gtdb100genus, rdp100genus
mismatch <- mismatch[, c(1, 4, 5, 6, 7, 8, 9, 11, 12, 15, 16)]
#remove isolates with only NA entries
mismatch <- mismatch[!is.na(mismatch$ids), ]
#save mismatch a txt
write.table(mismatch, "./03-Routputfiles/B_SangerSeq/B3_ListOfDoneIsolates_mismatches.txt", sep = ";")
################################################################################
#make a mismatch overview for Table S3
mismatchred <- mismatch
#paste Genus gtdb100genus and rdp100genus together separated by "/"
mismatchred$Genus <- paste(mismatchred$Genus, mismatchred$gtdb100genus,
                           mismatchred$rdp100genus, sep = "/")
#remove gtdb100genus and rdp100genus and ids
mismatchred <- mismatchred[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]
#count number of isolates per Genus
mismatchred <- mismatchred %>%
  group_by(Genus) %>%
  summarise(count = n())
mismatchred <- data.frame(mismatchred)
#save mismatchred as RDS to 03-Routputfiles/B_SangerSeq
saveRDS(mismatchred, "./03-Routputfiles/B_SangerSeq/B3_mismatchedGenera100CI_bac.rds")
################################################################################
#EXKURS
#statistic for bac with sequences with seqlength < 500
dfbac2 <- dfbac1[dfbac1$seqlength < 500, ]
#make an empty dataframe with colnames "parameter" and "value"
sdfbactaxshort <- data.frame("parameter" = 0, "value" = 0)
#count number of isolates
sdfbactaxshort[1, ] <- c("number of isolates with seqlength < 500", nrow(dfbac2))
#count number of isolates per genus that are not NA
sdfbactaxshort[2, ] <- c("number of isolates in silva100 that are not NA",
                         nrow(dfbac2[!is.na(dfbac2$Genus), ]))
#count number of isolates per gtdb100genus that are not NA
sdfbactaxshort[3, ] <- c("number of isolates in gtdb100genus that are not NA",
                         nrow(dfbac2[!is.na(dfbac2$gtdb100genus), ]))
#count number of isolates per gtdb100species that are not NA
sdfbactaxshort[4, ] <- c("number of isolates in gtdb100species that are not NA",
                         nrow(dfbac2[!is.na(dfbac2$gtdb100species), ]))
#count number of isolates per rdp100genus that are not NA
sdfbactaxshort[5, ] <- c("number of isolates in rdp100genus that are not NA",
                         nrow(dfbac2[!is.na(dfbac2$rdp100genus), ]))
#count number of isolates per rdp100species that are not NA
sdfbactaxshort[6, ] <- c("number of isolates in rdp100species that are not NA",
                         nrow(dfbac2[!is.na(dfbac2$rdp100species), ]))
print(sdfbactaxshort)
#make empty dataframe df3 with colnames "parameter", "value", "nNA_matching", "nNA_not_matching"
df3 <- data.frame("parameter" = 0, "value" = 0,
                  "nNA_matching" = 0, "nNA_not_matching" = 0)
#count NA per row and add to dfbac2
dfbac2$NAcount <- rowSums(is.na(dfbac2[, c(9, 11, 15)]))
#filter isolates per sum of NA in Genus, gtdb100genus and rdp100genus
dfbac20 <- dfbac2[dfbac2$NAcount == 0, ]
dfbac21 <- dfbac2[dfbac2$NAcount == 1, ]
dfbac22 <- dfbac2[dfbac2$NAcount == 2, ]
dfbac23 <- dfbac2[dfbac2$NAcount == 3, ]
#counts of 0 NA isolates
df3[1, ] <-
  c("number of isolates with 0 NA", nrow(dfbac20),
    nrow(dfbac20[dfbac20$Genus == dfbac20$gtdb100genus & dfbac20$Genus == dfbac20$rdp100genus, ]),
    nrow(dfbac20[dfbac20$Genus != dfbac20$gtdb100genus | dfbac20$Genus != dfbac20$rdp100genus, ]))
#counts of 1 NA isolates
df3[2, ] <- c("number of isolates with 1 NA", nrow(dfbac21),
              nrow(dfbac21[dfbac21$Genus == dfbac21$gtdb100genus & dfbac21$Genus == dfbac21$rdp100genus, ]),
              nrow(dfbac21) - nrow(dfbac21[dfbac21$Genus == dfbac21$gtdb100genus & dfbac21$Genus == dfbac21$rdp100genus, ]))
#counts of 2 NA isolates
df3[3, ] <- c("number of isolates with 2 NA",
              nrow(dfbac22), nrow(dfbac22), 0)
#counts of 3 NA isolates
df3[4, ] <- c("number of isolates with 3 NA",
              nrow(dfbac23), nrow(dfbac23), 0)
#show the isolates where the strings in Genus, gtdb100genus and rdp100genus are not the same
mismatchshort <- dfbac2[dfbac2$Genus != dfbac2$gtdb100genus | dfbac2$Genus != dfbac2$rdp100genus | dfbac2$gtdb100genus != dfbac2$rdp100genus, ]
#reduce mismatch to the columns ids, Family, Genus, gtdb100genus, rdp100genus
mismatchshort <- mismatchshort[, c(1, 4, 5, 6, 7, 8, 9, 11, 12, 15, 16)]
#remove isolates with only NA entries
mismatchshort <- mismatchshort[!is.na(mismatchshort$ids), ]
#save mismatchshort a txt
write.table(mismatchshort,
            "./03-Routputfiles/B_SangerSeq/B3_ListOfDoneIsolates_mismatches_shortseqs.txt",
            sep = ";")
#END EXKURS
###############################################################################################
#Generate structure of Table S1
#prepare a list for the supplements that contains counts per Genus
#generate make dfbac2 to dfbac3
dfbac3 <- dfbac1
#make an empty column named mismatching
dfbac3$mismatching <- 0
#fill mismatching with 1 if id is present in mismatch
dfbac3$mismatching[dfbac3$ids %in% mismatch$ids] <- 1
#make an empty column named newGenus
dfbac3$newGenus <- NA
#fill newGenus with Genus if NAcount == 0 or 1, if Genus is NA fill newGenus with
# gtdb100genus if gtdb100genus is NA fill newGenus with rdp100genus
dfbac3$newGenus[dfbac3$NAcount == 0 & dfbac3$mismatching == 0] <-
  dfbac3$Genus[dfbac3$NAcount == 0 & dfbac3$mismatching == 0]
dfbac3$newGenus[dfbac3$NAcount == 1 & dfbac3$mismatching == 0] <-
  dfbac3$Genus[dfbac3$NAcount == 1 & dfbac3$mismatching == 0]
dfbac3$newGenus[is.na(dfbac3$newGenus) & dfbac3$NAcount == 1 & dfbac3$mismatching == 0] <-
  dfbac3$gtdb100genus[is.na(dfbac3$newGenus) & dfbac3$NAcount == 1 & dfbac3$mismatching == 0]
dfbac3$newGenus[is.na(dfbac3$newGenus) & dfbac3$NAcount == 1 & dfbac3$mismatching == 0] <-
  dfbac3$rdp100genus[is.na(dfbac3$newGenus) & dfbac3$NAcount == 1 & dfbac3$mismatching == 0]
#make an empty column named newGenusweak
dfbac3$newGenusweak <- NA
#fill newGenus weak if NAcount == 2, if Genus is NA fill newGenusweak with gtdb100genus if
# gtdb100genus is NA fill newGenusweak with rdp100genus
dfbac3$newGenusweak[dfbac3$NAcount == 2 & dfbac3$mismatching == 0] <-
  dfbac3$Genus[dfbac3$NAcount == 2 & dfbac3$mismatching == 0]
dfbac3$newGenusweak[is.na(dfbac3$newGenusweak) & dfbac3$NAcount == 2 & dfbac3$mismatching == 0] <-
  dfbac3$gtdb100genus[is.na(dfbac3$newGenusweak) & dfbac3$NAcount == 2 & dfbac3$mismatching == 0]
dfbac3$newGenusweak[is.na(dfbac3$newGenusweak) & dfbac3$NAcount == 2 & dfbac3$mismatching == 0] <-
  dfbac3$rdp100genus[is.na(dfbac3$newGenusweak) & dfbac3$NAcount == 2 & dfbac3$mismatching == 0]
#save dfbac for later
saveRDS(dfbac3, "./03-Routputfiles/B_SangerSeq/B3_ListOfDoneIsolatesRevised.rds")

#count number of isolates per newGenus that are not NA
genuscounts <- dfbac3[!is.na(dfbac3$newGenus), ] %>%
  group_by(newGenus) %>%
  summarise(count = n())
genuscounts <- data.frame(genuscounts)
#count number  of isolates per newGenusweak that are not NA
genuscountsweak <- dfbac3[!is.na(dfbac3$newGenusweak), ] %>%
  group_by(newGenusweak) %>%
  summarise(count = n())
genuscountsweak <- data.frame(genuscountsweak)
#merge genuscounts and genuscountsweak
genuscounts <- merge(genuscounts, genuscountsweak,
                     by.x = "newGenus", by.y = "newGenusweak", all = TRUE)
#rename genuscounts$newGenus to genuscounts$th
colnames(genuscounts)[1] <- "th"
#make a dataframe named taxonomy helper that contains the columns Domain, Phylum,
# Class, Order, Family, newGenus, newGenusweak and keep only unique rows
taxonomyhelper <- dfbac3 %>%
  select(Kingdom, Phylum, Class, Order, Family, newGenus, newGenusweak) %>%
  distinct()
#create a column in taxonomyhelper named th that contains newGenus and if newGenus
# is NA it contains newGenusweak
taxonomyhelper$th <- taxonomyhelper$newGenus
taxonomyhelper$th[is.na(taxonomyhelper$th)] <- taxonomyhelper$newGenusweak[is.na(taxonomyhelper$th)]
#join genuscounts and taxonomyhelper
genuscounts <- plyr::join(genuscounts, taxonomyhelper, by = "th", type = "left")
#paste count.x and count.y together to counts
genuscounts$counts <- ifelse(is.na(genuscounts$count.y),
                             genuscounts$count.x,
                             ifelse(is.na(genuscounts$count.x),
                                    paste0("(", genuscounts$count.y, ")"),
                                    paste0(genuscounts$count.x,
                                           " (", genuscounts$count.y, ")")))
#remove duplicates of Genus in genuscounts
genuscounts <- genuscounts[!duplicated(genuscounts$th), ]
#for internal check calculate the sums of genuscounts$count.x and genuscounts$count.y
sum(genuscounts$count.x, na.rm = TRUE)
sum(genuscounts$count.y, na.rm = TRUE)
#make a nice structured table from genuscounts that contains Phylum, Class, Family, th and counts
finaltable <- genuscounts %>%
  select(Phylum, Class, Family, th, counts, count.x, count.y) %>%
  distinct()
#rename finaltable$th to finaltable$Genus
colnames(finaltable)[4] <- "Genus"
#sort by Phylum, Class, count.x (descending) and count.y
finaltable <- finaltable[order(finaltable$Phylum, finaltable$Class, -finaltable$count.x, -finaltable$count.y), ]
#for internal check calculate the sums of finaltable$count.x and finaltable$count.y
sum(finaltable$count.x, na.rm = TRUE)
sum(finaltable$count.y, na.rm = TRUE)
#save finaltable as RDS
saveRDS(finaltable, "./03-Routputfiles/B_SangerSeq/B3_Isolates100CI.rds")
###############################################################################
#Part 2: minboot 90
#create a table with the isolates that had no match on Genus level on minBoot
# 100 with any of the three databases
nomatch100 <- dfbac13
#prep Part2:
#count NA per row and add to nomatch100
nomatch100$NAcount <- rowSums(is.na(nomatch100[, c(10, 13, 17)]))
#make an empty dataframe with colnames "parameter" and "value"
df4 <- data.frame("parameter" = 0, "value" = 0,
                  "nNA_matching" = 0, "nNA_not_matching" = 0)
#filter isolates per sum of NA in silva90, gtdb90genus and rdp90genus
nomatch1000 <- nomatch100[nomatch100$NAcount == 0, ]
nomatch1001 <- nomatch100[nomatch100$NAcount == 1, ]
nomatch1002 <- nomatch100[nomatch100$NAcount == 2, ]
nomatch1003 <- nomatch100[nomatch100$NAcount == 3, ]
###############################################################################
#stats Part2:
#counts of 0 NA isolates
df4[1, ] <- c("number of isolates with 0 NA",
              nrow(nomatch1000),
              nrow(nomatch1000[nomatch1000$silva90 == nomatch1000$gtdb90genus & nomatch1000$silva90 == nomatch1000$rdp90genus, ]),
              nrow(nomatch1000[nomatch1000$silva90 != nomatch1000$gtdb90genus | nomatch1000$silva90 != nomatch1000$rdp90genus, ]))
#counts of 1 NA isolates
df4[2, ] <- c("number of isolates with 1 NA",
              nrow(nomatch1001),
              nrow(nomatch1001[nomatch1001$silva90 == nomatch1001$gtdb90genus & nomatch1001$silva90 == nomatch1001$rdp90genus, ]),
              nrow(nomatch1001) - nrow(nomatch1001[nomatch1001$silva90 == nomatch1001$gtdb90genus & nomatch1001$silva90 == nomatch1001$rdp90genus, ]))
#counts of 2 NA isolates
df4[3, ] <- c("number of isolates with 2 NA",
              nrow(nomatch1002), nrow(nomatch1002), 0)
#counts of 3 NA isolates
df4[4, ] <- c("number of isolates with 3 NA",
              nrow(nomatch1003), nrow(nomatch1003), 0)
print(df4)
###############################################################################
#mismatches Part2:
#show the isolates where the strings in silva90, gtdb90genus and rdp90genus are not the same
mismatch90 <- nomatch100[nomatch100$silva90 != nomatch100$gtdb90genus | nomatch100$silva90 != nomatch100$rdp90genus | nomatch100$gtdb90genus != nomatch100$rdp90genus, ]
#reduce mismatch90 to the columns ids, Family, silva90, gtdb90genus, rdp90genus
mismatch90 <- mismatch90[, c(1, 4, 5, 6, 7, 8, 10, 13, 14, 17, 18)]
#remove isolates where the whole row is NA
mismatch90 <- mismatch90[!is.na(mismatch90$ids), ]
#save mismatch90 a txt
write.table(mismatch90, "./03-Routputfiles/B_SangerSeq/B3_mismatchedGenera90CI.txt",
            sep = ";")
###############################################################################
#make a mismatch overview for Table S3
mismatchred90 <- mismatch90
#paste silva90 gtdb90genus and rdp90genus together separated by "/"
mismatchred90$silva90 <- paste(mismatchred90$silva90, mismatchred90$gtdb90genus,
                               mismatchred90$rdp90genus, sep = "/")
#remove gtdb90genus and rdp90genus and ids
mismatchred90 <- mismatchred90[, c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)]
#count number of isolates per silva90
mismatchred90 <- mismatchred90 %>%
  group_by(silva90) %>%
  summarise(count = n())
mismatchred90 <- data.frame(mismatchred90)
#rename mismatchred90$silva90 to mismatchred90$Genus
colnames(mismatchred90)[1] <- "Genus"
#save mismatchred90 as tsv to 03-Routputfiles/B_SangerSeq
saveRDS(mismatchred90,
        "./03-Routputfiles/B_SangerSeq/B3_mismatchedGenera90CI_bac.rds")
###############################################################################
#Generate structure of Table S2
#prepare a list for the supplements that contains counts per Genus
nomatch100 <- dfbac13
#count NA per row and add to nomatch100
nomatch100$NAcount <- rowSums(is.na(nomatch100[, c(10, 13, 17)]))
#make an empty column named mismatching
nomatch100$mismatching <- 0
#fill mismatching with 1 if id is present in mismatch
nomatch100$mismatching[nomatch100$ids %in% mismatch90$ids] <- 1
#make an empty column named newGenus
nomatch100$newGenus <- NA
#fill newGenus with Genus if NAcount == 0 or 1, if Genus is NA fill newGenus with gtdb100genus if gtdb100genus is NA fill newGenus with rdp100genus
nomatch100$newGenus[nomatch100$NAcount == 0 & nomatch100$mismatching == 0] <- nomatch100$silva90[nomatch100$NAcount == 0 & nomatch100$mismatching == 0]
nomatch100$newGenus[nomatch100$NAcount == 1 & nomatch100$mismatching == 0] <- nomatch100$silva90[nomatch100$NAcount == 1 & nomatch100$mismatching == 0]
nomatch100$newGenus[is.na(nomatch100$newGenus) & nomatch100$NAcount == 1 & nomatch100$mismatching == 0] <- nomatch100$gtdb90genus[is.na(nomatch100$newGenus) & nomatch100$NAcount == 1 & nomatch100$mismatching == 0]
nomatch100$newGenus[is.na(nomatch100$newGenus) & nomatch100$NAcount == 1 & nomatch100$mismatching == 0] <- nomatch100$rdp90genus[is.na(nomatch100$newGenus) & nomatch100$NAcount == 1 & nomatch100$mismatching == 0]
#make an empty column named newGenusweak
nomatch100$newGenusweak <- NA
#fill newGenus weak if NAcount == 2, if Genus is NA fill newGenusweak with gtdb100genus if gtdb100genus is NA fill newGenusweak with rdp100genus
nomatch100$newGenusweak[nomatch100$NAcount == 2 & nomatch100$mismatching == 0] <- nomatch100$silva90[nomatch100$NAcount == 2 & nomatch100$mismatching == 0]
nomatch100$newGenusweak[is.na(nomatch100$newGenusweak) & nomatch100$NAcount == 2 & nomatch100$mismatching == 0] <- nomatch100$gtdb90genus[is.na(nomatch100$newGenusweak) & nomatch100$NAcount == 2 & nomatch100$mismatching == 0]
nomatch100$newGenusweak[is.na(nomatch100$newGenusweak) & nomatch100$NAcount == 2 & nomatch100$mismatching == 0] <- nomatch100$rdp90genus[is.na(nomatch100$newGenusweak) & nomatch100$NAcount == 2 & nomatch100$mismatching == 0]
#count number of isolates per newGenus that are not NA
genuscounts90 <- nomatch100[!is.na(nomatch100$newGenus), ] %>%
  group_by(newGenus) %>%
  summarise(count = n())
genuscounts90 <- data.frame(genuscounts90)
#count number  of isolates per newGenusweak that are not NA
genuscountsweak90 <- nomatch100[!is.na(nomatch100$newGenusweak), ] %>%
  group_by(newGenusweak) %>%
  summarise(count = n())
genuscountsweak90 <- data.frame(genuscountsweak90)
#merge genuscounts and genuscountsweak
genuscounts90 <- merge(genuscounts90, genuscountsweak90,
                       by.x = "newGenus", by.y = "newGenusweak", all = TRUE)
#rename genuscounts$newGenus to genuscounts$th
colnames(genuscounts90)[1] <- "th"
#make a dataframe named taxonomy helper that contains the columns Domain, Phylum,
# Class, Order, Family, newGenus, newGenusweak and keep only unique rows
taxonomyhelper90 <- nomatch100 %>%
  select(Kingdom, Phylum, Class, Order, Family, newGenus, newGenusweak) %>%
  distinct()
#create a column in taxonomyhelper named th that contains newGenus and if
# newGenus is NA it contains newGenusweak
taxonomyhelper90$th <- taxonomyhelper90$newGenus
taxonomyhelper90$th[is.na(taxonomyhelper90$th)] <- taxonomyhelper90$newGenusweak[is.na(taxonomyhelper90$th)]
#join genuscounts and taxonomyhelper
genuscounts90 <- plyr::join(genuscounts90, taxonomyhelper90, by = "th", type = "left")
#paste count.x and count.y together to counts
genuscounts90$counts <- ifelse(is.na(genuscounts90$count.y),
                               genuscounts90$count.x,
                               ifelse(is.na(genuscounts90$count.x),
                                      paste0("(", genuscounts90$count.y, ")"),
                                      paste0(genuscounts90$count.x,
                                             " (", genuscounts90$count.y, ")")))
#remove duplicates of Genus in genuscounts
genuscounts90 <- genuscounts90[!duplicated(genuscounts90$th), ]
#for internal check calculate the sums of genuscounts$count.x and genuscounts$count.y
sum(genuscounts90$count.x, na.rm = TRUE)
sum(genuscounts90$count.y, na.rm = TRUE)
#make a nice structured table from genuscounts that contains Phylum, Class, Family, th and counts
finaltable90 <- genuscounts90 %>%
  select(Phylum, Class, Family, th, counts, count.x, count.y) %>%
  distinct()
#rename finaltable$th to finaltable$Genus
colnames(finaltable90)[4] <- "Genus"
#sort by Phylum, Class, count.x (descending) and count.y
finaltable90 <- finaltable90[order(finaltable90$Phylum, finaltable90$Class, -finaltable90$count.x, -finaltable90$count.y), ]
#for internal check calculate the sums of finaltable$count.x and finaltable$count.y
sum(finaltable90$count.x, na.rm = TRUE)
sum(finaltable90$count.y, na.rm = TRUE)
#save finaltable as RDS
saveRDS(finaltable90, "./03-Routputfiles/B_SangerSeq/B3_Isolates90CI.rds")
###############################################################################
#make structure of Table S3
#prepare mismatchred and mismatchred90 for the manuscript
#add column to mismatchred called separator and fill it with mismatchre
mismatchred$separator <- "mismatchred"
#add column to mismatchred90 called separator and fill it with mismatchred90
mismatchred90$separator <- "mismatchred90"
#combine mismatchred and mismatchred90
mismatchredall <- rbind(mismatchred, mismatchred90)
#save mismatchredall as RDS
saveRDS(mismatchredall, "./03-Routputfiles/B_SangerSeq/B3_mismatchedGenera.rds")
###############################################################################
# manually checked mismatches (using blast and check for new nomenclature)
## Hafnia-Obesumbacterium vs Hafnia
# accepted as matching -> Hafnia-Obesumbacterium is changed to Hafnia in column
# Genus
## Myroides vs Flavobacterium
# check sequence with blast -> Myroides is correct probably correct
## Enemella vs Ponticoccus
#Ponticoccus is correct -> it is a heterotypic synonym of Enemella, but they
# have differen lineage descriptions in the blast db
## Mammaliicoccus vitulinus vs Staphylococcus vitulinus
# Mammaliicoccus viutlinus is the more recent name -> change Staphylococcus
# vitulinus to Mammaliicoccus vitulinus
###############################################################################


####statistics for fungi isolates#############################
#General statsitics for fungi isolates, how many isolates are there, how many
# have a genus classification in silva100, how many have a genus classification
# in silva90
#undock dffun for new analysis
dffun1 <- dffun
#make an empty dataframe with colnames "parameter" and "value"
sdffuntax <- data.frame("parameter" = 0, "value" = 0)
#count number of isolates
sdffuntax[1, ] <- c("number of isolates", nrow(dffun1))
#count number of isolates per genus that are not NA
sdffuntax[2, ] <- c("number of isolates in silva100 that are not NA",
                    nrow(dffun1[!is.na(dffun1$Genus), ]))
#count number of isolates per silva90 that are not NA
sdffuntax[3, ] <- c("number of isolates in silva90 that are not NA",
                    nrow(dffun1[!is.na(dffun1$silva90), ]))
print(sdffuntax)
####EXKURS########################################################
#statistic for dffun with sequences with seqlength < 500
dffun2 <- dffun1[dffun1$seqlength < 500, ]
#make an empty dataframe with colnames "parameter" and "value"
sdffuntaxshort <- data.frame("parameter" = 0, "value" = 0)
#count number of isolates
sdffuntaxshort[1, ] <- c("number of isolates with seqlength < 500",
                         nrow(dffun2))
#count number of isolates per genus that are not NA
sdffuntaxshort[2, ] <- c("number of isolates in silva100 that are not NA",
                         nrow(dffun2[!is.na(dffun2$Genus), ]))
#count number of isolates per silva90 that are not NA
sdffuntaxshort[3, ] <- c("number of isolates in silva90 that are not NA",
                         nrow(dffun2[!is.na(dffun2$silva90), ]))
print(sdffuntaxshort)
##############################################################
#produce tables for supplements, containing classified taxa
#statistic for dfbac to check the classification metrics
#make dffun2 from dffun1
dffun2 <- dffun1
#keep only the columns Kingdom to silva90
dffun2 <- dffun2[, c(1, 4:10)]
#make dffun2100 from dffun2 that only contains isolates with an entry in Genus
dffun2100 <- dffun2[!is.na(dffun2$Genus), ]
#make dffun290 from dffun2 that only countains isolates that are NA in Genus
# but have an entry in silva90
dffun290 <- dffun2[is.na(dffun2$Genus) & !is.na(dffun2$silva90), ]
#make dffun2lowres from dffun2 that only contains isolates that are NA in
# Genus and silva90
dffun2lowres <- dffun2[is.na(dffun2$Genus) & is.na(dffun2$silva90), ]
#summarize the statistics of dffun2100, dffun290 and dffun2lowres
#make an empty dataframe with colnames "parameter" and "value"
df5 <- data.frame("parameter" = 0, "value" = 0)
#count number of isolates with an entry in Genus
df5[1, ] <- c("number of isolates with an entry in Genus", nrow(dffun2100))
#count number of isolates with an entry in silva90
df5[2, ] <- c("number of isolates with an entry in silva90", nrow(dffun290))
#count number of isolates without entry in Genus or silva90
df5[3, ] <- c("number of isolates without Genus classification",
              nrow(dffun2lowres))
print(df5)
###############################################################################
#create structure for Table S4
#count number of isolates per Genus in dffun2100
dffun2100counts <- dffun2100 %>%
  group_by(Genus) %>%
  summarise(count = n())
dffun2100counts <- data.frame(dffun2100counts)
#add taxonomic levels to dffun2100counts
dffun2100counts <- plyr::join(dffun2100counts, dffun2100[, 3:7],
                              by = "Genus", type = "left")
#reorder columns in dffun2100counts to Kingdom, Phylum, Class, Order, Family,
# Genus, count
dffun2100counts <- dffun2100counts[, c(3, 4, 5, 6, 1, 2)]
#count number of isolates per Genus in dffun290
dffun290counts <- dffun290 %>%
  group_by(silva90) %>%
  summarise(count = n())
dffun290counts <- data.frame(dffun290counts)
#add taxonomic levels to dffun290counts
dffun290counts <- plyr::join(dffun290counts, dffun290[, c(3:6, 8)],
                             by = "silva90", type = "left")
#rename count to weakcount
colnames(dffun290counts)[2] <- "weakcount"
#rename silva90 to Genus
colnames(dffun290counts)[1] <- "Genus"
#join dffun2100counts and dffun290counts
dffun2100counts <- plyr::join(dffun2100counts, dffun290counts,
                              by = "Genus", type = "full")
#remove duplicates
dffun2100counts <- dffun2100counts[!duplicated(dffun2100counts$Genus), ]
#order by Phylum, Class, Order, Family, count (descending), weakcount (descending)
dffun2100counts <- dffun2100counts[order(dffun2100counts$Phylum,
                                         dffun2100counts$Class,
                                         dffun2100counts$Order,
                                         dffun2100counts$Family,
                                         -dffun2100counts$count,
                                         -dffun2100counts$weakcount), ]

#remove NA rows from dffun2lowres
dffun2lowres <- dffun2lowres[!is.na(dffun2lowres$Kingdom), ]
#remove Genus
dffun2lowres <- dffun2lowres[, 1:7]
#find the column with the last NA entry counting from left
dffun2lowres$firstNA <- apply(dffun2lowres, 1, function(x) min(which(is.na(x))))
#give back a list of the columns that are indext in dffun2lowres$lastNA
dffun2lowres$lastNAcol <- apply(dffun2lowres, 1,
                                function(x) names(x)[min(which(is.na(x)) - 1)])
#count number of isolates per lastNAcol
dffun2lowrescounts <- dffun2lowres %>%
  group_by(lastNAcol) %>%
  summarise(count = n())

print(dffun2lowrescounts)

#add column to dffun2lowres named newtax that contains "unclassified " and the
# content of the column that is indexed in lastNAcol
dffun2lowres$newtax <- paste0("unclassified ",
                              apply(dffun2lowres, 1,
                                    function(x) x[min(which(is.na(x)) - 1)]))

#count number of isolates per newtax
dffun2lowrescounts2 <- dffun2lowres %>%
  group_by(newtax) %>%
  summarise(count = n())
dffun2lowrescounts2 <- data.frame(dffun2lowrescounts2)

#join dffun2lowrescounts2 and dffun2lowres by new tax
dffun2lowrescounts2 <- plyr::join(dffun2lowrescounts2, dffun2lowres[, 8:10],
                                  by = "newtax", type = "left")
#remove duplicates
dffun2lowrescounts2 <- dffun2lowrescounts2[!duplicated(dffun2lowrescounts2$newtax), ]
#sort by firstNA (descending) and count (descending)
dffun2lowrescounts2 <- dffun2lowrescounts2[order(-dffun2lowrescounts2$firstNA,
                                                 -dffun2lowrescounts2$count), ]


#remove rows with Kingdom in lastNAcol
dffun2lowrescounts2red <- dffun2lowrescounts2[!dffun2lowrescounts2$lastNAcol == "Kingdom", ]
#remove column firstNA
dffun2lowrescounts2red <- dffun2lowrescounts2red[, c(1, 2, 4)]

print(dffun2lowrescounts2red)


#####additional way to add the unclassified taxa to the table of classified taxa
#undock dffun2lowres for new analysis
funlowres <- dffun2lowres
#put newtax from funlowres to the column that is indexed in firstNA
for (i in 1:nrow(funlowres)){
  funlowres[i, funlowres$firstNA[i] - 1] <- funlowres$newtax[i]
}
#count number of isolates per newtax
funlowrescounts <- funlowres %>%
  group_by(newtax) %>%
  summarise(count = n())

#join funlowrescounts and funlowres by newtax
funlowres <- plyr::join(funlowres, funlowrescounts,
                        by = "newtax", type = "left")
#remove duplicates
funlowres <- funlowres[!duplicated(funlowres$newtax), ]

#keep only Phylum, Class, Order, Family, Genus, count from funlowres
funlowres <- funlowres[, c(3:7, 11)]
#add weakcount to funlowres with the values from count
funlowres$weakcount <- funlowres$count
#fill count with NA
funlowres$count <- NA

#combine dfun2100counts and funlowres
dffunall <- rbind(dffun2100counts, funlowres)
#order by Phylum, Class, Order, Family, count (descending), weakcount (descending)
dffunall <- dffunall[order(dffunall$Phylum,
                           dffunall$Class,
                           dffunall$Order,
                           dffunall$Family,
                           -dffunall$count,
                           -dffunall$weakcount), ]

#remove rows with NA in Phylum
dffunall <- dffunall[!is.na(dffunall$Phylum), ]
#rename count with counts 100CI
colnames(dffunall)[6] <- "counts 100% CI"
#rename weakcount with counts 90CI
colnames(dffunall)[7] <- "counts 90% CI"

#replace all "_" with " " in dffunall
dffunall <- data.frame(apply(dffunall, 2, function(x) gsub("_", " ", x)))

#save dffunall as RDS
saveRDS(dffunall, "./03-Routputfiles/B_SangerSeq/B3_IsolatesFungi.rds")
###############################################################################

###analysis on species level
#make dfbacspec from dfbac
dfbacspec <- dfbac
#keep only the columns ids, Kingdom, Phylum, Class, Order, Family, gtdb100genus,
# gtdb100species,  rdp100genus, rdp100species
dfbacspec <- dfbacspec[, c(1, 4:8, 11, 12, 15, 16)]
#keep only rows that have an entry in gtdb100species or rdp100species
dfbacspec <- dfbacspec[!is.na(dfbacspec$gtdb100species) | !is.na(dfbacspec$rdp100species), ]

#create vector with paste of gtdb100species and rdp100species
dfbacspec$species <- paste(dfbacspec$gtdb100species,
                           dfbacspec$rdp100species,
                           sep = "/")
#replace "_" with " " in dfbacspec$species
dfbacspec$species <- gsub("_", " ", dfbacspec$species)
#remove "NA/" and "/NA" from dfbacspec$species
dfbacspec$species <- gsub("NA/", "", dfbacspec$species)
dfbacspec$species <- gsub("/NA", "", dfbacspec$species)

#count number of isolates per species
dfbacspeccounts <- dfbacspec %>%
  group_by(species) %>%
  summarise(count = n())
#make dfbacspeccounts a dataframe
dfbacspeccounts <- data.frame(dfbacspeccounts)
#order by count (descending)
dfbacspeccounts <- dfbacspeccounts[order(-dfbacspeccounts$count), ]
#save as RDS
saveRDS(dfbacspeccounts, "./03-Routputfiles/B_SangerSeq/B3_IsolatesSpecies.rds")
###############################################################################
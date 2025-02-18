################################################################
################################################################
###
### RUN ON SERVER !!!!!!!!!!!!!!
###
################################################################
################################################################
# D2v1_Isolates_Figure_P2.r creates a tree uses 
## 1. a muscle alignment 
## 2. Gblocks to remove gaps
## 3. the DECIPHER::TreeLine for tree generation
# 
## Problem with this approach is that the alpha bacteria cluster 
#  closer to the Bacillota than to the pseudomonata
################################################################
################################################################
### FIGURE 1 ###
## Figure 1 shows the following:
## 1. A tree with one representative of each Genus based on the 
#     16S/18S rRNA gene sequences, representative is isolate with 
#     the longest sequence in the most prevalent Cluster within a 
#     Genus
## 2. Additional information about the isolates and their 
#     diversity within a Genus
##   2a. column for each batch with points
##   2b. points are sized by the grouped number of 16S/18S rRNA 
#        gene genotypes
##   2c. points are coloured by the last slice the isolate was 
#        found in only sequences >1000bp with taxonomy >=2 
#        matching databases were used
###############################################################
###############################################################
### THIS IS PART 2 - TREE GENERATION ###
###############################################################
library(ggplot2)
library(DECIPHER)
library(dplyr)
library(Biostrings)
library(stringr)
library(ape)
library(phangorn)
library(msa)

################################################################
#### FIGURE 1A - BACTERIA ####
### Load data ###
ssbrep<-readRDS("./03-Routputfiles/C_IsolateTree/C1_RepresentativeIsolates_bac.RDS")
## use newGenus from metadata to name the sequences in ssbrep
names(ssbrep)<-metadata(ssbrep)$newGenus

################################################################
# alignment with muscle
alignment <- msa(ssbrep, method = "Muscle")
# write alignment as DNAStringSet
alstring <- as(alignment, "DNAStringSet")
# save aligned sequences as fasta
writeXStringSet(alstring, file="./03-Routputfiles/C_IsolateTree/C2_muscleAlignedRepresentatives.fasta")
# save aligned sequences as rds
saveRDS(alstring, file="./03-Routputfiles/C_IsolateTree/C2_muscleAlignedRepresentatives.rds")

################################################################
################################################################
###
### NON R CODE ####
###
# run Gblocks from a conda environment called phyml
system("/home/roch/miniconda3/envs/phyml/bin/Gblocks ./03-Routputfiles/C_IsolateTree/C2_muscleAlignedRepresentatives.fasta -t=d -b5=n") #-b5=n removes gaps

################################################################
################################################################
# read in Gblocks output and make it a DNAStringSet
gblock<-readDNAStringSet("./03-Routputfiles/C_IsolateTree/C2_muscleAlignedRepresentatives.fasta-gb", format="fasta")

################################################################
# run TreeLine on the Gblocked alignment
set.seed(123)
bactree<-TreeLine(gblock, maxTime = 1, processors = 100)
# save bactree as RDS and txt
saveRDS(bactree,"./03-Routputfiles/C_IsolateTree/C2_bactree.RDS")
WriteDendrogram(bactree,file="./03-Routputfiles/C_IsolateTree/C2_bactree.txt")
# get tree attributes
attributes(bactree)


################################################################
#### FIGURE 1B - FUNGI ####
### Load data ###
ssfrep<-readRDS("./03-Routputfiles/C_IsolateTree/C1_RepresentativeIsolates_fun.RDS")
## use newtax from metadata to name the sequences in ssfrep
names(ssfrep)<-metadata(ssfrep)$newtax

################################################################
# alignment with muscle
alignment <- msa(ssfrep, method = "Muscle")
# write alignment as DNAStringSet
alstring <- as(alignment, "DNAStringSet")
# save aligned sequences as fasta
writeXStringSet(alstring, file="./03-Routputfiles/C_IsolateTree/C2_muscleAlignedRepresentatives_fun.fasta")
# save aligned sequences as rds
saveRDS(alstring, file="./03-Routputfiles/C_IsolateTree/C2_muscleAlignedRepresentatives_fun.rds")

################################################################
################################################################
###
### NON R CODE ####
###
# run Gblocks from a conda environment called phyml
system("/home/roch/miniconda3/envs/phyml/bin/Gblocks ./03-Routputfiles/C_IsolateTree/C2_muscleAlignedRepresentatives_fun.fasta -t=d -b5=n") #-b5=n removes gaps

################################################################
################################################################
# read in Gblocks output and make it a DNAStringSet
gblock<-readDNAStringSet("./03-Routputfiles/C_IsolateTree/C2_muscleAlignedRepresentatives_fun.fasta-gb", format="fasta")

################################################################
# run TreeLine on the Gblocked alignment
set.seed(123)
funtree<-TreeLine(gblock, maxTime = 1, processors = 100)
# save funtree as RDS and txt
saveRDS(funtree,"./03-Routputfiles/C_IsolateTree/C2_funtree.RDS")
WriteDendrogram(funtree,file="./03-Routputfiles/C_IsolateTree/C2_funtree.txt")
# get tree attributes
attributes(funtree)


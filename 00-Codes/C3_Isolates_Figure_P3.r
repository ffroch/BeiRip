############################################################################################################
### FIGURE 1 ###
#Figure 1 shows the following:
#1. A tree with one representative of each Genus based on the 16S/18S rRNA gene sequences, representative is isolate with the longest sequence in the most prevalent Cluster within a Genus
#2. Additional information about the isolates and their diversity within a Genus
##2a. column for each batch with points
##2b. points are sized by the grouped number of 16S/18S rRNA gene genotypes
##2c. points are coloured by the last slice the isolate was found in
#only sequences >1000bp with taxonomy >=2 matching databases were used
############################################################################################################
### THIS IS PART 3 - CREATE METADATA FOR 2. ###
############################################################################################################
library(ggplot2)
library(DECIPHER)
library(dplyr)
library(forcats)
library(ggtree)
library(phylobase)
library(stringr)
library(phylotools)
library(ape)
library(RColorBrewer)
library(viridis)
library(cowplot)
############################################################################################################
#### FIGURE 1A - BACTERIA ####
### Load data ###
data<-readRDS("./03-Routputfiles/C_IsolateTree/C1_IsolatesUsedForFig1bac.RDS")
############################################################################################################
### Get batch and slice information ###
#split ids to get batch and slice for the isolates
data$batch<-sapply(strsplit(as.character(data$ids), "_"), "[", 1)
data$slice<-sapply(strsplit(as.character(data$ids), "_"), "[", 2)
#remove S in data$slice and make it numeric
data$slice<-as.numeric(gsub("S", "", data$slice))
#reduce data to relevant columns
datared<-data%>%
          select(ids, seqlength, Phylum, Class, Order, Family, newGenus, cluster, count, diversity, batch, slice)
#paste the batch and the Genus together
datared$Genusbatch<-paste(datared$batch, datared$newGenus, sep="_")
#keep only the isolate with the highest value in slice for each unique data$Genusbatch
dataredlate<-datared[order(datared$slice, decreasing=TRUE),]
dataredlate<-dataredlate[!duplicated(dataredlate$Genusbatch),]
# keep only the isolate with the lowest value in slice for each unique data$Genusbatch
dataredearly<-datared[order(datared$slice, decreasing=FALSE),]
dataredearly<-dataredearly[!duplicated(dataredearly$Genusbatch),]

### layout transformation ###
#expand the object datared to have all Genus for each batch
dataex<-expand.grid(batch=unique(datared$batch), newGenus=unique(datared$newGenus))
dataex$Genusbatch<-paste(dataex$batch, dataex$newGenus, sep="_")
#remove duplicates in the Genus from datared[,1:6]
taxo<-datared[,3:7]
taxo<-taxo[!duplicated(taxo$newGenus),]
# join datex and taxo
dataex<-plyr::join(dataex, taxo, by="newGenus", type="left")
# add dataredlate and dataredearly
dataex<-plyr::join(dataex, dataredlate[,c(10,12,13)], by="Genusbatch", type="left")
dataex<-plyr::join(dataex, dataredearly[,c(12,13)], by="Genusbatch", type="left")
# rename slice columns to latesttp and earliesttp
colnames(dataex)[9]<-"latesttp"
colnames(dataex)[10]<-"earliesttp"
# add presence column
dataex$presence<-ifelse(is.na(dataex$latesttp), 0, 1)
#add 0.5 to the diversity of isolates with diversity=0
dataex$diversity[is.na(dataex$diversity)]<-0.5
############################################################################################################
#save dataex as RDS
saveRDS(dataex, "./03-Routputfiles/C_IsolateTree/C3_bacmeta.RDS")


############################################################################################################
#### FIGURE 1A - FUNGI ####
### Load data ###
dataf<-readRDS("./03-Routputfiles/C_IsolateTree/C1_IsolatesUsedForFig1fun.RDS")
#rename first column to ids
colnames(dataf)[1]<-"ids"
############################################################################################################
### Get batch and slice information ###
#split ids to get batch and slice for the isolates
dataf$batch<-sapply(strsplit(as.character(dataf$ids), "_"), "[", 1)
dataf$slice<-sapply(strsplit(as.character(dataf$ids), "_"), "[", 2)
#remove S in data$slice and make it numeric
dataf$slice<-as.numeric(gsub("S", "", dataf$slice))
#reduce data to relevant columns
datafred<-dataf%>%
          select(ids, lengthfun, Phylum, Class, Order, Family, newtax, cluster, count, diversity, batch, slice)
#paste the batch and the Genus together
datafred$taxbatch<-paste(datafred$batch, datafred$newtax, sep="_")
#keep only the isolate with the highest value in slice for each unique data$Genusbatch
datafredlate<-datafred[order(datafred$slice, decreasing=TRUE),]
datafredlate<-datafredlate[!duplicated(datafredlate$taxbatch),]

#keep only the isolate with the lowest value in slice for each unique data$Genusbatch
datafredearly<-datafred[order(datafred$slice, decreasing=FALSE),]
datafredearly<-datafredearly[!duplicated(datafredearly$taxbatch),]

### layout transformation ###
# expand the object dataex to have all Genus for each batch
datafex<-expand.grid(batch=unique(datafred$batch), newtax=unique(datafred$newtax))
datafex$taxbatch<-paste(datafex$batch, datafex$newtax, sep="_")
#remove duplicates in the Genus from dataredlate
taxof<-datafred[,3:7]
taxof<-taxof[!duplicated(taxof$newtax),]
# join datafex and taxof
datafex<-plyr::join(datafex, taxof, by="newtax", type="left")
# add datafredlate and datafredearly
datafex<-plyr::join(datafex, datafredlate[,c(10,12,13)], by="taxbatch", type="left")
datafex<-plyr::join(datafex, datafredearly[,c(12,13)], by="taxbatch", type="left")
# rename slice columns to latesttp and earliesttp
colnames(datafex)[9]<-"latesttp"
colnames(datafex)[10]<-"earliesttp"
# add presence column
datafex$presence<-ifelse(is.na(datafex$latesttp), 0, 1)
#add 0.5 to the diversity of isolates with diversity=0
datafex$diversity[is.na(datafex$diversity)]<-0.5
############################################################################################################
#save dataex as RDS
saveRDS(datafex, "./03-Routputfiles/C_IsolateTree/C3_funmeta.RDS")
############################################################################################################

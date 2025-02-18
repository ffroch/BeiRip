setwd("/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip")

#RUN ON SERVER!!!
library(Biostrings)
library(DECIPHER)
library(dada2)

#load DNA string set of trimmed isolate sequences
stringset <- readDNAStringSet("./03-Routputfiles/B_SangerSeq/B1_seqstrimclean_35_15_100bp_bac.fa","fasta")
length(stringset)


#SILVA-genus
#ref_fasta<-"C:/Users/rochf/Desktop/veggiemeat/veggiemeat/01-metadata/silva_nr99_v138.1_train_set.fa.gz"
#ref_fasta <- "/data/Unit_LMM/selberherr-group/roch/00_databases/dada2_pretrained_classifier/silva_nr99_v138.1_train_set.fa.gz"
ref_fasta <- "D:/0000Databases/dada2_pretrained_classifier/silva_nr99_v138.2_toGenus_trainset.fa.gz" #https://zenodo.org/records/14169026
#assign Taxonomy
taxtab1 <- assignTaxonomy(stringset, refFasta = ref_fasta,
                          multithread = 10, minBoot = 100)
#write taxonomy table
write.table(taxtab1,
            "./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_silva_genus_CI100.txt",
            sep = ";")
#assign Taxonomy
taxtab1a <- assignTaxonomy(stringset, refFasta = ref_fasta,
                           multithread = 10, minBoot = 90)
#write taxonomy table
write.table(taxtab1a,
            "./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_silva_genus_CI90.txt",
            sep = ";")#change here filename


#GTDB-Species
#ref_fasta<-"C:/Users/rochf/Desktop/Pretrained Classifier/GTDB_bac120_arc53_ssu_r207_fullTaxo.fa.gz"
#ref_fasta <- "/data/Unit_LMM/selberherr-group/roch/00_databases/dada2_pretrained_classifier/GTDB_bac120_arc53_ssu_r207_fullTaxo.fa.gz" #from zenodo https://zenodo.org/record/6655692
ref_fasta <- "D:/0000Databases/dada2_pretrained_classifier/GTDB_bac120_arc53_ssu_r220_fullTaxo.fa.gz" #https://zenodo.org/records/10403693
#assign Taxonomy
taxtab2 <- assignTaxonomy(stringset, refFasta = ref_fasta,
                          multithread = 10, minBoot = 100)
#write taxonomy table
write.table(taxtab2,
            "./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_GTDB_species_CI100.txt",
            sep = ";")
#assign Taxonomy
taxtab2a <- assignTaxonomy(stringset, refFasta = ref_fasta,
                           multithread = 10, minBoot = 90)
#write taxonomy table
write.table(taxtab2a,
            "./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_GTDB_species_CI90.txt",
            sep = ";")
###############################################################################
# create classification for database
taxtab2b <- assignTaxonomy(stringset, refFasta = ref_fasta,
                           multithread = 10, outputBootstraps = TRUE)
nn <- names(stringset)
seq <- rownames(data.frame(taxtab2b))
taxtab2b <- data.frame(IsolateID = nn, Sequence = seq, data.frame(taxtab2b))
write.table(taxtab2b,
            "./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_GTDB_species_inklCI.txt",
            sep = ";", row.names = FALSE)

#RDP-Species
#ref_fasta<-"C:/Users/rochf/Desktop/Pretrained Classifier/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz"
#ref_fasta <- "/data/Unit_LMM/selberherr-group/roch/00_databases/dada2_pretrained_classifier/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz"
ref_fasta <- "D:/0000Databases/dada2_pretrained_classifier/RefSeq_16S_6-11-20_RDPv16_fullTaxo.fa.gz" #https://zenodo.org/records/10403693
#assign Taxonomy
taxtab3<-assignTaxonomy(stringset, refFasta = ref_fasta,
                        multithread = 10, minBoot = 100)
#write taxonomy table
write.table(taxtab3,
            "./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_RDP_species_CI100.txt",
            sep = ";")#change here filename
#assign Taxonomy
taxtab3a <- assignTaxonomy(stringset, refFasta = ref_fasta,
                           multithread = 10, minBoot = 90)
#write taxonomy table
write.table(taxtab3a,
            "./03-Routputfiles/B_SangerSeq/B1_taxonomy_35_15_100_bac_RDP_species_CI90.txt",
            sep = ";")#change here filename

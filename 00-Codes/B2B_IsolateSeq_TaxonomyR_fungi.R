#RUN ON SERVER!!!
#load DNA string set of trimmed isolate sequences
stringset <- readDNAStringSet("./03-Routputfiles/B_SangerSeq/B2_seqstrimclean_35_15_100bp_fun.fa","fasta")
length(stringset)

#SILVA-genus
#load reference data
#ref_fasta <- "C:/Users/rochf/Desktop/BeiRip/02-rawdata/silva_nr_v138_train_set.fa.gz"
#ref_fasta <- "/data/Unit_LMM/selberherr-group/roch/00_databases/dada2_pretrained_classifier/silva_nr_v138_train_set.fa.gz" #https://zenodo.org/records/3731176 v.1
ref_fasta <- "D:/0000Databases/dada2_pretrained_classifier/silva_nr_v138_2_train_set.fa.gz"
#assign Taxonomy
taxtab <- assignTaxonomy(stringset,refFasta = ref_fasta,
                         multithread = 10, minBoot = 100)
#write taxonomy table
write.table(taxtab,
            "./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_silva_genus_CI100.txt",
            sep = ";")

#assign Taxonomy
taxtaba <- assignTaxonomy(stringset, refFasta = ref_fasta,
                          multithread = 10, minBoot = 90)
#write taxonomy table
write.table(taxtaba,
            "./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_silva_genus_CI90.txt",
            sep = ";")


# create classification for database
taxtab2b <- assignTaxonomy(stringset, refFasta = ref_fasta,
                           multithread = 10, outputBootstraps = TRUE)
nn <- names(stringset)
# remove ".EK528" from nn
nn <- sub(".EK528", "", nn)
seq <- rownames(data.frame(taxtab2b))
taxtab2c <- data.frame(IsolateID = nn, Sequence = seq,
                       data.frame(taxtab2b)[, 1:6],
                       tax.Species = NA,
                       data.frame(taxtab2b)[, 7:12],
                       boot.Species = NA,
                       Method = "18SSeq", Database = "PR2")
write.table(taxtab2c,
            "./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_silva_inklCI.txt",
            sep = ";", row.names = FALSE)


# comparison with PR2
ref_fasta <- "D:/0000Databases/dada2_pretrained_classifier/pr2_version_5.0.0_SSU_dada2.fasta.gz" #PR2 database v5.0.0
#assign Taxonomy
taxtab <- assignTaxonomy(stringset,refFasta = ref_fasta,
                         multithread = 10,minBoot = 100)
#write taxonomy table
write.table(taxtab,
            "./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_pr2_species_CI100.txt",
            sep = ";")

#assign Taxonomy
taxtaba <- assignTaxonomy(stringset, refFasta = ref_fasta,
                          multithread = 10, minBoot = 90)
#write taxonomy table
write.table(taxtaba,
            "./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_pr2_species_CI90.txt",
            sep = ";")


# create classification for database
taxtab2b <- assignTaxonomy(stringset, refFasta = ref_fasta,
                           multithread = 10, outputBootstraps = TRUE)
nn <- names(stringset)
# remove ".EK528" from nn
nn <- sub(".EK528", "", nn)
seq <- rownames(data.frame(taxtab2b))
taxtab2c <- data.frame(IsolateID = nn, Sequence = seq,
                       data.frame(taxtab2b)[, 1:6],
                       tax.Species = NA,
                       data.frame(taxtab2b)[, 7:12],
                       boot.Species = NA,
                       Method = "18SSeq", Database = "PR2")
write.table(taxtab2c,
            "./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_pr2_species_inklCI.txt",
            sep = ";", row.names = FALSE)

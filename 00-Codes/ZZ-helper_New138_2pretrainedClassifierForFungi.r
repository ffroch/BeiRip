# code modified from https://zenodo.org/records/3731176 v.1
library(here)
library(dada2)
packageVersion("dada2")


#Create the location where we will download the precursor files,
dir.create(here("precursors"))


#Download and extract the Mothur-formatted taxonomy files (see https://mothur.org/wiki/silva_reference_files/),
download.file(
  "https://mothur.s3.us-east-2.amazonaws.com/wiki/silva.nr_v138_2.tgz",
  here("precursors", "silva.nr_v138_2.tgz")
)

command <- paste("tar -xvf", 
  here("precursors", "silva.nr_v138_2.tgz"),
  "-C", 
  here("precursors")
  )
system(command)
list.files(here("precursors"))


#Download the file `SILVA_138_2_SSURef_tax_silva.fasta.gz` from the Silva website,
download.file(
  "https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_tax_silva.fasta.gz",
  here("precursors", "SILVA_138.2_SSURef_tax_silva.fasta.gz")
)


#Create the DADA2-formatted taxonomy database file,
dada2:::makeTaxonomyFasta_Silva(
  here("precursors", "silva.nr_v138_2.align"), 
  here("precursors", "silva.nr_v138_2.tax"), 
  here("silva_nr_v138_2_train_set.fa.gz")
)


#Create the DADA2-formatted species database file,
dada2:::makeSpeciesFasta_Silva(
  here("precursors", "SILVA_138.2_SSURef_tax_silva.fasta.gz"),
  here("silva_species_assignment_v138_2.fa.gz")
)


## Check new database assignments on test set own data compared to version 138

library(phyloseq)
library(dplyr)

output138 <- read.table("./03-Routputfiles/B_SangerSeq/oldclassifications/B2_taxonomy_35_15_100_fun_silva_genus_CI100.txt", header = TRUE, sep = "\t")
output138_2 <- read.table("./03-Routputfiles/B_SangerSeq/B2_taxonomy_35_15_100_fun_silva_genus_CI100.txt", header = TRUE, sep = "\t")

output138 <- output138 %>% select(seq, genus_old = Genus)
output138_2 <- output138_2 %>% select(seq, genus_new = Genus)

# check if order is identical
identical(output138$seq, output138_2$seq)

# Combine the two data frames
df_combined <- cbind(output138, output138_2)
colnames(df_combined) <- c("seq", "genus_old", "seq2", "genus_new")


# Separate into matching and non-matching genera
matching_genera <- df_combined %>%
  filter(genus_old == genus_new)

non_matching_genera <- df_combined %>%
  filter(genus_old != genus_new)

# Output the two data frames
nrow(matching_genera)
nrow(non_matching_genera)

non_matching_genera[, c("genus_old", "genus_new")]

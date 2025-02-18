########################################################################################################
prodir <- getwd()
# set working directory
setwd("./05-qiime/HTS_18S/")
########################################################################################################
# load removed features
removed <- read.table("bostaurushitstoremovefinal.txt", header = TRUE, sep = "\t")
colnames(removed)[1] <- "Feature.ID"
# load taxonomy
tax <- read.table("./06-taxonomy/01-exported-taxonomy-raw/taxonomy.tsv", header = TRUE, sep = "\t")
# combine removed features and taxonomy
comb <- plyr::join(removed, tax, by = "Feature.ID", type = "left")
# save comb as txt
write.table(comb, file = "./04-asvtable/02-mamfilt/blastcheck.txt", sep = "\t", row.names = FALSE)
########################################################################################################
setwd(prodir)

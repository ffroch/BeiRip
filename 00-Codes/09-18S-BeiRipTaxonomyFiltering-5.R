###################################################################
prodir <- getwd()
# set working directory
setwd("./05-qiime/HTS_18S/")
###################################################################
# load blast results
blast <- read.table("./08-blastdb/comparison.mamfilt.blastout",
                    sep = "\t", header = FALSE)
colnames(blast) <- c("qseqid", "sseqid", "pident", "qcovs",
                     "length", "mismatch", "gapopen", "qstart",
                     "qend", "sstart", "send", "evalue", "bitscore")
# keep only most similar hit per feature
# sort blast by pident
blast <- blast[order(blast$pident, decreasing = TRUE), ]
# remove duplicates
blastred <- blast[!duplicated(blast$sseqid), ]
# load taxonomy
tax <- read.table("./06-taxonomy/02-taxonomymamfilt/taxonomy.tsv",
                  header = TRUE, sep = "\t")
colnames(tax) <- c("sseqid", "Taxon", "Confidence")
# combine removed features and taxonomy
taxblast <- plyr::join(tax, blastred, by = "sseqid")
taxblast <- taxblast[, c(1:5)]
# Keep only rows that does not contain ";" in the Taxon column
taxblastred <- taxblast[!grepl(";", taxblast$Taxon), ]
# save taxblast as txt
write.table(taxblastred, "./08-blastdb/taxblast.mamfilt.txt",
            sep = "\t", row.names = FALSE, quote = FALSE)
###################################################################
setwd(prodir)

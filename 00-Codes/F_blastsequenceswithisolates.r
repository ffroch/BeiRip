blast <- read.table(("./05-qiime/HTS_16S/08-blastdb/comparison.blastout"), sep = "\t", header = FALSE)

colnames(blast) <- c("qseqid", "sseqid", "pident", "dontknow", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score")

# keep only rows with 100% identity
matches <- blast[blast$pident == 100,]
# keep only rows with alignment length >= 350
matches <- matches[matches$alignment_length >= 350,]
# count the number of matches for each query sequence
counts <- data.frame(table(matches$sseqid))
colnames(counts)[1] <- "sseqid"

# unique query sequences
matchesred <- data.frame("sseqid" = unique(matches$sseqid))
matchesred <- merge(matchesred, counts, by = "sseqid")

# save matchred as txt file
write.table(matches, file = "./03-Routputfiles/F_blastn_bac_seqvsiso.txt")
write.table(matchesred, file = "./03-Routputfiles/F_blastn_bac_seqvsiso_frequencies.txt", sep = "\t", row.names = FALSE, quote = FALSE)

####################################################################################################
# same for fungi
# import data of isolates
blast <- read.table("./05-qiime/HTS_18S/08-blastdb/comparison.mamfilt.blastout", sep = "\t", header = FALSE)
colnames(blast) <- c("qseqid", "sseqid", "pident", "dontknow", "alignment_length", "mismatches", "gap_opens", "q_start", "q_end", "s_start", "s_end", "evalue", "bit_score")

# keep only rows with 100% identity
matches <- blast[blast$pident == 100,]
# keep only rows with alignment length >= 250
matches <- matches[matches$alignment_length >= 250,]
# count the number of matches for each query sequence
counts <- data.frame(table(matches$sseqid))
colnames(counts)[1] <- "sseqid"

# unique query sequences
matchesred <- data.frame("sseqid" = unique(matches$sseqid))
matchesred <- merge(matchesred, counts, by = "sseqid")

# save matchred as txt file
write.table(matches, file = "./03-Routputfiles/F_blastn_fun_seqvsiso.txt", sep = "\t", row.names = FALSE, quote = FALSE)
write.table(matchesred, file = "./03-Routputfiles/F_blastn_fun_seqvsiso_frequencies.txt", sep = "\t", row.names = FALSE, quote = FALSE)

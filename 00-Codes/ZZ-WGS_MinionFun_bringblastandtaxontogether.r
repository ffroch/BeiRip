taxname <- read.table("taxid_to_name.txt", header=FALSE, sep="\t")

samplelist <- c("0367", "0368", "0991", "1109")

for (sample in samplelist) {
    acctax <- read.table(paste0("sample", sample, ".matched_taxids.txt"), header=FALSE, sep="\t")
    blast <- read.table(paste0("sample", sample, ".blastn"), header=FALSE, sep="\t")
    blastc <- dplyr::left_join(blast, acctax, by=c("V2"="V2"))
    blastcc <- dplyr::left_join(blastc, taxname, by=c("V3.y"="V1"))
    colnames(blastcc) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sseqid2", "taxid", "whatever", "taxname")
    write.table(blastcc, paste0("sample", sample, ".blastn.taxname"), quote=FALSE, row.names=FALSE, col.names=TRUE, sep="\t")
}
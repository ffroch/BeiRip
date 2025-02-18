setwd("/data/Unit_LMM/selberherr-group/roch/02_projects/BeiRip")

samplelist <- c("0367", "0368", "0991", "1109")

for (i in seq_along(samplelist)) {
    kraken_out <- read.table(paste0("./07-MinionFungi/08-kraken2/sample", samplelist[i], ".kraken2.out"), header = FALSE, sep = "\t")
    kraken_out <- kraken_out[, c(1, 2, 3, 4)]
    colnames(kraken_out) <- c("Status", "Contig", "TaxonomyID", "LengthAfterPolishing")

    kraken_report <- read.table(paste0("./07-MinionFungi/08-kraken2/sample", samplelist[i], ".report"), header = FALSE, sep = "\t")
    colnames(kraken_report) <- c("Percentage", "Reads", "taxonlevelreads", "Level", "TaxonomyID", "Taxon")

    flye_report <- read.table(paste0("./07-MinionFungi/05-assembly/sample", samplelist[i], "/assembly_info.txt"), header = FALSE, sep = "\t")
    colnames(flye_report) <- c("Contig", "Length", "Coverage", "circular", "repeat", "multiplicity", "alt_group", "graph_path")

    kraken_combined <- merge(kraken_out, kraken_report[, c("TaxonomyID", "Taxon")], by = "TaxonomyID")

    final_data <- merge(flye_report, kraken_combined, by.x = "Contig", by.y = "Contig")
    final_data <- final_data[order(final_data$Length, decreasing = TRUE), ]
    final_data$Taxon <- trimws(final_data$Taxon)

    write.table(final_data, paste0("./07-MinionFungi/08-kraken2/sample", samplelist[i], "_combined_output.txt"), sep = "\t", row.names = FALSE)
}

# remove unwanted contigs from the assembly

library(dplyr)
library(Biostrings)
library(tidyverse)
library(cowplot)
samplelist <- c("0367", "0368", "0991", "1109")
plotlist <- list()
for (i in seq_along(samplelist)) {
    data <- read.table(paste0("./07-MinionFungi/08-kraken2/sample", samplelist[i], "_combined_output.txt"), sep = "\t", header = TRUE)
    unwanted_taxa <- c("Lactococcus raffinolactis" , "Lactococcus lactis")
    unwanted_contigs <- data %>%
        filter(Taxon %in% unwanted_taxa) %>%
        pull(Contig)
    fasta_file <- paste0("./07-MinionFungi/06-polished2/sample", samplelist[i], "/consensus.fasta")
    output_file <- paste0("./07-MinionFungi/06-polished2/sample", samplelist[i], "/consensus_filtered.fasta")

    fasta <- readDNAStringSet(fasta_file)

    filtered_fasta <- fasta[!names(fasta) %in% unwanted_contigs]

    writeXStringSet(filtered_fasta, filepath = output_file)
    cat("Filtered FASTA file saved as:", output_file, "\n")

    # additionally calculate the gc content of the contigs
    removed_contigs <- fasta[names(fasta) %in% unwanted_contigs]
    removed_contigs_gc <- letterFrequency(removed_contigs, "GC", as.prob = TRUE)
    removed_contigs_gc <- data.frame(removed_contigs_gc)
    removed_contigs_gc$status <- "removed"
    removed_contigs_gc$name <- names(removed_contigs)
    kept_contigs <- fasta[!names(fasta) %in% unwanted_contigs]
    kept_contigs_gc <- letterFrequency(kept_contigs, "GC", as.prob = TRUE)
    kept_contigs_gc <- data.frame(kept_contigs_gc)
    kept_contigs_gc$status <- "kept"
    kept_contigs_gc$name <- names(kept_contigs)
    gc_data <- rbind(removed_contigs_gc, kept_contigs_gc)
    p <- ggplot(gc_data, aes(x = G.C, fill = status)) +
        geom_histogram() +
        theme_bw() +
        scale_fill_manual(values = c("darkred", "grey")) +
        labs(x = "GC content", y = "removed contigs") +
        theme(legend.position = c(0.05, 0.95), legend.justification = c(0, 1),
        legend.background = element_rect(fill = "white", color = "black"))
    plotlist[[i]] <- p
    ggsave(paste0("./07-MinionFungi/06-polished2/sample", samplelist[i], "/gc_content.png"), p)
}
plot <- do.call(plot_grid, c(plotlist, ncol = 2, labels = "auto"))
ggsave("./07-MinionFungi/06-polished2/gc_content.png", plot, width = 10, height = 10)
setwd(owd)

# get the GC content of "Lactococcus raffinolactis" that was sequenced on the same flowcell
lactostring <- readDNAStringSet("./06-Minion/08-polished/barcode04/consensus.fasta")
letterFrequency(lactostring, "GC", as.prob = TRUE)
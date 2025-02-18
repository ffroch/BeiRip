library(ggplot2)
library(sangeranalyseR)

#find all ab1 files from bacterial isolates (27F) and make list
filelist <- intersect(list.files(path = "./02-rawdata/SangerSeqData_bacteria",
                                 pattern = "*.ab1", full.names = TRUE),
                      list.files(path = "./02-rawdata/SangerSeqData_bacteria",
                                 pattern = "*27F", full.names = TRUE))

#skip lines steps if file StopRunningCompleteB1.txt exists
if (!file.exists("./02_StopRunningCompleteB1.txt")) {
  #make dirtory for outputfiles
  dir.create("./03-Routputfiles/B_SangerSeq/B1_trimmedSangerfasta_bac/")

  #generate DNA strings from all the ab1 files incl. trimming
  # and generate fasta files
  emptylist <- vector(mode = "list", length = length(filelist))
  i <- 0
  for (file in filelist){
    vegprod1 <- SangerRead(readFeature = "Forward Read",
                           readFileName = file,
                           TrimmingMethod = "M2",
                           M2CutoffQualityScore = 35,
                           M2SlidingWindowSize = 15,
                           M1TrimmingCutoff = NULL)
    if (vegprod1@QualityReport@trimmedMeanQualityScore > 35) {
      string <- vegprod1@primarySeqRaw
      trims <- vegprod1@QualityReport@trimmedStartPos
      trime <- vegprod1@QualityReport@trimmedFinishPos
      ifelse(trims < 20, trimseq <- subseq(string, start = 20, end = trime),
             trimseq <- subseq(string, start = trims, end = trime))
      print(trimseq)
      i <- i + 1
      emptylist[[i]] <- trimseq
      writeFastaSR(vegprod1,
                   outputDir = "./03-Routputfiles/B_SangerSeq/B1_trimmedSangerfasta_bac/")
    }
  }
  emptylist
  emptylist2 <- emptylist[-which(sapply(emptylist, is.null))]
  saveRDS(emptylist2, file = "./03-Routputfiles/B_SangerSeq/B1_SangerReadsList_bac.RData")
  #create an empty file called StopMakingDataStructure.txt
  text <- "This file stops B1_IsolateSeq_Import_Trimming_bacteria.R from running completely, if it exists. Delete this file to run A_MakeDataStructure.R"
  write(text, file = "./02_StopRunningCompleteB1.txt")
}
emptylist2 <- readRDS("./03-Routputfiles/B_SangerSeq/B1_SangerReadsList_bac.RData")


#make a list of all the fasta files
filelist2 <- list.files(path = "./03-Routputfiles/B_SangerSeq/B1_trimmedSangerfasta_bac/",
                        pattern = "*.fa", full.names = TRUE)
#generate namelist for the DNA strings
emptynamelist <- vector(mode = "list", length = length(filelist2))
j <- 0
for (file in filelist2) {
  j <- j + 1
  temporary <- str_remove(file, "./03-Routputfiles/B_SangerSeq/B1_trimmedSangerfasta_bac/")
  emptynamelist[[j]] <- str_remove(temporary, ".fa")
}
emptynamelist

# generate DNA stringset from the sequence list and name them with isolate ID
seqstrim <- DNAStringSet(unlist(emptylist2))
seqstrim
names(seqstrim) <- emptynamelist
seqstrim
seqstrim <- OrientNucleotides(seqstrim)


#check for Ns in the stringset, for later
# in the addSpecies Function Ns are not allowed
#check stringset for Ns
strings_with_n <- seqstrim[grepl("N", seqstrim)]
strings_with_n


#find the longest contiguous sequence parts in the stringsets and keep only that
find_longest_contiguous_part <- function(sequence) {
  parts <- unlist(strsplit(as.character(sequence), "N"))
  longest_part <- parts[which.max(nchar(parts))]
  return(longest_part)
}

longest_parts <- sapply(seqstrim, find_longest_contiguous_part)
longest_parts <- DNAStringSet(longest_parts)

#check how much of the sequence get lost, when
before <- Biostrings::width(seqstrim)
after <- Biostrings::width(longest_parts)
perc_loss <- 100 / before * (before - after)

#decision made:
# keep sequences without N as they are
# in sequences with N keep the longest contiguous part as
# long as the sequence is not shortend by more then 10
# percent of the original length
# in sequences with N and a large loss of the sequence when keeping
# only the longest sequence keep the complete sequence, but mark
# them in the trees and exclude them vom addSpecies

mcols(seqstrim)$perc_loss <- perc_loss
mcols(longest_parts)$perc_loss <- perc_loss

part1 <- longest_parts[mcols(longest_parts)$perc_loss < 10]
part2 <- seqstrim[mcols(seqstrim)$perc_loss >= 10]

seqstrim <- c(part1, part2)

writeXStringSet(seqstrim, "./03-Routputfiles/B_SangerSeq/B1_seqstrim_35_15_100bp_bac.fa")

#filter for trimmed sequences >100 bp
seqstrimclean <- seqstrim[nchar(seqstrim) >= 100]

writeXStringSet(seqstrimclean,
                "./03-Routputfiles/B_SangerSeq/B1_seqstrimclean_35_15_100bp_bac.fa")

#filter for trimmed sequences >100 bp
part1clean <- part1[nchar(part1) >= 100]
part2clean <- part2[nchar(part2) >= 100]
saveRDS(part1clean, "./03-Routputfiles/B_SangerSeq/B1_seqstrimclean_part1_bac.RDS")
saveRDS(part2clean, "./03-Routputfiles/B_SangerSeq/B1_seqstrimclean_part2_bac.RDS")

#make table with relevant infos from this step (summary data frame)
bsdf_c1 <- data.frame("parameter" = 0, "value" = 0)
bsdf_c1[1, ] <- c("number of ab1 files", length(filelist))
bsdf_c1[2, ] <- c("sequences passed quality filtering", length(emptylist2))
bsdf_c1[3, ] <- c("number of sequences", length(seqstrim))
bsdf_c1[4, ] <- c("number of sequences with N", length(strings_with_n))
bsdf_c1[5, ] <- c("number of sequences with N and >10% sequence length loss",
                  length(part2))
bsdf_c1[6, ] <- c("number of sequences with N and <10% sequence length loss",
                  length(part1))
bsdf_c1[7, ] <- c("number of sequences >100 bp", length(seqstrimclean))
bsdf_c1[8, ] <- c("mean sequence length of seqstrim",
                  round(mean(Biostrings::width(seqstrim)), digits = 0))
bsdf_c1[9, ] <- c("mean sequence length of seqstrimclean",
                  round(mean(Biostrings::width(seqstrimclean)), digits = 0))
bsdf_c1[10, ] <- c("median sequence length of seqstrim",
                   round(median(Biostrings::width(seqstrim)), digits = 0))
bsdf_c1[11, ] <- c("median sequence length of seqstrimclean",
                   round(median(Biostrings::width(seqstrimclean)), digits = 0))

seqslength <- data.frame("length" = Biostrings::width(seqstrimclean))
seqslength$length <- as.numeric(seqslength$length)
ggplot(seqslength, aes(x = length)) +
  geom_histogram()

bsdf_c1[12, ] <- c("sequences with <1000bp",
                   length(seqslength[seqslength$length < 1000, ]))

#make bsdf_c1$value numeric
bsdf_c1$value <- as.numeric(bsdf_c1$value)
#add column named index for latex DTLfetch
bsdf_c1$index <- c(1:12)

print(bsdf_c1)

#save bsdf_c1 as txt
write.table(bsdf_c1, file = "./03-Routputfiles/B_SangerSeq/B1_statisticalsummary_bac.txt",
            sep = "\t", row.names = FALSE, col.names = TRUE)
saveRDS(bsdf_c1, file = "./ZZ-FilesForManuscript/B1_statisticalsummary_bac.rds")

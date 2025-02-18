library(ggplot2)
library(sangeranalyseR)


#find all ab1 files from bacterial isolates (27F) and make list
filelist <- intersect(list.files(path = "./02-rawdata/SangerSeqData_fungi",
                                 pattern = "*.ab1", full.names = TRUE),
                      list.files(path = "./02-rawdata/SangerSeqData_fungi",
                                 pattern = "*EK528", full.names = TRUE))



#skip lines steps if file StopRunningCompleteC1.txt exists
if (!file.exists("./02_StopRunningCompleteC1.txt")) {
  #make dirtory for outputfiles
  dir.create("./03-Routputfiles/B_SangerSeq/B2_trimmedSangerfasta_fun/")

  # generate DNA strings from all the ab1 files incl. trimming and
  # generate fasta files
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
                   outputDir =  "./03-Routputfiles/B_SangerSeq/B2_trimmedSangerfasta_fun/")
    }
  }
  emptylist
  emptylist2 <- emptylist[-which(sapply(emptylist, is.null))]
  saveRDS(emptylist2, file = "./03-Routputfiles/B_SangerSeq/B2_SangerReadsList_fun.RData")
  #create an empty file called StopMakingDataStructure.txt
  text <- "This file stops B2_IsolateSeq_Import_Trimming_fungi.R from running completely, if it exists. Delete this file to run A_MakeDataStructure.R"
  write(text, file = "./03_StopRunningCompleteC2.txt")
}

emptylist2 <- readRDS("./03-Routputfiles/B_SangerSeq/B2_SangerReadsList_fun.RData")

#make a list of all the fasta files
filelist2 <- list.files(path = "./03-Routputfiles/B_SangerSeq/B2_trimmedSangerfasta_fun/",
                        pattern = "*.fa", full.names = TRUE)
#generate namelist for the DNA strings
emptynamelist <- vector(mode = "list", length = length(filelist2))
j <- 0
for (file in filelist2) {
  j <- j + 1
  temporary <- str_remove(file, "./03-Routputfiles/B_SangerSeq/B2_trimmedSangerfasta_fun/")
  emptynamelist[[j]] <- str_remove(temporary, ".fa")
}
emptynamelist

# generate DNA stringset from the sequence list and name them with isolate ID
seqstrim <- DNAStringSet(unlist(emptylist2))
seqstrim
names(seqstrim) <- emptynamelist
seqstrim
seqstrim <- OrientNucleotides(seqstrim)


#check for Ns in the stringset, for later in the addSpecies
# Function Ns are not allowed
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
# in sequences with N keep the longest contiguous part as long as the
# sequence is not shortend by more then 10 percent of the originall length
# in sequences with N and a large loss of the sequence when keeping only
# the longest sequence keep the complete sequence, but mark them in the
# trees and exclude them vom addSpecies

mcols(seqstrim)$perc_loss <- perc_loss
mcols(longest_parts)$perc_loss <- perc_loss

part1 <- longest_parts[mcols(longest_parts)$perc_loss < 10]
part2 <- seqstrim[mcols(seqstrim)$perc_loss >= 10]

seqstrim <- c(part1, part2)

writeXStringSet(seqstrim, "./03-Routputfiles/B_SangerSeq/B2_seqstrim_35_15_100bp_fun.fa")


#filter for trimmed sequences >100 bp
seqstrimclean <- seqstrim[nchar(seqstrim) >= 100]

writeXStringSet(seqstrimclean,
                "./03-Routputfiles/B_SangerSeq/B2_seqstrimclean_35_15_100bp_fun.fa")

#filter for trimmed sequences >100 bp
part1clean <- part1[nchar(part1) >= 100]
part2clean <- part2[nchar(part2) >= 100]
saveRDS(part1clean, "./03-Routputfiles/B_SangerSeq/B2_seqstrimclean_part1_fun.RDS")
saveRDS(part2clean, "./03-Routputfiles/B_SangerSeq/B2_seqstrimclean_part2_fun.RDS")

#make table with relevant infos from this step (summary data frame)
fsdf_c2 <- data.frame("parameter" = 0, "value" = 0)
fsdf_c2[1, ] <- c("number of ab1 files", length(filelist))
fsdf_c2[2, ] <- c("sequences passed quality filtering", length(emptylist2))
fsdf_c2[3, ] <- c("number of sequences", length(seqstrim))
fsdf_c2[4, ] <- c("number of sequences with N", length(strings_with_n))
fsdf_c2[5, ] <- c("number of sequences with N and >10% sequence length loss",
                  length(part2))
fsdf_c2[6, ] <- c("number of sequences with N and <10% sequence length loss",
                  length(part1))
fsdf_c2[7, ] <- c("number of sequences >100 bp", length(seqstrimclean))
fsdf_c2[8, ] <- c("mean sequence length of seqstrim",
                  round(mean(Biostrings::width(seqstrim)), digits = 0))
fsdf_c2[9, ] <- c("mean sequence length of seqstrimclean",
                  round(mean(Biostrings::width(seqstrimclean)), digits = 0))
fsdf_c2[10, ] <- c("median sequence length of seqstrim",
                   round(median(Biostrings::width(seqstrim)), digits = 0))
fsdf_c2[11, ] <- c("median sequence length of seqstrimclean",
                   round(median(Biostrings::width(seqstrimclean)), digits = 0))

seqslength <- data.frame("length" = Biostrings::width(seqstrimclean))
seqslength$length <- as.numeric(seqslength$length)
ggplot(seqslength, aes(x = length)) +
  geom_histogram()

fsdf_c2[12, ] <- c("sequences with <900bp",
                   length(seqslength[seqslength$length < 900, ]))

#make fsdf_c2$value numeric
fsdf_c2$value <- as.numeric(fsdf_c2$value)

print(fsdf_c2)
#add column named index for latex DTLfetch
fsdf_c2$index <- c(1:12)

#save fsdf_c2 as txt
write.table(fsdf_c2,
            file = "./03-Routputfiles/B_SangerSeq/B2_statisticalsummary_fun.txt",
            sep = "\t",
            row.names = FALSE, col.names = TRUE)
#save fsdf_c2 as csv
saveRDS(fsdf_c2, file = "./ZZ-FilesForManuscript/B2_statisticalsummary_fun.rds")

#Run R script A_MakeDataStructure.R if file StopMakingDataStructure.txt does not exist	
 if (!file.exists("./StopMakingDataStructure.txt")){
   dir.create("03-Routputfiles")
#create an empty file called StopMakingDataStructure.txt
text<-"This file stops A_MakeDataStructure.R from running, if it exists. Delete this file to run A_MakeDataStructure.R"
write(text, file = "./01_StopMakingDataStructure.txt")
}

#Run R script A1_ColonyFormingUnits.R
# script summarizes the plate count data of the bacteria and fungi
# further it generates Supplementary Figure 1
source("./00-Codes/A1_ColonyFormingUnits.R")

#Run R script B1_IsolateSeq_Import_Trimming_bacteria.R
# script imports the sanger sequences of the bacteria and trims them
source("./00-Codes/B1A_IsolateSeq_Import_Trimming_bacteria.R")

#Run R script B1B_IsolateSeq_TaxonomyR_bacteria.R
# assigns taxonomy to the trimmed sequences of the bacteria
source("./00-Codes/B1B_IsolateSeq_Import_TaxonomyR_bacteria.R")

#Run R script B2_IsolateSeq_Import_Trimming_fungi.R
# script imports the sanger sequences of the fungi and trims them
source("./00-Codes/B2A_IsolateSeq_Import_Trimming_fungi.R")

#Run R script B2B_IsolateSeq_TaxonomyR_fungi.R
# assigns taxonomy to the trimmed sequences of the fungi
source("./00-Codes/B2B_IsolateSeq_TaxonomyR_fungi.R")

#Run R script B3_GenerateListOfDoneIsolates.R
# script generates a list of the isolates that have been trimmed and assigned taxonomy
# further it is the bases for the Supplementary Tables 1 - 4
source("./00-Codes/B3_GenerateListOfDoneIsolates.R")

#Run R script C1_Isolates_Figure_P1 - P4 .R
# script generates the data for the isolate tree and the Supplementary Figure 8
source("./00-Codes/C1_Isolates_Figure_P1.R")
source("./00-Codes/C2_Isolates_Figure_P2.R")
source("./00-Codes/C3_Isolates_Figure_P3.R")
source("./00-Codes/C4_Isolates_Figure_P4.R")

# Run R script D1A_PhyloObjectsBac.r
# imports the raw qiime data of the bacteria and creates a phyloseq object
# final phyloobject contains the raw reads inlcuding the spike in bacteria. 
# two files are created, one with and one without mock community samples and negative controls
source("./00-Codes/D1A_PhyloObjectsBac.r")

# Run R script D1B_PhyloObjectsBac.r
# takes the phyloseq object from D1A_PhyloObjectsBac.r and calculates the relative abundances. and calculates the absolute cell numbers per sample based on the spike in bacteria
# final phyloobject contains the relative abundance of the 16S rRNA gene sequencing
# two files are created, one with and one without mock community samples and negative controls
source("./00-Codes/D1B_PhyloObjectsBac.r")

# Run R script D1C_PhyloObjectsBac.r
# takes the phyloseq object from D1A_PhyloObjectsBac.r and calculates the copy number corrected version of the 16S rRNA gene sequencing
# final phyloobject contains the copy number corrected version of the 16S rRNA gene sequencing
# two files are created, one with and one without mock community samples and negative controls
source("./00-Codes/D1C_PhyloObjectsBac.r")

# Run R script D1D_PhyloObjectsBac.r
# takes the phyloseq object from D1B_PhyloObjectsBac.r and calculates the absolute abundances based on Imtechella spike in  cells
# final phyloobject contains the absolute abundance of the 16S rRNA gene sequencing
# two files are created, on with the spike in bacteria and one without
source("./00-Codes/D1D_PhyloObjectsBac.r")

# Run R script D1E_PhyloObjectsBac.r
# takes the phyloseq object from D1D_PhyloObjectsBac.r and calculates the copy number corrected version of the 16S rRNA gene sequencing
# final phyloobject contains the copy number corrected version of the 16S rRNA gene sequencing
source("./00-Codes/D1E_PhyloObjectsBac.r")

# Run R script D1F_PhyloObjectsBac.r
# takes the phyloseq object fomr D1D_PhyloObjectsBac.r, calculates the copy number corrected version and aggregates the data to the species and genus level
source("./00-Codes/D1F_PhyloObjectsBac.r")


# Run R script D2A_PhyloObjectsFun.r
# imports the raw qiime data of the fungi and creates a phyloseq object
# final phyloobject contains the raw reads
# two files are created, one with and one without negative controls
source("./00-Codes/D2A_PhyloObjectsFun.r")

# Run R script D2B_PhyloObjectsFun.r
# takes the phyloseq object from D2A_PhyloObjectsFun.r and calculates the relative abundances
# final phyloobject contains the relative abundance of the 18S rRNA gene sequencing
# two files are created, one with and one without negative controls
source("./00-Codes/D2B_PhyloObjectsFun.r")

# Run R script Z-18S-qPCR-primerdesign-sequencingBased.r
# checks the primers from primer blast 
# creates Supplementary Figure 37
source("./00-Codes/Z-18S-qPCR-primerdesign-sequencingBased.r")

# Run R script Z-18S-qPCR-results.r
# imports the qPCR data and prepares it for the calculation of the absolute abundances
source("./00-Codes/Z-18S-qPCR-results.r")

# Run R script D2C_PhyloObjectsFun.r
# takes the phyloseq object from D2A_PhyloObjectsFun.r and calculates the absolute abundances based on qPCR data
# final phyloobject contains the absolute abundance of the 18S rRNA gene sequencing
# two files are created, one with and one without negative controls
source("./00-Codes/D2C_PhyloObjectsFun.r")

# Run R script D2D_PhyloObjectsFun.r
# takes the phyloseq object from D2C_PhyloObjectsFun.r and aggregates the data to the genus level
source("./00-Codes/D2D_PhyloObjectsFun.r")

# Run R script D3A_PhyloObjectsComrel.r
# loads the bacterial and fungal absolute abundance phyloseq objects
# transform it into relative abundance data
# rescale the values using Z-score scaling
# combine the data into one phyloseq object
# and save the data for CoNet and Cytoscape, which is used for Figure 2a
source("./00-Codes/D3A_PhyloObjectsComrel.r")

# Run r script E_CountsVsGenecopies.r
# compares the counts and growth curves of the CFU data with the gene copy numbers of the 16S and 18S rRNA gene sequencing
# creates Supplementary Figures 2 and 3
source("./00-Codes/E_CountsVsGenecopies.r")

# run R script F_balstsequenceswithisolates.r
# combines the isolate sequences with the 16S and 18S rRNA gene sequencing data
sourece("./00-Codes/F_blastsequenceswithisolates.r")

# run R script ZZ-helper_BarPlotcolorDefiner.r
# creates the color definitions for the bar plots
source("./00-Codes/ZZ-helper_BarPlotcolorDefiner.r")

# run R script G1_TaxaPlotsBac.r
# creates the taxa plots of the 16S rRNA gene sequencing data
# creates Figure 1a, and Supplementary Figure 4
source("./00-Codes/G1_TaxaPlotsBac.r")

# run R script G2_TaxaPlotsFun.r
# creates the taxa plots of the 18S rRNA gene sequencing data
# creates Figure 1c, and Supplementary Figure 9
source("./00-Codes/G2_TaxaPlotsFun.r")

# run R script H1_GrowthCurvesBac.r
# creates the growth curves of the 16S rRNA gene sequencing data
# creates heatmaps of the 16S rRNA gene sequencing data
# creates Supplementary Figures 5 - 7
source("./00-Codes/H1_GrowthCurvesBac.r")

# run R script H2_GrowthCurvesFun.r
# creates the growth curves of the 18S rRNA gene sequencing data
# creates heatmaps of the 18S rRNA gene sequencing data
# create Supplementary Figures 10 and 11
source("./00-Codes/H2_GrowthCurvesFun.r")

# run R script I1-AlphadiversityBac.r
# calculates the alpha diversity of the 16S rRNA gene sequencing data
# creates Supplementary Figure 12 and Fig 1 d
source("./00-Codes/I1-AlphadiversityBac.r")

# run R script I2-AlphadiversityFun.r
# calculates the alpha diversity of the 18S rRNA gene sequencing data
# creates Supplementary Figure 13 and Fig 1 e
source("./00-Codes/I2-AlphadiversityFun.r")

# run R script J-AvgDist_PCoA.r
# calculates the beta diversity of the 16S and 18S rRNA gene sequencing data
# creates the PCoA plots of the 16S and 18S rRNA gene sequencing data
# creates Figure 1b
source("./00-Codes/J-AvgDist_PCoA.r")

# run R script S_qPCR_CoCulture.r
# calculates the qPCR data of the co-culture experiments
# creates figure 2b
source("./00-Codes/S_qPCR_CoCulture.r")



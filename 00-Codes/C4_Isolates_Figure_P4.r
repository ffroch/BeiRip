############################################################################################################
### FIGURE 1 ###
#Figure 1 shows the following:
#1. A tree with one representative of each Genus based on the 16S/18S rRNA gene sequences, representative is isolate with the longest sequence in the most prevalent Cluster within a Genus
#2. Additional information about the isolates and their diversity within a Genus
##2a. column for each batch with points
##2b. points are sized by the grouped number of 16S/18S rRNA gene genotypes
##2c. points are coloured by the last slice the isolate was found in
#only sequences >1000bp with taxonomy >=2 matching databases were used
############################################################################################################
### THIS IS PART 4 - PLOTTING ###
############################################################################################################
library(ggplot2)
library(DECIPHER)
library(dplyr)
library(forcats)
library(ggtree)
library(phylobase)
library(stringr)
library(phylotools)
library(ape)
library(RColorBrewer)
library(viridis)
library(cowplot)
################################################################
#### FIGURE 1A - BACTERIA ####
### Load data ###
# tree data
bactreephylo<-read.tree("./03-Routputfiles/C_IsolateTree/C2_bactree.txt")#if phyml does not work

# meta data
dataex<-readRDS("./03-Routputfiles/C_IsolateTree/C3_bacmeta.RDS")
### Prepare data #######################################
# make dataex pivot_wider
dataexred<-dataex[,c(1:2,4:ncol(dataex))]
dataexwide<-dataexred%>%
    tidyr::pivot_wider(names_from = batch, values_from = c(latesttp, earliesttp, presence, diversity))
dataexwide<-data.frame(dataexwide)
row.names(dataexwide)<-dataexwide$newGenus
# make diversity a ranged factor
dataexwide$diversity<-cut(dataexwide$diversity_C1, breaks=c(0,1,5,10,50), labels=c("<2", "2-5", "6-10", ">10"))
### combine tree and metadata ###########################
# remove \" from tip.label to have same names in tree and metadata
bactreephylo$tip.label<-gsub("\"","",bactreephylo$tip.label)
# make phylo4d object from tree 
g1<-as(bactreephylo, "phylo4")
# combine phylo4d object with metadata
bacset<-phylo4d(g1, dataexwide, missing.data="warn")
### make plot ##########################################
# make color brewer palette from viridis with 4 colors
pal <- viridis_pal(option="viridis")(4)
# add column aesthetics to plot
d1 <- data.frame(x = seq(1.7, 2.1, length.out = 3),
                lab = c("batch 1", "batch 2", "batch 3"))
#first draft to get node labels
win.graph()   
p1 <- ggtree(bacset, color = "grey")+ 
    geom_text(aes(label=node), hjust=-.3)+
    geom_tippoint(aes(col = diversity, size=latesttp_C1, alpha=presence_C1), x = d1$x[1], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C2, alpha=presence_C2), x = d1$x[2], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C3, alpha=presence_C3), x = d1$x[3], shape = 16) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C1, alpha=presence_C1), x = d1$x[1], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C2, alpha=presence_C2), x = d1$x[2], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C3, alpha=presence_C3), x = d1$x[3], shape = 21) + 
    scale_size_continuous(breaks=c(1,5,10), labels=c("day 0", "day 15", "day 85"), range = c(1,5), name="earliest (bordered)\nand latest isolation") + 
    scale_color_manual(values=pal,name="no. of 16S/18S rRNA\ngene genotypes") +
    scale_alpha(guide="none") +
    scale_shape_manual(values=c(16,21), name="")+
    geom_tiplab(offset = 0.2, fontface=3) + xlim(0, 4) + vexpand(0.05, -1)+ 
    theme(axis.text.x = element_blank(),legend.position = c(0.85,0.8))+
    guides(color = guide_legend(override.aes = list(size = 5)))
    
p1

# force alpha bacteria to higher node
print(bacset@edge)
bacset@edge[31, 1] <- 41
print(bacset@edge)

# final plot
p2 <- ggtree(bacset, color = "grey")+ 
    geom_tippoint(aes(col = diversity, size=latesttp_C1, alpha=presence_C1), x = d1$x[1], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C2, alpha=presence_C2), x = d1$x[2], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C3, alpha=presence_C3), x = d1$x[3], shape = 16) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C1, alpha=presence_C1), x = d1$x[1], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C2, alpha=presence_C2), x = d1$x[2], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C3, alpha=presence_C3), x = d1$x[3], shape = 21) +
    geom_cladelabel(node=67, label="Bacilli", offset=2, offset.text = 0.05, align=T, fontface = "italic") +
    geom_cladelabel(node=59, label="Actinobacteria", offset=2, offset.text = 0.05,extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=17, label="Deinococci", offset=2, offset.text = 0.05, extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=56, label="Alphaproteobacteria", offset=2, offset.text = 0.05,extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=42, label="Gammaproteobacteria", offset=2, offset.text = 0.05,extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=77, label="Bacteroidia", offset=2, offset.text = 0.05, extend = 0.4,align=T, fontface = "italic") +
    scale_size_continuous(breaks=c(1,5,10), labels=c("day 0", "day 15", "day 85"), range = c(1,5), name="earliest (bordered)\nand latest isolation") + 
    scale_color_manual(values=pal,name="no. of 16S rRNA\ngene genotypes") +
    scale_alpha(guide="none") +
    scale_shape_manual(values=c(16,21), name="")+
    geom_tiplab(offset = 0, fontface=3, align=TRUE) + xlim(0, 4) + vexpand(0.05, -1)+ 
    theme(axis.text.x = element_blank(),legend.position = c(0.9,0.8))+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    #geom_treescale(x=0.1, y=-0.3, offset = -1.1, color="grey")+
    #geom_text(aes(0.1, y=-0.3, label = "substitutions/site"), vjust = 2, hjust=-0.3, color = "grey")+
    geom_text(aes(x = x, y = -0.9, label = lab), data = d1, angle = 90)
p2
ggsave("./03-Routputfiles/C_IsolateTree/C4_bactreephylo_Fig1A.png", width=20, height=20, units="cm", dpi=300)
dev.off()

############################################################################################################
#### FIGURE 1B - FUNGI ####
### Load data ##########################################
# tree data
#funtreephylo<-read.tree("./03-Routputfiles/D2_funtree.txt")#if phyml does not work
funtreephylo<-read.tree("./03-Routputfiles/C_IsolateTree/C2_funtree.txt")
# meta data
datafex<-readRDS("./03-Routputfiles/C_IsolateTree/C3_funmeta.RDS")
### Prepare data #######################################
# make datafex pivot_wider
datafexred<-datafex[,c(1:2,4:ncol(datafex))]
datafexwide<-datafexred%>%
    tidyr::pivot_wider(names_from = batch, values_from = c(latesttp, earliesttp, presence, diversity))
datafexwide<-data.frame(datafexwide)
#replace unclassified in newtax with uncl.
datafexwide$newtax<-gsub("unclassified","uncl.",datafexwide$newtax)
row.names(datafexwide)<-datafexwide$newtax
# make diversity a ranged factor
datafexwide$diversity<-cut(datafexwide$diversity_C1, breaks=c(0,1,5,10,50), labels=c("<2", "2-5", "6-10", ">10"))
### combine tree and metadata ###########################
# remove \" from tip.label to have same names in tree and metadata
funtreephylo$tip.label<-gsub("\"","",funtreephylo$tip.label)
# add space between unclassified and taxon
funtreephylo$tip.label<-gsub("unclassified","uncl. ",funtreephylo$tip.label)
# make phylo4d object from tree
g2<-as(funtreephylo, "phylo4")
# combine phylo4d object with metadata
funset<-phylo4d(g2, datafexwide, missing.data="warn")
### make plot ##########################################
# make color brewer palette from viridis with 4 colors
pal <- viridis_pal(option="viridis")(4)

#first draft to get node labels
win.graph()
# add column aesthetics to plot
d2 <- data.frame(x = seq(1.8,2.1, length.out = 3),
                lab = c("batch 1", "batch 2", "batch 3"))
p3 <- ggtree(funset, color = "grey")+ 
    geom_text(aes(label=node), hjust=-.3)+
    geom_tippoint(aes(col = diversity, size=latesttp_C1, alpha=presence_C1), x = d2$x[1], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C2, alpha=presence_C2), x = d2$x[2], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C3, alpha=presence_C3), x = d2$x[3], shape = 16) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C1, alpha=presence_C1), x = d2$x[1], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C2, alpha=presence_C2), x = d2$x[2], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C3, alpha=presence_C3), x = d2$x[3], shape = 21) + 
    scale_size_continuous(breaks=c(1,5,10), labels=c("day 0", "day 15", "day 85"), range = c(1,5), name="earliest (bordered)\nand latest isolation") + 
    scale_color_manual(values=pal,name="no. of 16S/18S rRNA\ngene genotypes") +
    scale_alpha(guide="none") +
    scale_shape_manual(values=c(16,21), name="")+
    geom_tiplab(offset = 0.2, fontface=3) + xlim(0, 1) + vexpand(0.05, -1)+ 
    theme(axis.text.x = element_blank(),legend.position = c(0.85,0.8))+
    guides(color = guide_legend(override.aes = list(size = 5)))
p3    
# final plot 1B
p4 <- ggtree(funset, color = "grey")+ 
    geom_tippoint(aes(col = diversity, size=latesttp_C1, alpha=presence_C1), x = d2$x[1], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C2, alpha=presence_C2), x = d2$x[2], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C3, alpha=presence_C3), x = d2$x[3], shape = 16) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C1, alpha=presence_C1), x = d2$x[1], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C2, alpha=presence_C2), x = d2$x[2], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C3, alpha=presence_C3), x = d2$x[3], shape = 21) +
    geom_cladelabel(node=30, label="Ascomycota", offset=2, offset.text = 0.05, align=T, fontface = "italic") +
    geom_cladelabel(node=24, label="Basidiomycota", offset=2, offset.text = 0.05,extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=21, label="Mucormycota", offset=2, offset.text = 0.05, extend = 0.4, align=T, fontface = "italic") +
    scale_size_continuous(breaks=c(1,5,10), labels=c("day 0", "day 15", "day 85"), range = c(1,5), name="earliest (bordered)\nand latest isolation") + 
    scale_color_manual(values=pal,name="no. of 16S rRNA\ngene genotypes") +
    scale_alpha(guide="none") +
    scale_shape_manual(values=c(16,21), name="")+
    geom_tiplab(offset = 0, fontface=3, align=TRUE) + xlim(0, 4) + vexpand(0.05, -1)+ 
    theme(axis.text.x = element_blank(),legend.position = c(0.9,0.8))+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    geom_treescale(x=0.1, y=-0.3, width=0.1, offset = -0.7, color="grey")+
    geom_text(aes(0.1, y=-0.3, label = "substitutions/site"), vjust = 1.9, hjust=-0.3, color = "grey")+
    geom_text(aes(x = x, y = -0.8, label = lab), data = d2, angle = 90)
p4
ggsave("./03-Routputfiles/C_IsolateTree/C4_funtreephylo_Fig1B.png", width=20, height=17, units="cm", dpi=300)

################################################################################
#### FIGURE 1 - Plot adaptions for combined plotting ####
#combine plot for bacteria and fungi with cowplot
# add column aesthetics to plot 
d3 <- data.frame(x = seq(1.6,1.9, length.out = 3),
                lab = c("#1", "#2", "#3"))
# Plot 1A
p5 <- ggtree(bacset, color = "grey")+ 
    geom_tippoint(aes(col = diversity, size=latesttp_C1, alpha=presence_C1), x = d1$x[1], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C2, alpha=presence_C2), x = d1$x[2], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C3, alpha=presence_C3), x = d1$x[3], shape = 16) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C1, alpha=presence_C1), x = d1$x[1], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C2, alpha=presence_C2), x = d1$x[2], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C3, alpha=presence_C3), x = d1$x[3], shape = 21) +
    geom_cladelabel(node=67, label="Bacilli", offset=2, offset.text = 0.05, align=T, fontface = "italic") +
    geom_cladelabel(node=59, label="Actinobacteria", offset=2, offset.text = 0.05,extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=17, label="Deinococci", offset=2, offset.text = 0.05, extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=56, label="Alphaproteobacteria", offset=2, offset.text = 0.05,extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=42, label="Gammaproteobacteria", offset=2, offset.text = 0.05,extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=77, label="Bacteroidia", offset=2, offset.text = 0.05, extend = 0.4,align=T, fontface = "italic") +
    scale_size_continuous(breaks=c(1,5,10), labels=c("day 0", "day 15", "day 85"), range = c(1,5), name="earliest (bordered)\nand latest isolation") + 
    scale_color_manual(values=pal,name="no. of 16S rRNA\ngene genotypes") +
    scale_alpha(guide="none") +
    scale_shape_manual(values=c(16,21), name="")+
    geom_tiplab(offset = 0, fontface=3, align=TRUE) + xlim(0, 4) + vexpand(0.05, -1)+ 
    theme(axis.text.x = element_blank(),legend.position = c(0.9,0.8))+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    #geom_treescale(x=0.1, y=-0.3, width = 0.1, offset = -0.9, color="grey")+
    #geom_text(aes(0.1, y=-0.3, label = "substitutions/site"), vjust = 1.5, hjust=-0.3, color = "grey")+
    geom_text(aes(2.5, 39, label = "Class"))
# plot 1B
p6 <- ggtree(funset, color = "grey")+ 
    geom_tippoint(aes(col = diversity, size=latesttp_C1, alpha=presence_C1), x = d3$x[1], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C2, alpha=presence_C2), x = d3$x[2], shape = 16) + 
    geom_tippoint(aes(col = diversity, size=latesttp_C3, alpha=presence_C3), x = d3$x[3], shape = 16) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C1, alpha=presence_C1), x = d3$x[1], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C2, alpha=presence_C2), x = d3$x[2], shape = 21) + 
    geom_tippoint(aes(col = NULL, fill= NULL, size=earliesttp_C3, alpha=presence_C3), x = d3$x[3], shape = 21) +
    geom_cladelabel(node=30, label="Ascomycota", offset=1.95, offset.text = 0.05, align=T, fontface = "italic") +
    geom_cladelabel(node=24, label="Basidiomycota", offset=1.95, offset.text = 0.05,extend = 0.4, align=T, fontface = "italic") +
    geom_cladelabel(node=21, label="Mucormycota", offset=1.95, offset.text = 0.05, extend = 0.4, align=T, fontface = "italic") +
    scale_size_continuous(breaks=c(1,5,10), labels=c("day 0", "day 15", "day 85"), range = c(1,5), name="earliest (bordered)\nand latest isolation") + 
    scale_color_manual(values=pal,name="no. of 18S rRNA\ngene genotypes") +
    scale_alpha(guide="none") +
    scale_shape_manual(values=c(16,21), name="")+
    geom_tiplab(offset = 0, fontface=3, align=TRUE) + xlim(0, 4) + vexpand(0.05, -1)+ 
    theme(axis.text.x = element_blank(),legend.position = c(0.9,0.8))+
    guides(color = guide_legend(override.aes = list(size = 5)))+
    #geom_treescale(x=0.1, y=-0.3, width=0.1, offset = -0.7, color="grey")+
    #geom_text(aes(0.1, y=-0.3, label = "substitutions/site"), vjust = 1.5, hjust=-0.3, color = "grey")+
    geom_text(aes(1.3, y=-0.3, label = "batch"), vjust = 1.5, hjust=0, color = "black")+
    geom_text(aes(x = x, y = -0.75, label = lab), data = d3, angle = 0, size=4)+
    geom_text(aes(2.5, 22, label = "Phylum"))


plot_grid(p5, p6, ncol=1, align="v", rel_heights=c(1.4,1), labels = "auto")

ggsave("./03-Routputfiles/C_IsolateTree/C4_Figure1_AB.png", width=21, height=29.5, units="cm")
ggsave("./03-Routputfiles/C_IsolateTree/C4_Figure1_AB.svg", width=21, height=29.5, units="cm")

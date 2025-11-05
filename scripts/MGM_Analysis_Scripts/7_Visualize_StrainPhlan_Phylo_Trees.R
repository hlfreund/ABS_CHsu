#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("")

suppressPackageStartupMessages({ # load packages quietly
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  #library(scales)
  library(grid)
  library(ape)
  library(plyr)
  library(dplyr)
  library(viridis)
  library(readxl)
  library(metagenomeSeq)
  library(DESeq2)
  library(dplyr)
  library(magrittr)
  library(MASS)
  library(dendextend)
  library(tidyr)
  library(reshape)
  library(reshape2)
  library(wesanderson)
  library(shades)
  library(ALDEx2)
  library(rstatix)
  library(devtools)
  library(decontam)
  library(openxlsx)
  library(Maaslin2)
  library(stringr)
  library(forcats)
  library(NatParksPalettes)
  library("colorspace")
  library(ggtree)
})

#### Import ABS Project Metadata ####

abs.meta.basic<-as.data.frame(read_xlsx("SRA_submission/SRA_Metagenomes_metadata.xlsx", sheet="SRA_data"))
head(abs.meta.basic) # using the original ABS metadata Cynthia provided me to match up with sequencing ID names used in relab table
# ^^ match these data by patient ID

abs.meta.basic<-subset(abs.meta.basic, select=c("sample_name","library_ID","filename"))
abs.meta.basic$filename<-gsub("_R1_001.fastq.gz","",abs.meta.basic$filename)

abs.categories<-read.table("data/Metagenomes/Revisions_3.3.2025/StrainPhlan_Results/MGM_Condition_Metadata.txt",header=TRUE,sep="\t")

abs.meta.cats<-merge(abs.meta.basic,abs.categories,by="sample_name")

# create group color palette to merge with abs.metadata
grp.clrs2 = as.data.frame(t(data.frame("Flare"="#ffa600","Remission"="#ff6361","HHP"="#007561")))

grp.clrs2$condition<-rownames(grp.clrs2)
colnames(grp.clrs2)[which(names(grp.clrs2) == "V1")] <- "Group_Color"
grp.clrs2

abs.metadata<-merge(abs.meta.cats,grp.clrs2,by="condition")
head(abs.metadata)

#### Load Phylogenetic (Clade-Specific) Trees ####
## strainphlan used clade markers to create MSA --> tree (from metaphlan results)

# Ecoli clades tree
ecoli.tree<-read.tree("data/Metagenomes/Revisions_3.3.2025/StrainPhlan_Results/RAxML_result.t__SGB10068.StrainPhlAn4.tre")
ecoli.tree$tip.label
#ecoli.tree$tip.label<-gsub("_L001","",ecoli.tree$tip.label)
ggplot(ecoli.tree) + geom_tree() + theme_tree() + geom_tiplab()

# Kpneumoniae clades tree
kpneum.tree<-read.tree("data/Metagenomes/Revisions_3.3.2025/StrainPhlan_Results/RAxML_bestTree.t__SGB10115_group.StrainPhlAn4.tre")
#kpneum.tree$tip.label<-gsub("_L001","",kpneum.tree$tip.label)
ggplot(kpneum.tree) + geom_tree() + theme_tree() + geom_tiplab()

#### Attach Ecoli Metadata to Ecoli Tree ####

ecoli.meta<-tibble(abs.metadata[abs.metadata$filename %in% ecoli.tree$tip.label,])
head(ecoli.meta)
colnames(ecoli.meta)

# add on rows of NAs for next step
# ecoli.meta[nrow(ecoli.meta)+20,] <- NA

# renaming this column "label" to merge with tree
colnames(ecoli.meta)[colnames(ecoli.meta)=="filename"]<-"label" 
#ecoli.meta$label<-gsub("_L001","",ecoli.meta$label)

# turn tree into tibble for joining
ecoli.tree2<-as_tibble(ecoli.tree)
#ecoli.tree2$label<-gsub("_L001","",ecoli.tree2$label)
head(ecoli.tree2)

# full_join to merge the tbl_tree with metadata df
ecoli.all<-full_join(ecoli.tree2,ecoli.meta,by=join_by("label"))
#ecoli.tree$tip.label<-gsub("_L001","",ecoli.tree$tip.label)

#### Plot Ecoli Phylo Tree ####
# plot it all!
ecoli.tree.fig1<-ggtree(as.phylo(ecoli.tree2)) %<+% ecoli.all + geom_tiplab(size=6) + geom_tippoint(aes(color=condition),size=3) +
  coord_cartesian(clip = 'off') + 
  scale_color_manual(values=c(Flare="#ffa600",Remission="#ff6361",HHP="#007561")) +
  theme(plot.title = element_text(colour = "black", size=25),text = element_text(colour = "black", size=22)) +
  labs(title = "Phylo Tree - E.coli Strains Across Samples",legend.text = element_text(size = 20, colour = "black"),colour="Condition")

# may want to give the donor a different color ?
ggsave(ecoli.tree.fig1,filename = "figures/PhyloTrees/ABS_MGMs_Ecoli_CladeMarkers_samestr_phylotree.png", width = 4000,height = 2000,units = "px",dpi = 200,create.dir = TRUE)

ecoli.tree.fig2<-msaplot(ecoli.tree.fig1, fasta="data/Metagenomes/Revisions_3.3.2025/StrainPhlan_Results/t__SGB10068.StrainPhlAn4_concatenated.aln", offset=0.3, bg_line = TRUE,window=c(150, 500))+
  labs(subtitle="Multiple Sequence Alignment, 150 - 500 window")

ggsave(ecoli.tree.fig2,filename = "figures/PhyloTrees/ABS_MGMs_Ecoli_CladeMarkers_samestr_phylotree_MSA.png", width = 6000,height = 2000,units = "px",dpi = 200,create.dir = TRUE)


#### Attach Kpneumoniae Metadata to Kpneumoniae Tree ####

kpneum.meta<-tibble(abs.metadata[abs.metadata$filename %in% kpneum.tree$tip.label,])
head(kpneum.meta)
colnames(kpneum.meta)

# add on rows of NAs for next step
# kpneum.meta[nrow(kpneum.meta)+20,] <- NA

# renaming this column "label" to merge with tree
colnames(kpneum.meta)[colnames(kpneum.meta)=="filename"]<-"label" 

# turn tree into tibble for joining
kpneum.tree2<-as_tibble(kpneum.tree)

# full_join to merge the tbl_tree with metadata df
kpneum.all<-full_join(kpneum.tree2,kpneum.meta,by=join_by("label"))

#### Plot Kpneumoniae Phylo Tree ####
# plot it all!
kpneum.tree.fig1<-ggtree(as.phylo(kpneum.tree2)) %<+% kpneum.all + geom_tiplab(size=7) + geom_tippoint(aes(color=condition),size=3) +
  coord_cartesian(clip = 'off') + 
  theme_tree2(plot.margin=margin(6, 120, 6, 6)) + scale_color_manual(values=c(Flare="#ffa600",Remission="#ff6361",HHP="#007561")) +
  theme(plot.title = element_text(colour = "black", size=25),
        text = element_text(colour = "black", size=22)) +
  labs(title = "Phylo Tree - K.pneumoniae Strains Across Samples",legend.text = element_text(size = 20, colour = "black"),colour="Condition")
# may want to give the donor a different color ?

ggsave(kpneum.tree.fig1,filename = "figures/PhyloTrees/ABS_MGMs_Kpneumoniae_CladeMarkers_samestr_phylotree.png", width = 4000,height = 2000,units = "px",dpi = 200,create.dir = TRUE)

kpneum.tree.fig2<-msaplot(kpneum.tree.fig1, fasta="data/Metagenomes/Revisions_3.3.2025/StrainPhlan_Results/t__SGB10115_group.StrainPhlAn4_concatenated.aln", offset=0.3, bg_line = TRUE)+
  labs(subtitle="Multiple Sequence Alignment, 150 - 500 window")

ggsave(kpneum.tree.fig2,filename = "figures/PhyloTrees/ABS_MGMs_Kpneumoniae_CladeMarkers_samestr_phylotree_MSA.png", width = 6000,height = 2000,units = "px",dpi = 200,create.dir = TRUE)


#### Load FMT/Donor-Specific Phylogenetic (Clade-Specific) Trees ####
## strainphlan used clade markers to create MSA --> tree (from metaphlan results)

# Ecoli clades tree
ecoli.fmt.tree<-read.tree("data/Metagenomes/Revisions_3.3.2025/StrainPhlan_Results_FMT/RAxML_result.t__SGB10068.StrainPhlAn4.tre")
ecoli.fmt.tree$tip.label
#ecoli.fmt.tree$tip.label<-gsub("_L001","",ecoli.fmt.tree$tip.label)
ggplot(ecoli.fmt.tree) + geom_tree() + theme_tree() + geom_tiplab()

#### FMT - Attach Ecoli Metadata to Ecoli Tree ####

# ecoli.meta<-tibble(abs.metadata[abs.metadata$filename %in% ecoli.tree$tip.label,])
# head(ecoli.meta)
# colnames(ecoli.meta)

# add on rows of NAs for next step
# ecoli.meta[nrow(ecoli.meta)+20,] <- NA

# renaming this column "label" to merge with tree
# colnames(ecoli.meta)[colnames(ecoli.meta)=="filename"]<-"label" 
#ecoli.meta$label<-gsub("_L001","",ecoli.meta$label)

# turn tree into tibble for joining
ecoli.fmt.tree2<-as_tibble(ecoli.fmt.tree)
#ecoli.tree2$label<-gsub("_L001","",ecoli.tree2$label)
head(ecoli.fmt.tree2)

# full_join to merge the tbl_tree with metadata df
ecoli.fmt.all<-full_join(ecoli.fmt.tree2,ecoli.meta,by=join_by("label"))
#ecoli.tree$tip.label<-gsub("_L001","",ecoli.tree$tip.label)

#### Plot Ecoli Phylo Tree ####
# plot it all!
ecoli.fmt.tree.fig1<-ggtree(as.phylo(ecoli.fmt.tree2)) %<+% ecoli.fmt.all + geom_tiplab(size=6) + geom_tippoint(aes(color=condition),size=3) +
  coord_cartesian(clip = 'off') + 
  scale_color_manual(values=c(Flare="#ffa600",Remission="#ff6361",HHP="#007561")) +
  theme(plot.title = element_text(colour = "black", size=25),text = element_text(colour = "black", size=22)) +
  labs(title = "Phylo Tree - E.coli Strains Across FMT vs Donor Samples",legend.text = element_text(size = 20, colour = "black"),colour="Condition")

# may want to give the donor a different color ?
ggsave(ecoli.fmt.tree.fig1,filename = "figures/PhyloTrees/ABS_FMT_Only_MGMs_Ecoli_CladeMarkers_samestr_phylotree.png", width = 4000,height = 2000,units = "px",dpi = 200,create.dir = TRUE)

ecoli.fmt.tree.fig2<-msaplot(ecoli.fmt.tree.fig1, fasta="data/Metagenomes/Revisions_3.3.2025/StrainPhlan_Results_FMT/t__SGB10068.StrainPhlAn4_concatenated.aln", offset=0.3, bg_line = TRUE,window=c(150, 500))+
  labs(subtitle="Multiple Sequence Alignment, 150 - 500 window")

ggsave(ecoli.fmt.tree.fig2,filename = "figures/PhyloTrees/ABS_FMT_Only_MGMs_Ecoli_CladeMarkers_samestr_phylotree_MSA.png", width = 6000,height = 2000,units = "px",dpi = 200,create.dir = TRUE)

ecoli.fmt.tree.fig3<-msaplot(ecoli.fmt.tree.fig1, fasta="data/Metagenomes/Revisions_3.3.2025/StrainPhlan_Results_FMT/t__SGB10068.StrainPhlAn4_concatenated.aln", offset=0.3, bg_line = TRUE)+
  labs(subtitle="Multiple Sequence Alignment")

ggsave(ecoli.fmt.tree.fig3,filename = "figures/PhyloTrees/ABS_FMT_Only_MGMs_Ecoli_CladeMarkers_samestr_phylotree_MSA_complete.png", width = 7000,height = 2000,units = "px",dpi = 200,create.dir = TRUE)

#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  library(grid)
  library(ape)
  library(plyr)
  library(dplyr)
  library(ggbiplot)
  library(readxl)
  library(dplyr)
  library(tidyr)
  library(reshape)
  library(reshape2)
  library(NatParksPalettes)
  library(shades)
  library(ggvegan)
  library(corrplot)
  library(forcats)
  library(factoextra)
  library(stats)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/ITS2_DADA2/ABS_ITS2_DataReady.Rdata") # save global env to Rdata file

its2.ASV_table[1:4,1:4]
its2.ASV_table[(nrow(its2.ASV_table)-4):(nrow(its2.ASV_table)),(ncol(its2.ASV_table)-4):(ncol(its2.ASV_table))] # last 4 rows & cols
head(metadata)

# drop Donor sample
metadata.v2<-metadata[!grepl("DonorSamp", metadata$SampleID),]
head(metadata.v2)

# drop Donor sample from ASV table
its2.ASV_table.v2<-its2.ASV_table[!grepl("DonorSamp", its2.ASV_table$SampleID),]
its2.ASV_table.v2[,1:4]

#### List of Important R Objects for Reference ####

its2.ASV_table[1:5,1:5] 
# ^ this is the updated ASV table where singletons have been removed
## sample names were also updated
its2.ASV_table.v2[1:5,1:5] # same as ASV table above but FMT donor sample dropped

head(metadata) # updated metadata that includes Flare group colors
head(metadata.v2) # same as metadata above but FMT donor sample dropped
its2.ASV_meta[1:5,1:10] # combined metadata + ASV counts

head(its2.ASV_tax.clean) # updated taxa table with singleton ASVs removed
head(its2.AllTax.melt) # melted ASVs + taxonomy table -- for relative abundance calcs

head(its2.ALL.dat) # ALL data combined
# ^ metadata, ASVs + counts, & taxonomy data

#### CLR Transform All Comp Data ####
rownames(its2.ASV_table.v2)
its2.ASV_table.v2[1:4,1:4]

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
f.clr<-decostand(its2.ASV_table.v2[,-1],method = "clr", pseudocount = 1) #CLR transformation
f.clr[1:4,1:4]

#### PCA of ASV Composition - Beta Diversity ####

# check rownames of CLR transformed ASV data & metadata.v2
# f.clr has SampleIDs as rownames, columns are ASV IDs
# metadata.v2 is metadata excluding the Donor FMT sample

rownames(f.clr) %in% rownames(metadata.v2)
metadata.v2=metadata.v2[rownames(f.clr),] ## reorder metadata.v2 to match order of CLR data

# Use princomp() to run PCA of Euclidean distances of CLR data (i.e., still Aitchison distances) so we can get to species contributions!
## create PCA of CLR-transformed counts (this will use Euclidean distance and gives the same results as the PCoA; you can compare based on the PC axes %)
f.pca<-prcomp(f.clr)
summary(f.pca) # proportion of variance gives us variance explained by each axis; should match your pcoa results if you're using Aitchison distance
f.pca$rotation # gives the contributions of each ASV to the variance/similarity of samples across the Euclidean space

#### Visualize PCA & Contributing Variables ####

# visualize PCA with contribution variables and the samples
## color scale corresponds to variable (here ASVs) contributions to PC axes' variation
fviz_pca_var(f.pca,
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,angle=45),panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white"))
ggsave(filename = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_PCA_ContributingASVs_biplot.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

# PCA biplot that shows contributing variables (ASVs) and samples in different colors
fviz_pca_biplot(f.pca,
                col.var = "contrib", # Variables color
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                col.ind = "#696969") + theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,angle=45),panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white"))
ggsave(filename = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_PCA_ContributingASVs_Samples_biplot.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

# view the contributions of the variables; aka visualizing f.pca$rotations in a meaningful way
## ASV contributions for PC1
fviz_contrib(f.pca,choice="var",fill="lightslateblue",color=NA,top=50, axes = 1)+theme_minimal() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,angle=45),panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white"))+
  xlab("Top 50 Contributing ASVs - PC Axis 1") + ylab("Contributions (%)")
ggsave(filename = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_PCA_Top50ContributingASVs_PC1.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

## ASV contributions for PC2
fviz_contrib(f.pca,choice="var",fill="lightslateblue",color=NA,top=50, axes = 2)+theme_minimal() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,angle=45),panel.background = element_rect(fill = "white", colour = "white"),
        plot.background = element_rect(fill = "white"))+
  xlab("Top 50 Contributing ASVs - PC Axis 2") + ylab("Contributions (%)")
ggsave(filename = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_PCA_Top50ContributingASVs_PC2.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

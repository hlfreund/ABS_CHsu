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
  library(pairwiseAdonis)
  library(plotly)
  library(htmlwidgets)
})

#### Import ABS + HC Combined abs.all.basic.meta ####

# this abs.all.basic.meta is a file I created using the SRA submission file + added propensity-matched HCs to spreadsheet
abs.all.basic.meta<-as.data.frame(read_xlsx("data/Metagenomes/SRA_Metagenomes_metadata_basic.xlsx", sheet="AllSamples_List"))
head(abs.all.basic.meta)

abs.all.basic.meta$SampleID<-gsub("_L00[0-9]_R1_001.fastq.gz","",abs.all.basic.meta$filename)

# create group color palette to merge with abs.all.basic.meta
grp.clrs = as.data.frame(t(data.frame("HC"="#40557B","Flare"="#ffa600","Remission"="#ff6361","HHP"="#007561")))
grp.clrs

grp.clrs$Condition<-rownames(grp.clrs)
colnames(grp.clrs)[which(names(grp.clrs) == "V1")] <- "Group_Color"
grp.clrs

abs.all.basic.meta<-merge(abs.all.basic.meta,grp.clrs,by="Condition")
unique(abs.all.basic.meta$Group_Color)

rownames(abs.all.basic.meta)<-abs.all.basic.meta$SampleID

#### Import SameSTR Results ####

# samestr cooccurrences file
abs.samestr.cooccur <- as.data.frame(read.delim('data/Metagenomes/Revisions_3.3.2025/ABS_HC_SameSTR_Results/samestr_out_summarize/ABS_HC_sstr_cooccurrences.tsv', sep="\t"))
head(abs.samestr.cooccur)
abs.samestr.cooccur$row<-gsub("_L00[0-9]","",abs.samestr.cooccur$row) # remove seq file info from sample ID
abs.samestr.cooccur$col<-gsub("_L00[0-9]","",abs.samestr.cooccur$col) # remove seq file info from sample ID
abs.samestr.cooccur[1:4,]
colnames(abs.samestr.cooccur)[which(names(abs.samestr.cooccur) == "row")] <- "SampleID"
# drop A022_10_S12 sample - excluded from SRA upload
abs.samestr.cooccur<-subset(abs.samestr.cooccur,SampleID!="A022_10_S12" & col!="A022_10_S12")

# IMPORTANT NOTE FOR LATER
# SameSTR results are such that not all SampleIDs appear in either SampleID or col column (see below)
# aka, there is one sampleID in the rows that is not in the columns, and one sampleID in the column that is not in the rows
unique(abs.samestr.cooccur$SampleID) %in% unique(abs.samestr.cooccur$col) # one false in samestr cooccurrence data
unique(abs.samestr.cooccur$col) %in% unique(abs.samestr.cooccur$SampleID) # one false in samestr cooccurrence data
# this will matter when you need to generate PCoA, which considers the rows as samples

# samestr events file
abs.samestr.events <- as.data.frame(read.delim('data/Metagenomes/Revisions_3.3.2025/ABS_HC_SameSTR_Results/samestr_out_summarize/ABS_HC_sstr_strain_events.tsv', sep="\t"))
head(abs.samestr.events)
abs.samestr.events$row<-gsub("_L00[0-9]","",abs.samestr.events$row) # remove seq file info from sample ID
abs.samestr.events$col<-gsub("_L00[0-9]","",abs.samestr.events$col) # remove seq file info from sample ID
abs.samestr.events$clade<-gsub("t__","",abs.samestr.events$clade) # remove leading "t__" from clade names
colnames(abs.samestr.events)[which(names(abs.samestr.events) == "row")] <- "SampleID"
abs.samestr.events[1:4,]
# drop A022_10_S12 sample - excluded from SRA upload
abs.samestr.events<-subset(abs.samestr.events,SampleID!="A022_10_S12" & col!="A022_10_S12")

abs.cooccur.meta<-merge(abs.samestr.cooccur,abs.all.basic.meta,by="SampleID")
# sanity checks: note that $SampleID and $col do not all contain same sample IDs due to nature of dataframe and sample to sample comparisons
length(unique(abs.cooccur.meta$SampleID))
unique(abs.cooccur.meta$SampleID) %in% abs.all.basic.meta$SampleID
length(unique(abs.cooccur.meta$col))
unique(abs.cooccur.meta$col) %in% abs.all.basic.meta$SampleID

abs.events.meta<-merge(abs.samestr.events,abs.all.basic.meta,by="SampleID")
# sanity checks: note that $SampleID and $col do not all contain same sample IDs due to nature of dataframe and sample to sample comparisons

abs.cooccur.meta$Condition<-factor(abs.cooccur.meta$Condition,levels=c("HC","Flare","Remission","HHP"))
abs.events.meta$Condition<-factor(abs.events.meta$Condition,levels=c("HC","Flare","Remission","HHP"))

#### Visualize SameSTR Results ####

# heatmap with genus & species, faceted by condition
hm.test1<-ggplot(abs.cooccur.meta, aes(SampleID, col, fill= shared_strain)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=25,
                       breaks=c(0,10,20,30,40),
                       labels=c("0","10","20","30","40"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="", y="", title="Shared Strains by Sample + Condition",fill="# of Shared Strains") + facet_grid(Condition~., scales = "free",space="free") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(hm.test1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/SameSTR_Figs/ABS_MGMs_SameSTR_shared_strains_heatmap.png",width = 9800,height = 9500,units = "px",dpi = 200,create.dir = TRUE)

# #### Shared Strains within Clades per Sample ####
# head(abs.events.meta)
# abs.sharedstrains.meta<-abs.events.meta[abs.events.meta$event=="shared_strain",]
# head(abs.sharedstrains.meta)
# col1.strs<-unique(abs.sharedstrains.meta[,c(1,3)])
# head(col1.strs)
# 
# col2.strs<-unique(abs.sharedstrains.meta[,c(2,3)])
# colnames(col2.strs)[which(colnames(col2.strs)=="col")] <- "SampleID"
# 
# shared.clades<-unique(rbind(col1.strs,col2.strs))
# head(shared.clades)
# 
# shared.clades.table<-dcast(shared.clades,SampleID~clade,length)
# rownames(shared.clades.table)<-shared.clades.table$SampleID
# 
# shared.clades.meta<-merge(abs.all.basic.meta,shared.clades.table,by="SampleID")
# 
# # PCOA w/ Jaccard distance matrix (of binary data)
# #pcoa(vegdist(abs.cooccur.table.sym,na.rm = TRUE,binary=TRUE,method="jaccard"),correction="lingoes") # 
# abs.shared.jac.pcoa1 <- pcoa(dist(shared.clades.table[,-1],method="binary"),correction="lingoes") # 
# #save.image("data/Amplicon/FT_16SV3V4_CLR_EucDist_Ready.Rdata")
# 
# # The proportion of variances explained is in its element values$Relative_eig
# abs.shared.jac.pcoa1$values
# 
# # extract principal coordinates
# abs.shared.jac.pcoa1.vectors<-data.frame(abs.shared.jac.pcoa1$vectors)
# abs.shared.jac.pcoa1.vectors$SampleID<-rownames(abs.shared.jac.pcoa1$vectors)
# 
# # merge pcoa coordinates w/ abs.all.basic.meta
# abs.shared.jac.pcoa1.meta<-merge(abs.shared.jac.pcoa1.vectors, abs.all.basic.meta, by.x="SampleID", by.y="SampleID")
# 
# head(abs.shared.jac.pcoa1.meta)
# rownames(abs.shared.jac.pcoa1.meta)<-abs.shared.jac.pcoa1.meta$SampleID
# 
# head(abs.shared.jac.pcoa1$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# # PC1 = 4.74%, PC2 = 2.15%
# 
# #### Visualize PCoA ####
# 
# ## *** all figures that end in _PCOA1 come from the same single PCoA
# #data is just subsetted to understand how points are related to each other w/in & across timepoints
# 
# # create PCoA ggplot fig
# 
# pcoa1<-ggplot(abs.shared.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Condition)) +geom_point(size=4)+theme_bw()+
#   labs(color="Condition",title="PCoA of Presence/Absence of Shared Strains (SameSTR Results)",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
#         axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Condition",values=unique(abs.shared.jac.pcoa1.meta$Group_Color[order(abs.shared.jac.pcoa1.meta$Condition)])) +
#   scale_fill_manual(name ="Condition",values=unique(abs.shared.jac.pcoa1.meta$Group_Color[order(abs.shared.jac.pcoa1.meta$Condition)])) +
#   xlab("PC1 [4.74%]") + ylab("PC2 [2.15%]")+ annotate("text",x=0.40,y=0.5,label="PERMANOVA, p.adj = 0.0003")
# 
# ggsave(pcoa1,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedStrains_Binary_PCOA1.png", 
#        width = 2000,height = 1500,units = "px",
#        dpi = 200,create.dir = TRUE)
# 
# pcoa2<-ggplot(abs.shared.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Condition)) +geom_point(size=4)+theme_bw()+
#   labs(color="Condition",title="PCoA of Presence/Absence of Shared Strains (SameSTR Results)",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
#         axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Condition",values=unique(abs.shared.jac.pcoa1.meta$Group_Color[order(abs.shared.jac.pcoa1.meta$Condition)])) +
#   scale_fill_manual(name ="Condition",values=unique(abs.shared.jac.pcoa1.meta$Group_Color[order(abs.shared.jac.pcoa1.meta$Condition)])) +
#   xlab("PC1 [4.74%]") + ylab("PC2 [2.15%]")+ stat_ellipse(alpha = 0.5) + annotate("text",x=0.80,y=0.55,label="PERMANOVA, p.adj = 0.0003")
# 
# ggsave(pcoa2,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedStrains_Binary_Ellipses_PCOA2.png", 
#        width = 2000,height = 1500,units = "px",
#        dpi = 200,create.dir = TRUE)
# 
# pcoa3<-ggplot(abs.shared.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Condition)) +geom_text(label=abs.shared.jac.pcoa1.meta$SampleID,position="identity")+theme_bw()+
#   labs(color="Condition",title="PCoA of Presence/Absence of Shared Strains (SameSTR Results)",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
#         axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Condition",values=unique(abs.shared.jac.pcoa1.meta$Group_Color[order(abs.shared.jac.pcoa1.meta$Condition)])) +
#   scale_fill_manual(name ="Condition",values=unique(abs.shared.jac.pcoa1.meta$Group_Color[order(abs.shared.jac.pcoa1.meta$Condition)])) +
#   xlab("PC1 [4.74%]") + ylab("PC2 [2.15%]")+ annotate("text",x=0.40,y=0.5,label="PERMANOVA, p.adj = 0.0003")
# 
# ggsave(pcoa3,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedStrains_Binary_Labeled_by_SampleID_PCOA3.png", 
#        width = 2000,height = 1500,units = "px",
#        dpi = 200,create.dir = TRUE)
# 
# # 3D PCoA
# 
# pltly.all.a<-plot_ly(abs.shared.jac.pcoa1.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~Condition, colors = c(unique(abs.shared.jac.pcoa1.meta$GroupName_Color[order(abs.shared.jac.pcoa1.meta$Condition)])),
#                      symbol=~Condition) %>%
#   layout(scene = list(xaxis = list(title = 'PC1 17.61%'),
#                       yaxis = list(title = 'PC2 13.91%'),
#                       zaxis = list(title = 'PC3 12.07%')))
# 
# 
# #### Generate PCoA with SameSTR Results ####
# # my logic behind this PCoA, since we are comparing samples to each other...
# # we use PCoAs to determine how similar samples are to each other based on shared taxa or functions...
# # here, we are saying that who they share strain's with (i.e., each other) can be represented as taxa...
# # so if Samples A, B, C, and D share strains with E, they will be closer to each other in PCoA
# # if samples share a lot of strains with each other, they will be closer to each other on the PCoA
# # self-to-self comparisons are assigned 0s as to not falsely weight a sample being similar to itself
# 
# head(abs.cooccur.meta)
# abs.cooccur.meta$shared_strain<-as.numeric(abs.cooccur.meta$shared_strain)
# abs.cooccur.table1<-dcast(abs.cooccur.meta,col~SampleID,value.var = "shared_strain")
# abs.cooccur.table1[1:5,1:5] # sanity check
# rownames(abs.cooccur.table1)<-abs.cooccur.table1$col
# 
# abs.cooccur.table2<-dcast(abs.cooccur.meta,SampleID~col,value.var = "shared_strain")
# abs.cooccur.table2[1:5,1:5] # sanity check
# rownames(abs.cooccur.table2)<-abs.cooccur.table2$SampleID
# 
# abs.cooccurr.melt1<-melt(abs.cooccur.table1,by="col")
# head(abs.cooccurr.melt1)
# names(abs.cooccurr.melt1)[which(names(abs.cooccurr.melt1)=="col")]<-"SampleID"
# 
# abs.cooccurr.melt2<-melt(abs.cooccur.table2,by="SampleID")
# head(abs.cooccurr.melt2)
# #names(abs.cooccurr.melt2)[which(names(abs.cooccurr.melt2)=="variable")]<-"col"
# 
# abs.cooccur.rbind<-rbind(abs.cooccurr.melt1,abs.cooccurr.melt2)
# abs.cooccur.rbind[is.na(abs.cooccur.rbind)]<-0
# # drop Donor now for sanity's sake
# abs.cooccur.rbind<-subset(abs.cooccur.rbind,SampleID!="Donor41_33_S77" & variable!="Donor41_33_S77")
# abs.cooccur.final.table<-dcast(abs.cooccur.rbind,SampleID~variable,value.var = "value",sum)
# rownames(abs.cooccur.final.table)<-abs.cooccur.final.table$SampleID
# abs.cooccur.final.table[1:5,1:5]
# # 
# # # To continue, we have to make a symmetric table...
# # # aka, there is one sampleID in the rows that is not in the columns, and one sampleID in the column that is not in the rows
# # unique(abs.cooccur.meta$SampleID) %in% unique(abs.cooccur.meta$col) # one false in samestr cooccurrence data
# # unique(abs.cooccur.meta$col) %in% unique(abs.cooccur.meta$SampleID) # one false in samestr cooccurrence data
# # 
# # # problem: abs.cooccur.table is not symmetric, causing us to lose data points in PCOA...
# # ## code to the rescue from here: https://stackoverflow.com/questions/12591708/from-asymmetric-matrix-or-dataframe-into-a-symmetric-square-matrix-with-r
# # 
# # # use newly created cooccurr table to get the names as a vector
# # NAMES <- sort(unique(c(rownames(abs.cooccur.table[,-1]), colnames(abs.cooccur.table[,-1]))))
# # # create temp data frame using stack() to keep the rows we would lose with dcast to make table
# # temp <- data.frame(rows = rownames(abs.cooccur.table[,-1]), stack(abs.cooccur.table[,-1]))
# # 
# # # factor the row and column names using NAMES vector
# # temp$rows <- factor(temp$rows, NAMES)
# # temp$ind <- factor(temp$ind, NAMES)
# # 
# # # Use xtabs to get your desired output of symmetric table. Wrap it in as.data.frame.matrix to get a data.frame as output
# # abs.cooccur.table.sym<-as.data.frame.matrix(xtabs(values ~ rows + ind, temp))
# # abs.cooccur.table.sym[1:5,1:5]
# 
# # drop rows where all elements are NA
# ## rowSums(is.na()) counts # of NAs, then !=ncol() sees if the # of NAs is less than the total # of columns
# # abs.cooccur.table2<-abs.cooccur.table[,-1][rowSums(is.na(abs.cooccur.table[,-1])) != ncol(abs.cooccur.table[,-1]),]
# # abs.cooccur.table[is.na(abs.cooccur.table)]<-0
# 
# # PCOA w/ Jaccard distance matrix (of binary data)
# #pcoa(vegdist(abs.cooccur.table.sym,na.rm = TRUE,binary=TRUE,method="jaccard"),correction="lingoes") # 
# abs.jac.pcoa1 <- pcoa(dist(abs.cooccur.final.table2[,-1],method="binary"),correction="lingoes") # 
# #save.image("data/Amplicon/FT_16SV3V4_CLR_EucDist_Ready.Rdata")
# 
# # The proportion of variances explained is in its element values$Relative_eig
# abs.jac.pcoa1$values
# 
# # extract principal coordinates
# abs.jac.pcoa1.vectors<-data.frame(abs.jac.pcoa1$vectors)
# abs.jac.pcoa1.vectors$SampleID<-rownames(abs.jac.pcoa1$vectors)
# 
# # merge pcoa coordinates w/ abs.all.basic.meta
# abs.jac.pcoa1.meta<-merge(abs.jac.pcoa1.vectors, abs.all.basic.meta, by.x="SampleID", by.y="SampleID")
# 
# head(abs.jac.pcoa1.meta)
# rownames(abs.jac.pcoa1.meta)<-abs.jac.pcoa1.meta$SampleID
# 
# head(abs.jac.pcoa1$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# # PC1 = 4.74%, PC2 = 2.15%
# 
# #### Homogeneity of Variance & PERMANOVA tests - Shared Strains & Clades - Composition by Groups ####
# ## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# # multivariate analogue to Levene's test of homogeneity of variances
# # program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA
# 
# #While PERMANOVA tests differences in group means (analogous to MANOVA),
# ## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
# #(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
# ## * need a distance matrix!
# 
# head(abs.all.basic.meta)
# head(abs.cooccur.table[,-1])
# 
# # are both dfs in the same order? Need this for stats below
# shared.clade.str.meta<-abs.all.basic.meta[which(rownames(abs.all.basic.meta) %in% rownames(shared.clades.table)),]
# rownames(shared.clade.str.meta) %in% rownames(shared.clades.table)
# shared.clade.str.meta=shared.clade.str.meta[rownames(shared.clades.table),] ## reorder metadata to match order of shared strain-clade data
# 
# # drop Donor from this specific PCoA
# shared.clades.table1<-subset(shared.clades.table,shared.clades.table$SampleID!="Donor41_33_S77")
# shared.clade.str.meta1<-subset(shared.clade.str.meta,shared.clade.str.meta$SampleID!="Donor41_33_S77")
# 
# # first by compare dispersions by sampling date
# abs.disper1<-betadisper((dist(shared.clades.table1,method="binary")), shared.clade.str.meta1$Condition)
# abs.disper1
# 
# ## Significant differences in homogeneities can be tested using either parametric or permutational tests,
# ##and parametric post hoc contrasts can also be investigated:
# 
# permutest(abs.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
# #Pairwise comparisons:
# #  (Observed p-value below diagonal, permuted p-value above diagonal)
# #           Flare      HC     HHP Remission
# # Flare             0.51900 0.12700     0.511
# # HC        0.53095         0.74400     0.920
# # HHP       0.14247 0.74052             0.615
# # Remission 0.54058 0.91787 0.61171
# 
# anova(abs.disper1) # p = 0.8153 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across conditions
# # ANOVA adjusted p-value
# aov.beta.p1<-anova(abs.disper1)[["Pr(>F)"]] # get p values from ANOVA
# p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))
# 
# TukeyHSD(abs.disper1) # tells us which category's dispersion MEANS are significantly different than each other
# #                     diff       lwr       upr     p adj
# # HC-Flare        -0.026802964 -0.12978475 0.07617882 0.9053709
# # HHP-Flare       -0.041564030 -0.15541467 0.07228661 0.7775292
# # Remission-Flare -0.022086501 -0.13131517 0.08714217 0.9525217
# # HHP-HC          -0.014761066 -0.11888080 0.08935867 0.9827436
# # Remission-HC     0.004716463 -0.09432824 0.10376117 0.9993152
# # Remission-HHP    0.019477529 -0.09082467 0.12977973 0.9675878
# 
# # If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# # If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# # Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself
# ## adonis2 will fail if there are samples with rowSums = 0, let's drop them
# abs.cooccur.table.sym1[is.na(abs.cooccur.table.sym1)]<-0
# abs.cooccur.table.sym2<-subset(abs.cooccur.table.sym1,rowSums(abs.cooccur.table.sym1)!=0)
# abs.all.basic.meta2<-subset(abs.all.basic.meta1,rownames(abs.all.basic.meta1) %in% rownames(abs.cooccur.table.sym2))
# # went from 130 --> 90 samples, so 38 samples had rowSums = 0 aka no shared strains with other samples identified by SameStr
# 
# # use the pvalue below for PCOA that excludes Donor
# pnova1<-adonis2(shared.clades.table1[,-1] ~ Condition,data=shared.clade.str.meta1,method = "jaccard",by="terms",permutations= 10000)
# pnova1 # p-value = 9.999e-05
# p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# # 0.00029997
# 
# ##one issue with adonis is that it doesn't do multiple comparisons *******
# # tells us that something is different, but what is different? Which sample/plot/location?
# ## our four provinces differ, but do all of them differ,or just one?
# 
# ##random person on the internet to the rescue!
# #install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
# 
# abs.shrstr.dist = dist(shared.clades.table1[,-1],method= "binary") #Jaccard distance matrix
# pair.mod1<-pairwise.adonis(abs.shrstr.dist,shared.clade.str.meta1$Condition, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
# pair.mod1
# #        pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# # 1        Flare vs HC  1  1.557065 3.966016 0.09681234   1e-04     0.0006  **
# # 2       Flare vs HHP  1  1.254836 2.935974 0.06124782   2e-04     0.0012   *
# # 3 Flare vs Remission  1  1.622496 3.781513 0.07302824   1e-04     0.0006  **
# # 4          HC vs HHP  1  1.500746 3.899418 0.09306616   1e-04     0.0006  **
# # 5    HC vs Remission  1  1.346660 3.453688 0.07769181   1e-04     0.0006  **
# # 6   HHP vs Remission  1  1.240358 2.937042 0.05655005   1e-04     0.0006  **
# 
# # Visualize dispersions
# par(mar=c(1,1,1,1))
# png('figures/SameSTR_Figs/ABS_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
# plot(abs.disper1,main = "Centroids and Dispersion based on Jaccard Distance", col=grp.clrs$Group_Color)
# dev.off()
# 
# par(mar=c(1,1,1,1))
# png('figures/SameSTR_Figs/ABS_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
# boxplot(abs.disper1,xlab="By GroupName", main = "Distance to Centroid by Category", sub="Based on Jaccard Distance", col=grp.clrs$GroupName_Color)
# dev.off()
# 
# #### Homogeneity of Variance & PERMANOVA tests - Composition by Groups ####
# ## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# # multivariate analogue to Levene's test of homogeneity of variances
# # program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA
# 
# #While PERMANOVA tests differences in group means (analogous to MANOVA),
# ## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
# #(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
# ## * need a distance matrix!
# 
# head(abs.all.basic.meta)
# head(abs.cooccur.final.table[,-1])
# 
# # are both dfs in the same order? Need this for stats below
# abs.sub.meta<-abs.all.basic.meta[which(rownames(abs.all.basic.meta) %in% rownames(abs.cooccur.final.table)),]
# rownames(abs.sub.meta) %in% rownames(abs.cooccur.final.table)
# abs.sub.meta=abs.sub.meta[rownames(abs.cooccur.final.table),] ## reorder metadata to match order of shared strain-clade data
# 
# # first by compare dispersions by sampling date
# abs.disper1<-betadisper((dist(abs.cooccur.final.table[,-1],method="binary")), abs.sub.meta$Condition)
# abs.disper1
# 
# ## Significant differences in homogeneities can be tested using either parametric or permutational tests,
# ##and parametric post hoc contrasts can also be investigated:
# 
# permutest(abs.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
# #Pairwise comparisons:
# #  (Observed p-value below diagonal, permuted p-value above diagonal)
# #           Flare      HC     HHP Remission
# # Flare             0.51900 0.12700     0.511
# # HC        0.53095         0.74400     0.920
# # HHP       0.14247 0.74052             0.615
# # Remission 0.54058 0.91787 0.61171          
# 
# anova(abs.disper1) # p = 0.8153 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across conditions
# # ANOVA adjusted p-value
# aov.beta.p1<-anova(abs.disper1)[["Pr(>F)"]] # get p values from ANOVA
# p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))
# 
# TukeyHSD(abs.disper1) # tells us which category's dispersion MEANS are significantly different than each other
# #                     diff       lwr       upr     p adj
# # HC-Flare        -0.026802964 -0.12978475 0.07617882 0.9053709
# # HHP-Flare       -0.041564030 -0.15541467 0.07228661 0.7775292
# # Remission-Flare -0.022086501 -0.13131517 0.08714217 0.9525217
# # HHP-HC          -0.014761066 -0.11888080 0.08935867 0.9827436
# # Remission-HC     0.004716463 -0.09432824 0.10376117 0.9993152
# # Remission-HHP    0.019477529 -0.09082467 0.12977973 0.9675878
# 
# # If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# # If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# # Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself
# ## adonis2 will fail if there are samples with rowSums = 0, let's drop them
# #abs.cooccur.final.table[is.na(abs.cooccur.final.table)]<-0
# abs.cooccur.final.table2<-subset(abs.cooccur.final.table[,-1],rowSums(abs.cooccur.final.table[,-1])!=0)
# abs.sub.meta2<-subset(abs.sub.meta,rownames(abs.sub.meta) %in% rownames(abs.cooccur.final.table2))
# # went from 129 --> 97 samples, meaning 32 samples didn't share strains with any other samples
# 
# abs.PA.final.table<-ifelse(abs.cooccur.final.table2[,-1]>0,1,0)
# # use the pvalue below for PCOA that excludes Donor
# pnova1<-adonis2(abs.PA.final.table ~ Condition,data=abs.sub.meta2,method = "jaccard",by="terms",permutations= 10000)
# pnova1 # p-value = 9.999e-05
# p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# # 0.00029997
# 
# ##one issue with adonis is that it doesn't do multiple comparisons *******
# # tells us that something is different, but what is different? Which sample/plot/location?
# ## our four provinces differ, but do all of them differ,or just one?
# 
# ##random person on the internet to the rescue!
# #install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
# 
# abs.cooccur.dist = dist(abs.cooccur.final.table2[,-1],method= "binary") #Jaccard distance matrix
# pair.mod1<-pairwise.adonis(abs.cooccur.dist,abs.sub.meta2$Condition, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
# pair.mod1
# #        pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# # 1        Flare vs HC  1  1.557065 3.966016 0.09681234   1e-04     0.0006  **
# # 2       Flare vs HHP  1  1.254836 2.935974 0.06124782   2e-04     0.0012   *
# # 3 Flare vs Remission  1  1.622496 3.781513 0.07302824   1e-04     0.0006  **
# # 4          HC vs HHP  1  1.500746 3.899418 0.09306616   1e-04     0.0006  **
# # 5    HC vs Remission  1  1.346660 3.453688 0.07769181   1e-04     0.0006  **
# # 6   HHP vs Remission  1  1.240358 2.937042 0.05655005   1e-04     0.0006  **
# 
# # Visualize dispersions
# par(mar=c(1,1,1,1))
# png('figures/SameSTR_Figs/ABS_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
# plot(abs.disper1,main = "Centroids and Dispersion based on Jaccard Distance")
# dev.off()
# 
# par(mar=c(1,1,1,1))
# png('figures/SameSTR_Figs/ABS_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
# boxplot(abs.disper1,xlab="By GroupName", main = "Distance to Centroid by Category", sub="Based on Jaccard Distance", col=grp.clrs$GroupName_Color)
# dev.off()
# 
# 
# #### Visualize PCoA ####
# 
# ## *** all figures that end in _PCOA1 come from the same single PCoA
# #data is just subsetted to understand how points are related to each other w/in & across timepoints
# 
# # create PCoA ggplot fig
# 
# pcoa1<-ggplot(abs.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Condition)) +geom_point(size=4)+theme_bw()+
#   labs(color="Condition",title="PCoA of Presence/Absence of Shared Strains (SameSTR Results)",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
#         axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Condition",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Condition)])) +
#   scale_fill_manual(name ="Condition",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Condition)])) +
#   xlab("PC1 [4.74%]") + ylab("PC2 [2.15%]")+ annotate("text",x=0.40,y=0.5,label="PERMANOVA, p.adj = 0.0003")
# 
# ggsave(pcoa1,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedStrains_Binary_PCOA1.png", 
#        width = 2000,height = 1500,units = "px",
#        dpi = 200,create.dir = TRUE)
# 
# pcoa2<-ggplot(abs.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Condition)) +geom_point(size=4)+theme_bw()+
#   labs(color="Condition",title="PCoA of Presence/Absence of Shared Strains (SameSTR Results)",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
#         axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Condition",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Condition)])) +
#   scale_fill_manual(name ="Condition",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Condition)])) +
#   xlab("PC1 [4.74%]") + ylab("PC2 [2.15%]")+ stat_ellipse(alpha = 0.5) + annotate("text",x=0.80,y=0.55,label="PERMANOVA, p.adj = 0.0003")
# 
# ggsave(pcoa2,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedStrains_Binary_Ellipses_PCOA2.png", 
#        width = 2000,height = 1500,units = "px",
#        dpi = 200,create.dir = TRUE)
# 
# pcoa3<-ggplot(abs.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Condition)) +geom_text(label=abs.jac.pcoa1.meta$SampleID,position="identity")+theme_bw()+
#   labs(color="Condition",title="PCoA of Presence/Absence of Shared Strains (SameSTR Results)",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
#         axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Condition",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Condition)])) +
#   scale_fill_manual(name ="Condition",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Condition)])) +
#   xlab("PC1 [4.74%]") + ylab("PC2 [2.15%]")+ annotate("text",x=0.40,y=0.5,label="PERMANOVA, p.adj = 0.0003")
# 
# ggsave(pcoa3,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedStrains_Binary_Labeled_by_SampleID_PCOA3.png", 
#        width = 2000,height = 1500,units = "px",
#        dpi = 200,create.dir = TRUE)

#### Import FMT Metadata ####

fmt.metadata<-as.data.frame(read_xlsx("/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/data/ABS_FMT_Metadata_CLH.xlsx", sheet="Sheet1"))
fmt.metadata$SampleID<-gsub("_L00[0-9]_R1_001.fastq.gz","",fmt.metadata$filename)
unique(fmt.metadata$SampleID)

sample.order<-unique(fmt.metadata$SampleID)

rownames(fmt.metadata)<-fmt.metadata$SampleID
# 
# # create group color palette to merge with abs.all.basic.meta
# grp.clrs = as.data.frame(t(data.frame("HC"="#40557B","Flare"="#ffa600","Remission"="#ff6361","HHP"="#007561")))
# grp.clrs
# 
# grp.clrs$Condition<-rownames(grp.clrs)
# colnames(grp.clrs)[which(names(grp.clrs) == "V1")] <- "Group_Color"
# grp.clrs

# merge FMT metadata with group colors object
fmt.meta.clrs<-merge(fmt.metadata,grp.clrs,by="Condition")
unique(fmt.meta.clrs$Group_Color)

rownames(fmt.meta.clrs)<-fmt.meta.clrs$SampleID

# subset out only FMT/FMT-HHP/Donor - Cooccurr data
fmt.cooccur.meta<-abs.cooccur.meta[which(abs.cooccur.meta$SampleID %in% fmt.meta.clrs$SampleID & abs.cooccur.meta$col %in% fmt.meta.clrs$SampleID),]
head(fmt.cooccur.meta)
unique(fmt.cooccur.meta$SampleID) %in% unique(fmt.meta.clrs$SampleID) # sanity check, should all be TRUE

# subset out only FMT/FMT-HHP/Donor - Events data
head(abs.events.meta)
fmt.events.meta<-abs.events.meta[which(abs.events.meta$SampleID %in% fmt.meta.clrs$SampleID & abs.events.meta$col %in% fmt.meta.clrs$SampleID),]
head(fmt.events.meta)
unique(fmt.events.meta$SampleID) %in% unique(fmt.meta.clrs$SampleID) # sanity check, should all be TRUE

# merge FMT/FMT-HHP/Donor samestr data with FMT abs.all.basic.meta
fmt.cooccur.2meta<-merge(fmt.cooccur.meta,fmt.meta.clrs[,-c(3:4)],by=c("SampleID","Condition","Patient","Group_Color"))
head(fmt.cooccur.2meta)

# turn fmt.cooccur.2meta$FMT_Label & Condition columns into factors for ordering
unique(fmt.cooccur.2meta$FMT_Label)
fmt.cooccur.2meta$FMT_Label<-factor(fmt.cooccur.2meta$FMT_Label,
                                levels=c("Pre-FMT","Abx-1",
                                         "Post-FMT1","Abx-2","Post-FMT2","Donor"))

unique(fmt.cooccur.2meta$Condition)
#fmt.cooccur.2meta$Condition<-factor(fmt.cooccur.2meta$Condition,levels=c("HC","Flare","Remission","HHP"))

fmt.cooccur.2meta$Patient<-factor(fmt.cooccur.2meta$Patient,levels=c("A020","A021","Donor"))

fmt.cooccur.2meta$SampleID<-factor(fmt.cooccur.2meta$SampleID,levels=c(sample.order))

# order by sampleID and Patient since they are in chronological order
fmt.cooccur.2meta<-fmt.cooccur.2meta[order(fmt.cooccur.2meta$SampleID,fmt.cooccur.2meta$Patient),]
head(fmt.cooccur.2meta)

#### Visualize FMT Shared Strain Data ####

# heatmap faceted by condition
hm.test3<-ggplot(fmt.cooccur.2meta, aes(SampleID, col, fill= shared_strain)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=25,
                       breaks=c(0,10,20,30,40),
                       labels=c("0","10","20","30","40"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="", y="", title="Shared Strains by Sample + Condition - FMT Only",fill="# of Shared Strains") + facet_grid(Condition~., scales = "free",space="free") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(hm.test3,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/SameSTR_Figs/ABS_FMT_MGMs_SameSTR_shared_strains_heatmap.png",width = 6000,height = 5000,units = "px",dpi = 200,create.dir = TRUE)
# 
# hm.test4<-ggplot(fmt.cooccur.2meta, aes(SampleID, col, fill= shared_strain_percent)) +geom_tile(color="white")+
#   scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=0.85,
#                        breaks=c(0,0.5,1,1.5,1.95),
#                        labels=c("0","0.5","1","1.5","1.95"))+
#   theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
#                         axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
#                         legend.title = element_text(hjust=0.5,size=30),
#                         legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
#                         strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
#   labs(x="", y="", title="Shared Strains by Sample + Condition - FMT Only",fill="% of Shared Strains") + facet_grid(Condition~., scales = "free",space="free") +
#   guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()
# 
# ggsave(hm.test4,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/SameSTR_Figs/ABS_FMT_MGMs_SameSTR_shared_strains_percent_heatmap.png",width = 6000,height = 5000,units = "px",dpi = 200,create.dir = TRUE)

#### Visualize FMT Shared Strain Results without HHP A021 ####

head(fmt.cooccur.2meta)
fmt.NoHHP.cooccur.2meta<-fmt.cooccur.2meta[!grepl("A021_",fmt.cooccur.2meta$col) & !grepl("A021",fmt.cooccur.2meta$Patient),]

sample.order.noHHP<-sample.order[!grepl("A021",sample.order)]

fmt.NoHHP.cooccur.2meta$SampleID<-factor(fmt.NoHHP.cooccur.2meta$SampleID,levels=c(sample.order.noHHP))
fmt.NoHHP.cooccur.2meta$col<-factor(fmt.NoHHP.cooccur.2meta$col,levels=c(sample.order.noHHP))

# order by sampleID and Patient since they are in chronological order
fmt.NoHHP.cooccur.2meta<-fmt.NoHHP.cooccur.2meta[order(fmt.NoHHP.cooccur.2meta$SampleID,fmt.NoHHP.cooccur.2meta$col),]

# heatmap faceted by condition
hm.test5<-ggplot(fmt.NoHHP.cooccur.2meta, aes(SampleID, col, fill= shared_strain)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=25,
                       breaks=c(0,10,20,30,40),
                       labels=c("0","10","20","30","40"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="", y="", title="Shared Strains by Sample + Condition - FMT Only, NO HHP",fill="# of Shared Strains") + facet_grid(Condition~., scales = "free",space="free") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand = c(0,0))+coord_flip()

ggsave(hm.test5,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/SameSTR_Figs/ABS_FMT_NoHHP_MGMs_SameSTR_shared_strains_heatmap.png",width = 6000,height = 5000,units = "px",dpi = 200,create.dir = TRUE)

# #### Generate FMT PCoA with SameSTR Results ####
# head(fmt.cooccur.2meta)
# fmt.cooccur.2meta$shared_strain<-as.numeric(fmt.cooccur.2meta$shared_strain)
# fmt.cooccur.table<-dcast(fmt.cooccur.2meta,col~SampleID,value.var = "shared_strain")
# fmt.cooccur.table[1:5,1:5] # sanity check
# rownames(fmt.cooccur.table)<-fmt.cooccur.table$col
# 
# # To continue, we have to make a symmetric table...
# # problem: fmt.cooccur.table is not symmetric, causing us to lose data points in PCOA...
# ## code to the rescue from here: https://stackoverflow.com/questions/12591708/from-asymmetric-matrix-or-dataframe-into-a-symmetric-square-matrix-with-r
# 
# # use newly created cooccurr table to get the names as a vector
# fmt.NAMES <- sort(unique(c(rownames(fmt.cooccur.table[,-1]), colnames(fmt.cooccur.table[,-1]))))
# # create temp data frame using stack() to keep the rows we would lose with dcast to make table
# fmt.temp <- data.frame(rows = rownames(fmt.cooccur.table[,-1]), stack(fmt.cooccur.table[,-1]))
# 
# # factor the row and column names using NAMES vector
# fmt.temp$rows <- factor(fmt.temp$rows, fmt.NAMES)
# fmt.temp$ind <- factor(fmt.temp$ind, fmt.NAMES)
# 
# # Use xtabs to get your desired output of symmetric table. Wrap it in as.data.frame.matrix to get a data.frame as output
# fmt.cooccur.table.sym<-as.data.frame.matrix(xtabs(values ~ rows + ind, fmt.temp))
# fmt.cooccur.table.sym[1:5,1:5]
# 
# # drop rows where all elements are NA
# ## rowSums(is.na()) counts # of NAs, then !=ncol() sees if the # of NAs is less than the total # of columns
# # fmt.cooccur.table2<-fmt.cooccur.table[,-1][rowSums(is.na(fmt.cooccur.table[,-1])) != ncol(fmt.cooccur.table[,-1]),]
# # fmt.cooccur.table[is.na(fmt.cooccur.table)]<-0
# 
# # PCOA w/ Jaccard distance matrix (of binary data)
# #pcoa(vegdist(fmt.cooccur.table.sym,na.rm = TRUE,binary=TRUE,method="jaccard"),correction="lingoes") # 
# fmt.jac.pcoa <- pcoa(dist(fmt.cooccur.table.sym,method="binary"),correction="lingoes") # 
# #save.image("data/Amplicon/FT_16SV3V4_CLR_EucDist_Ready.Rdata")
# 
# # The proportion of variances explained is in its element values$Relative_eig
# fmt.jac.pcoa$values
# 
# # extract principal coordinates
# fmt.jac.pcoa.vectors<-data.frame(fmt.jac.pcoa$vectors)
# fmt.jac.pcoa.vectors$SampleID<-rownames(fmt.jac.pcoa$vectors)
# 
# # merge pcoa coordinates w/ fmt.meta.clrs
# fmt.jac.pcoa.meta<-merge(fmt.jac.pcoa.vectors, fmt.meta.clrs, by.x="SampleID", by.y="SampleID")
# 
# head(fmt.jac.pcoa.meta)
# rownames(fmt.jac.pcoa.meta)<-fmt.jac.pcoa.meta$SampleID
# 
# head(fmt.jac.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# # PC1 = 19.5%, PC2 = 11.69%
# 
# #### Homogeneity of Variance & PERMANOVA tests - FMT Only ####
# ## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# # multivariate analogue to Levene's test of homogeneity of variances
# # program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA
# 
# #While PERMANOVA tests differences in group means (analogous to MANOVA),
# ## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
# #(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
# ## * need a distance matrix!
# 
# head(fmt.meta.clrs)
# 
# # are both dfs in the same order? Need this for stats below
# rownames(fmt.meta.clrs) %in% rownames(fmt.cooccur.table.sym)
# 
# # first by compare dispersions by sampling date
# fmt.disper1<-betadisper((dist(fmt.cooccur.table.sym,method="binary")), fmt.meta.clrs$Condition)
# fmt.disper1
# 
# ## Significant differences in homogeneities can be tested using either parametric or permutational tests,
# ##and parametric post hoc contrasts can also be investigated:
# 
# permutest(fmt.disper1, pairwise=TRUE,permutations=9999) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
# #Pairwise comparisons:
# #  (Observed p-value below diagonal, permuted p-value above diagonal)
# #             Flare HC     HHP Remission
# # Flare                0.95720    0.3894
# # HC                                    
# # HHP       0.95718               0.4485
# # Remission 0.38827    0.45488          
# 
# anova(fmt.disper1) # p = 0.501 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across conditions
# # ANOVA adjusted p-value
# aov.beta.p1<-anova(fmt.disper1)[["Pr(>F)"]] # get p values from ANOVA
# p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))
# 
# TukeyHSD(fmt.disper1) # tells us which category's dispersion MEANS are significantly different than each other
# #                     diff       lwr       upr     p adj
# # HC-Flare        -0.44751504 -1.3550392 0.4600091 0.5352657
# # HHP-Flare       -0.01044016 -0.4642022 0.4433219 0.9999055
# # Remission-Flare -0.11690094 -0.5148766 0.2810747 0.8489651
# # HHP-HC           0.43707488 -0.4704493 1.3445990 0.5544421
# # Remission-HC     0.33061410 -0.5503426 1.2115708 0.7307638
# # Remission-HHP   -0.10646078 -0.5044364 0.2915149 0.8808193
# 
# # If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# # If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# # Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself
# 
# ## adonis2 will fail if there are samples with rowSums = 0, let's drop them
# fmt.cooccur.table.sym[is.na(fmt.cooccur.table.sym)]<-0
# fmt.cooccur.table.sym2<-subset(fmt.cooccur.table.sym,rowSums(fmt.cooccur.table.sym)!=0)
# fmt.meta.clrs2<-subset(fmt.meta.clrs,rownames(fmt.meta.clrs) %in% rownames(fmt.cooccur.table.sym2))
# # went from 28 --> 24 samples, so 4 samples had rowSums = 0 aka no shared strains with other samples identified by SameStr
# 
# # use the pvalue below for PCOA 
# ## not enough observations to have strata="Patient"
# pnova2<-adonis2(fmt.cooccur.table.sym2 ~ Condition,data=fmt.meta.clrs2,method = "jaccard",by="terms",permutations= 10000)
# pnova2 # p-value = 0.002
# p.adjust(pnova2$`Pr(>F)`,method="bonferroni",n=length(pnova2$`Pr(>F)`)) # adjusted pval
# # 0.0059994
# 
# ##one issue with adonis is that it doesn't do multiple comparisons *******
# # tells us that something is different, but what is different? Which sample/plot/location?
# ## our four provinces differ, but do all of them differ,or just one?
# 
# ##random person on the internet to the rescue!
# #install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
# 
# fmt.sym2.dist = dist(fmt.cooccur.table.sym2,method= "binary") #Jaccard distance matrix
# pair.mod2<-pairwise.adonis(fmt.sym2.dist,fmt.meta.clrs2$Condition, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
# pair.mod2
# #                 pairs Df SumsOfSqs   F.Model      R2  p.value p.adjusted sig
# # 1        Flare vs HC  1 0.1773148 0.6242869 0.1350018  0.5000     1.0000    
# # 2       Flare vs HHP  1 0.8833462 3.3603557 0.2515169  0.0208     0.1248    
# # 3 Flare vs Remission  1 1.5369744 8.2130002 0.3697384  0.0005     0.0030   *
# # 4          HC vs HHP  1 0.3566569 1.4336839 0.1928632  0.3793     1.0000    
# # 5    HC vs Remission  1 0.3695513 2.4905098 0.1993922  0.0835     0.5010    
# # 6   HHP vs Remission  1 1.4592047 7.8439871 0.3289713  0.0007     0.0042   *
# 
# # Visualize dispersions
# ## order of colors comes from looking at fmt.disper1 order; HC has distance to median as 0 - helpful indicator here
# par(mar=c(1,1,1,1))
# png('figures/SameSTR_Figs/ABS_FMT_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
# plot(fmt.disper1,main = "Centroids and Dispersion based on Jaccard Distance", col=c("#ffa600","#40557B","#007561","#ff6361"))
# dev.off()
# 
# par(mar=c(1,1,1,1))
# png('figures/SameSTR_Figs/ABS_FMT_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
# boxplot(fmt.disper1,xlab="By GroupName", main = "Distance to Centroid by Category", sub="Based on Jaccard Distance", col=c("#ffa600","#40557B","#007561","#ff6361"))
# dev.off()
# 
# #### Visualize FMT-Only PCoA ####
# 
# ## *** all figures that end in _PCOA1 come from the same single PCoA
# #data is just subsetted to understand how points are related to each other w/in & across timepoints
# 
# # create PCoA ggplot fig
# # PC1 = 19.5%, PC2 = 11.69%
# 
# pcoa<-ggplot(fmt.jac.pcoa.meta, aes(x=Axis.1, y=Axis.2,color=Condition,shape=Patient)) +geom_point(size=4)+theme_bw()+
#   labs(color="Condition",title="PCoA of Presence/Absence of Shared Strains in FMT Samples (SameSTR Results)",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
#         axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Condition",values=unique(fmt.jac.pcoa.meta$Group_Color[order(fmt.jac.pcoa.meta$Condition)])) +
#   scale_fill_manual(name ="Condition",values=unique(fmt.jac.pcoa.meta$Group_Color[order(fmt.jac.pcoa.meta$Condition)])) +
#   xlab("PC1 [19.5%]") + ylab("PC2 [11.69%]")+ annotate("text",x=0.5,y=0.8,label="PERMANOVA, p.adj = 0.006")
# 
# ggsave(pcoa,filename = "figures/SameSTR_Figs/ABS_FMT_MGMs_SameSTR_SharedStrains_Binary_PCOA1.png", 
#        width = 2000,height = 1500,units = "px",
#        dpi = 200,create.dir = TRUE)
# 
# pcoa2<-ggplot(fmt.jac.pcoa.meta, aes(x=Axis.1, y=Axis.2,color=Condition)) +geom_point(size=4)+theme_bw()+
#   labs(color="Condition",title="PCoA of Presence/Absence of Shared Strains in FMT Samles (SameSTR Results)",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
#         axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Condition",values=unique(fmt.jac.pcoa.meta$Group_Color[order(fmt.jac.pcoa.meta$Condition)])) +
#   scale_fill_manual(name ="Condition",values=unique(fmt.jac.pcoa.meta$Group_Color[order(fmt.jac.pcoa.meta$Condition)])) +
#   xlab("PC1 [19.5%]") + ylab("PC2 [11.69%]")+ stat_ellipse(alpha = 0.5) + annotate("text",x=1.1,y=0.8,label="PERMANOVA, p.adj = 0.006")
# 
# ggsave(pcoa2,filename = "figures/SameSTR_Figs/ABS_FMT_MGMs_SameSTR_SharedStrains_Binary_Ellipses_PCOA2.png", 
#        width = 2000,height = 1500,units = "px",
#        dpi = 200,create.dir = TRUE)
# 
# pcoa3<-ggplot(fmt.jac.pcoa.meta, aes(x=Axis.1, y=Axis.2,color=Condition)) +geom_text(label=fmt.jac.pcoa.meta$SampleID,position="identity")+theme_bw()+
#   labs(color="Condition",title="PCoA of Presence/Absence of Shared Strains (SameSTR Results)",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
#         axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Condition",values=unique(fmt.jac.pcoa.meta$Group_Color[order(fmt.jac.pcoa.meta$Condition)])) +
#   scale_fill_manual(name ="Condition",values=unique(fmt.jac.pcoa.meta$Group_Color[order(fmt.jac.pcoa.meta$Condition)])) +
#   xlab("PC1 [19.5%]") + ylab("PC2 [11.69%]")+ annotate("text",x=0.5,y=0.8,label="PERMANOVA, p.adj = 0.006")
# 
# ggsave(pcoa3,filename = "figures/SameSTR_Figs/ABS_FMT_MGMs_SameSTR_SharedStrains_Binary_Labeled_by_SampleID_PCOA3.png", 
#        width = 2000,height = 1500,units = "px",
#        dpi = 200,create.dir = TRUE)

#### A020 vs FMT - Percentage of Shared Strains ####

# how many total shared strains are there across all sample comparisons?
length(unique(abs.samestr.events$clade)) # 293 clades considered
a020.clades<-abs.samestr.events[which(abs.samestr.events$SampleID %in% "PA_09_06_22_S71" | abs.samestr.events$col %in% "PA_09_06_22_S71"),]
length(unique(a020.clades$clade)) # 36 total clades in PA_09_06_22_S71
length(unique(a020.clades$clade[a020.clades$event=="shared_strain"])) # 25 total clades with shared_strain events in PA_09_06_22_S71

# sample original a020 had 36 clades identified -- denominator for FMT barplot panel 1

a020.patterns<-c("A020","PA")

a020.fmt.only<-subset(fmt.cooccur.2meta,Patient=="A020" & grepl(paste(a020.patterns, collapse='|'),fmt.cooccur.2meta$col))

a020.fmt.table<-dcast(a020.fmt.only,SampleID~col,value.var="shared_strain")
rownames(a020.fmt.table)<-a020.fmt.table$SampleID

a020.fmt1<-a020.fmt.table[which(colnames(a020.fmt.table)=="PA_09_06_22_S71")]
a020.fmt1<-subset(a020.fmt1,rownames(a020.fmt1)!="PA_09_06_22_S71")
a020.fmt1[is.na(a020.fmt1)]<-0
colnames(a020.fmt1)[which(colnames(a020.fmt1)=="PA_09_06_22_S71")]<-"A020.Count1"
head(a020.fmt1)

a020.fmt2<-as.data.frame(t(a020.fmt.table[a020.fmt.table$SampleID=="PA_09_06_22_S71",-1]))
a020.fmt2<-subset(a020.fmt2,rownames(a020.fmt2)!="PA_09_06_22_S71")
a020.fmt2[is.na(a020.fmt2)]<-0
colnames(a020.fmt2)[which(colnames(a020.fmt2)=="PA_09_06_22_S71")]<-"A020.Count2"
head(a020.fmt2)

rownames(a020.fmt1) %in% rownames(a020.fmt2)
rownames(a020.fmt1)[which(!rownames(a020.fmt1) %in% rownames(a020.fmt2))]

rownames(a020.fmt2) %in% rownames(a020.fmt1)
rownames(a020.fmt2)[which(!rownames(a020.fmt2) %in% rownames(a020.fmt1))]

a020.fmt1[nrow(a020.fmt1)+1,]<-0 # adding counts that are missing in a020.fmt2 but in a020.fmt2
rownames(a020.fmt1)[which(rownames(a020.fmt1)=="19")]<-"S_A020_28_S25"
a020.fmt1$SampleID<-rownames(a020.fmt1)

a020.fmt2[nrow(a020.fmt2)+1,]<-0 # adding counts that are missing in a020.fmt2 but in a020.fmt2
rownames(a020.fmt2)[which(rownames(a020.fmt2)=="19")]<-"A020_12_S8"
a020.fmt2$SampleID<-rownames(a020.fmt2)

a020.fmt.counts<-merge(a020.fmt1,a020.fmt2,by="SampleID")
a020.fmt.counts
a020.fmt.counts$Total.vs.A020<-a020.fmt.counts$A020.Count1+a020.fmt.counts$A020.Count2
a020.fmt.counts$SharedPercent<-as.numeric((a020.fmt.counts$Total.vs.A020/36)*100) #shared totals with self / total # clades found in PA_09_06_22_S71

# add baseline sample back to this df for plot
a020.fmt.counts[nrow(a020.fmt.counts) + 1,]<-c("PA_09_06_22_S71",0,0,36,100) #36 total clades aka 100%
a020.fmt.counts$SharedPercent<-as.numeric(a020.fmt.counts$SharedPercent)
#rownames(a020.fmt.counts)[which(rownames(a020.fmt.counts)=="20")]<-"PA_09_06_22_S71"
a020.fmt.counts # sanity check that this worked

# drop samples that are pre-FMT ("A020_8_S10","A020_12_S8","A020_13_S9")
a020.fmt.counts<-a020.fmt.counts[which(!a020.fmt.counts$SampleID %in% c("A020_8_S10","A020_12_S8","A020_13_S9")),]

# merge fmt.baslines with fmt.metadata that has color labels (fmt.meta.clrs)
a020.compare.meta<-merge(a020.fmt.counts[,c(1,4:5)],fmt.meta.clrs,by="SampleID")
head(a020.compare.meta)
rownames(a020.compare.meta)<-a020.compare.meta$SampleID

# arrange samples in chronological order
a020.order<-data.frame(SampleID=sample.order[sample.order %in% a020.compare.meta$SampleID])
a020.order$SampleID<-factor(a020.order$SampleID,levels=c(a020.order$SampleID))
a020.order$Timeline<-c("T-60","T-8","T-1","T+5","T+12","T+26","T+185","T+188",
                        "T+206","T+216","T+274","T+277","T+286","T+314","T+375","T+376","T+410")

a020.compare.meta2<-merge(a020.compare.meta,a020.order,by="SampleID")
a020.compare.meta2$SampleID<-factor(a020.compare.meta2$SampleID,levels=c(a020.order$SampleID))
a020.compare.meta2$Timeline<-factor(a020.compare.meta2$Timeline,levels=c(a020.order$Timeline))
a020.compare.meta2$FMT_Label<-factor(a020.compare.meta2$FMT_Label,
                                      levels=c("Pre-FMT","Abx-1",
                                               "Post-FMT1","Abx-2","Post-FMT2","A020"))

a020.baseline.plot1<-ggplot(a020.compare.meta2, aes(x = Timeline, y = SharedPercent,fill=Condition)) +
  geom_col(position = "stack",color="black",width=0.8)+
  theme_classic()+
  scale_y_continuous(name = "% shared with pre-FMT",labels = function(x) paste0(x*1, "%"),limits=c(0,100))+
  scale_fill_manual(name = "Flare",values=unique(a020.compare.meta$Group_Color[order(a020.compare.meta$Condition)])) +
  theme(axis.title.x = element_blank(),
        text = element_text(colour = "black", size=8),
        axis.text=element_text(colour="black", size=8),
        axis.text.x = element_text(angle=90, vjust=.93, hjust=1.01,size=8),
        strip.background=element_rect(colour="white"),
        strip.text.x= element_text(colour = "black", size=9)) +
  labs(title = "") + facet_grid(.~FMT_Label, scales = "free",space="free")

ggsave(a020.baseline.plot1,filename = "figures/SameSTR_Figs/ABS_A020_FMT_Shared_Strains_Percentage_Plot1_barplot.png", width = 1800,height = 400,units = "px",dpi = 200,create.dir = TRUE)

a020.baseline.plot2<-ggplot(a020.compare.meta2, aes(x = Timeline, y = SharedPercent,fill=Condition)) +
  geom_col(position = "stack",color="black",width=0.8)+
  theme_classic()+
  scale_y_continuous(name = "% shared with pre-FMT",labels = function(x) paste0(x*1, "%"),limits=c(0,100))+
  scale_fill_manual(name = "Flare",values=unique(a020.compare.meta$Group_Color[order(a020.compare.meta$Condition)])) +
  theme(axis.title.x = element_blank(),
        text = element_text(colour = "black", size=10),
        axis.text=element_text(colour="black", size=10),
        axis.text.x = element_blank(),
        strip.background=element_rect(colour="white"),
        strip.text.x= element_text(colour = "black", size=10),
        legend.text = element_text(size=10)) +
  labs(title = "") + facet_grid(.~FMT_Label, scales = "free",space="free")

ggsave(a020.baseline.plot2,filename = "figures/SameSTR_Figs/ABS_A020_FMT_Shared_Strains_Percentage_Plot1_barplot2.png", width = 1800,height = 400,units = "px",dpi = 200,create.dir = TRUE)

#### Donor vs FMT - Percentage of Shared Strains ####

# how many total shared strains are there across all sample comparisons?
length(unique(abs.samestr.events$clade)) # 293 clades considered
length(unique(abs.samestr.events$clade[abs.samestr.events$SampleID=="Donor41_33_S77"]))
donor.clades<-abs.events.meta[which(abs.events.meta$SampleID %in% "Donor41_33_S77" | abs.events.meta$col %in% "Donor41_33_S77"),]
length(unique(donor.clades$clade)) # 31 total clades in Donor
length(unique(donor.clades$clade[donor.clades$event=="shared_strain"])) # 19 clades with shared strain events Donor

# sample Donor had 31 clades identified -- denominator for FMT barplot panel 1

donor.patterns<-c("A020","PA","Donor")
patient.patterns<-c("A020","Donor")

donor.fmt.only<-subset(fmt.cooccur.2meta,grepl(paste(patient.patterns, collapse='|'),fmt.cooccur.2meta$Patient) & grepl(paste(donor.patterns, collapse='|'),fmt.cooccur.2meta$col))

donor.fmt.table<-dcast(donor.fmt.only,SampleID~col,value.var="shared_strain")
rownames(donor.fmt.table)<-donor.fmt.table$SampleID

donor.fmt1<-donor.fmt.table[which(colnames(donor.fmt.table)=="Donor41_33_S77")]
donor.fmt1<-subset(donor.fmt1,rownames(donor.fmt1)!="Donor41_33_S77")
donor.fmt1[is.na(donor.fmt1)]<-0
colnames(donor.fmt1)[which(colnames(donor.fmt1)=="Donor41_33_S77")]<-"Donor.Count1"
head(donor.fmt1)

donor.fmt2<-as.data.frame(t(donor.fmt.table[donor.fmt.table$SampleID=="Donor41_33_S77",-1]))
donor.fmt2<-subset(donor.fmt2,rownames(donor.fmt2)!="Donor41_33_S77")
donor.fmt2[is.na(donor.fmt2)]<-0
colnames(donor.fmt2)[which(colnames(donor.fmt2)=="Donor41_33_S77")]<-"Donor.Count2"
head(donor.fmt2)

rownames(donor.fmt1) %in% rownames(donor.fmt2)
rownames(donor.fmt1)[which(!rownames(donor.fmt1) %in% rownames(donor.fmt2))]

rownames(donor.fmt2) %in% rownames(donor.fmt1)
rownames(donor.fmt2)[which(!rownames(donor.fmt2) %in% rownames(donor.fmt1))]

donor.fmt1[nrow(donor.fmt1)+1,]<-0 # adding counts that are missing in donor.fmt2 but in donor.fmt2
rownames(donor.fmt1)[which(rownames(donor.fmt1)=="20")]<-"S_A020_28_S25"
donor.fmt1$SampleID<-rownames(donor.fmt1)

donor.fmt2[nrow(donor.fmt2)+1,]<-0 # adding counts that are missing in donor.fmt2 but in donor.fmt2
rownames(donor.fmt2)[which(rownames(donor.fmt2)=="20")]<-"A020_12_S8"
donor.fmt2$SampleID<-rownames(donor.fmt2)

donor.fmt.counts<-merge(donor.fmt1,donor.fmt2,by="SampleID")
donor.fmt.counts$Total.vs.Donor<-donor.fmt.counts$Donor.Count1+donor.fmt.counts$Donor.Count2
donor.fmt.counts$SharedPercent<-(donor.fmt.counts$Total.vs.Donor/31)*100 #shared totals with self / total # clades found in PA_09_06_22_S71

# add baseline sample back to this df for plot
#donor.baselines[nrow(donor.baselines) + 1,]<-c(28,0,28,100.0) #36 total clades aka 100%
#rownames(donor.baselines)[which(rownames(donor.baselines)=="21")]<-"Donor41_33_S77"
#donor.baselines # sanity check that this worked

# drop samples that are pre-FMT ("A020_8_S10","A020_12_S8","A020_13_S9")
donor.fmt.counts<-donor.fmt.counts[which(!donor.fmt.counts$SampleID %in% c("A020_8_S10","A020_12_S8","A020_13_S9")),]

# merge fmt.baslines with fmt.metadata that has color labels (fmt.meta.clrs)
donor.compare.meta<-merge(donor.fmt.counts[,c(1,4:5)],fmt.meta.clrs,by="SampleID")
head(donor.compare.meta)
rownames(donor.compare.meta)<-donor.compare.meta$SampleID

# arrange samples in chronological order
donor.order<-data.frame(SampleID=sample.order[sample.order %in% donor.compare.meta$SampleID])
donor.order$SampleID<-factor(donor.order$SampleID,levels=c(donor.order$SampleID))
donor.order$Timeline<-c("T-60","T-8","T-1","T+5","T+12","T+26","T+185","T+188",
                       "T+206","T+216","T+274","T+277","T+286","T+314","T+375","T+376","T+410")

donor.compare.meta2<-merge(donor.compare.meta,donor.order,by="SampleID")
donor.compare.meta2$SampleID<-factor(donor.compare.meta2$SampleID,levels=c(donor.order$SampleID))
donor.compare.meta2$Timeline<-factor(donor.compare.meta2$Timeline,levels=c(donor.order$Timeline))
donor.compare.meta2$FMT_Label<-factor(donor.compare.meta2$FMT_Label,
                                      levels=c("Pre-FMT","Abx-1",
                                               "Post-FMT1","Abx-2","Post-FMT2","Donor"))

# make the bar plots
donor.baseline.plot1<-ggplot(donor.compare.meta2, aes(x = Timeline, y = SharedPercent,fill=Condition)) +
  geom_col(position = "stack",color="black",width=0.8)+
  theme_classic()+
  scale_y_continuous(name = "% shared with Donor",labels = function(x) paste0(x*1, "%"),limits=c(0,100))+
  scale_fill_manual(name = "Flare",values=unique(donor.compare.meta$Group_Color[order(donor.compare.meta$Condition)])) +
  theme(axis.title.x = element_blank(),
        text = element_text(colour = "black", size=8),
        axis.text=element_text(colour="black", size=8),
        axis.text.x = element_text(angle=90, vjust=.93, hjust=1.01),
        strip.background=element_rect(colour="white"),
        strip.text.x= element_text(colour = "black", size=9)) +
  labs(title = "") + facet_grid(.~FMT_Label, scales = "free",space="free")

ggsave(donor.baseline.plot1,filename = "figures/SameSTR_Figs/ABS_Donor_FMT_Shared_Strains_Percentage_Plot2_barplot.png", width = 1800,height = 400,units = "px",dpi = 200,create.dir = TRUE)

donor.baseline.plot2<-ggplot(donor.compare.meta2, aes(x = Timeline, y = SharedPercent,fill=Condition)) +
  geom_col(position = "stack",color="black",width=0.8)+
  theme_classic()+
  scale_y_continuous(name = "% shared with Donor",labels = function(x) paste0(x*1, "%"),limits=c(0,100))+
  scale_fill_manual(name = "Flare",values=unique(donor.compare.meta$Group_Color[order(donor.compare.meta$Condition)])) +
  theme(axis.title.x = element_blank(),
        text = element_text(colour = "black", size=10),
        axis.text=element_text(colour="black", size=10),
        axis.text.x = element_blank(),
        strip.background=element_rect(colour="white"),
        strip.text.x = element_text(colour="white",size=10),
        legend.text = element_text(size=10)) +
  labs(title = "") + facet_grid(.~FMT_Label, scales = "free",space="free")

ggsave(donor.baseline.plot2,filename = "figures/SameSTR_Figs/ABS_Donor_FMT_Shared_Strains_Percentage_Plot2_barplot2.png", width = 1800,height = 400,units = "px",dpi = 200,create.dir = TRUE)

#### Engraftment - Percentage of Shared Strains ####
# NOTE
## this is the number of shared strains between post-FMT and the donor excluding the strains shared between pre-FMT and the donor samples. 
## So if a strain is shared between post-FMT and donor, but already existed pre-FMT, we would not call that engraftment, as it already existed
head(fmt.events.meta)

unique(a020.clades$clade)
unique(c(a020.clades$SampleID,a020.clades$col))
unique(donor.clades$clade)
unique(c(donor.clades$SampleID,donor.clades$col))

donor.patterns

pre.fmt.samps<-c("PA_09_06_22_S71","PA_10_28_22_S72","PA_11_04_22_S73")

pre.fmt.clades<-subset(abs.events.meta,grepl(paste(pre.fmt.samps, collapse='|'),abs.events.meta$SampleID) | grepl(paste(pre.fmt.samps, collapse='|'),abs.events.meta$col))
#pre.fmt.clades2<-subset(pre.fmt.clades1,grepl(paste(donor.patterns, collapse='|'),pre.fmt.clades1$SampleID) | grepl(paste(donor.patterns, collapse='|'),pre.fmt.clades1$col))

unique(pre.fmt.clades$clade)
unique(pre.fmt.clades$clade)[which(unique(pre.fmt.clades$clade) %in% unique(donor.clades$clade))] # check this way
unique(pre.fmt.clades$clade[which(pre.fmt.clades$clade %in% donor.clades$clade)]) # another sanity check

head(donor.clades)

# first drop non-FMT-clades and the pre-FMT samples

# subset all strain events by comparisons with Donor
engraftment.clades1<-donor.clades[which(!donor.clades$clade %in% pre.fmt.clades$clade),]

# subset all Donor-compared strain events with events just for FMT patient
engraftment.clades2<-subset(engraftment.clades1,grepl(paste(donor.patterns, collapse='|'),engraftment.clades1$col) & 
                              grepl(paste(donor.patterns, collapse='|'),engraftment.clades1$Patient))

# subset other events that were not specifically shared strain events
engraftment.clades3<-subset(engraftment.clades2,engraftment.clades2$event=="shared_strain")

# sanity checks
unique(engraftment.clades3$clade) %in% unique(pre.fmt.clades$clade) # should be FALSE

# counting shared strain events between sample pairs, excluding pre-FMT clades
engraftment.table<-dcast(engraftment.clades3,SampleID~col,value.var="event")
rownames(engraftment.table)<-engraftment.table$SampleID

engraft1<-engraftment.table[which(colnames(engraftment.table)=="Donor41_33_S77")]
engraft1<-subset(engraft1,rownames(engraft1)!="Donor41_33_S77")
#engraft1[is.na(engraft1)]<-0
engraft1

engraft2<-as.data.frame(t(engraftment.table[engraftment.table$SampleID=="Donor41_33_S77",-1]))
engraft2<-subset(engraft2,rownames(engraft2)!="Donor41_33_S77")
#engraft2[is.na(engraft2)]<-0
engraft2

rownames(engraft1) %in% rownames(engraft2)
engraft.counts<-rbind(engraft1,engraft2) # bind the engraftment counts together

engraft.counts[nrow(engraft.counts)+4,]<-0 # adding back the pre-FMT samples for plot purposes ("PA_09_06_22_S71","PA_10_28_22_S72","PA_11_04_22_S73")
engraft.counts[is.na(engraft.counts)]<-0
engraft.counts

# add the pre-FMT samples + Abx sample with no shared strains back into df for plotting
rownames(engraft.counts)[which(rownames(engraft.counts)=="14")]<-"PA_09_06_22_S71"
rownames(engraft.counts)[which(rownames(engraft.counts)=="15")]<-"PA_10_28_22_S72"
rownames(engraft.counts)[which(rownames(engraft.counts)=="16")]<-"PA_11_04_22_S73"
rownames(engraft.counts)[which(rownames(engraft.counts)=="17")]<-"A020_2023_08_09_S296"

engraft.counts$SharedPercent<-(engraft.counts$Donor41_33_S77/31)*100 #shared totals with self / total # clades found in PA_09_06_22_S71
engraft.counts$SampleID<-rownames(engraft.counts)
engraft.counts

# drop samples that are pre-FMT ("A020_8_S10","A020_12_S8","A020_13_S9")
# engraft.counts<-engraft.counts[which(!engraft.counts$SampleID %in% c("A020_8_S10","A020_12_S8","A020_13_S9")),]

# merge fmt.baslines with fmt.metadata that has color labels (fmt.meta.clrs)
engraft.compare.meta<-merge(engraft.counts,fmt.meta.clrs,by="SampleID")
head(engraft.compare.meta)
rownames(engraft.compare.meta)<-engraft.compare.meta$SampleID

# arrange samples in chronological order
engraft.order<-data.frame(SampleID=sample.order[sample.order %in% engraft.compare.meta$SampleID])
engraft.order$SampleID<-factor(engraft.order$SampleID,levels=c(engraft.order$SampleID))
engraft.order$Timeline<-c("T-60","T-8","T-1","T+5","T+12","T+26","T+185","T+188",
                        "T+206","T+216","T+274","T+277","T+286","T+314","T+375","T+376","T+410")

engraft.compare.meta2<-merge(engraft.compare.meta,engraft.order,by="SampleID")
engraft.compare.meta2$SampleID<-factor(engraft.compare.meta2$SampleID,levels=c(engraft.order$SampleID))
engraft.compare.meta2$Timeline<-factor(engraft.compare.meta2$Timeline,levels=c(engraft.order$Timeline))
engraft.compare.meta2$FMT_Label<-factor(engraft.compare.meta2$FMT_Label,
                                      levels=c("Pre-FMT","Abx-1",
                                               "Post-FMT1","Abx-2","Post-FMT2","Donor"))

engraft.plot1<-ggplot(engraft.compare.meta2, aes(x = Timeline, y = SharedPercent,fill=Condition)) +
  geom_col(position = "stack",color="black",width=0.8)+
  theme_classic()+
  scale_y_continuous(name = "% engraftment",labels = function(x) paste0(x*1, "%"),limits=c(0,100))+
  scale_fill_manual(name = "Flare",values=unique(engraft.compare.meta$Group_Color[order(engraft.compare.meta$Condition)])) +
  theme(axis.title.x = element_blank(),
        text = element_text(colour = "black", size=8),
        axis.text=element_text(colour="black", size=8),
        axis.text.x = element_text(angle=90, vjust=.93, hjust=1.01),
        strip.background=element_rect(colour="white"),
        strip.text.x= element_text(colour = "black", size=9)) +
  labs(title = "") + facet_grid(.~FMT_Label, scales = "free",space="free")

ggsave(engraft.plot1,filename = "figures/SameSTR_Figs/ABS_Engraftment_FMT_Shared_Strains_Percentage_Plot3_barplot.png", width = 1800,height = 400,units = "px",dpi = 200,create.dir = TRUE)

engraft.plot2<-ggplot(engraft.compare.meta2, aes(x = Timeline, y = SharedPercent,fill=Condition)) +
  geom_col(position = "stack",color="black",width=0.8)+
  theme_classic()+
  scale_y_continuous(name = "% engraftment",labels = function(x) paste0(x*1, "%"),limits=c(0,100))+
  scale_fill_manual(name = "Flare",values=unique(engraft.compare.meta$Group_Color[order(engraft.compare.meta$Condition)])) +
  theme(axis.title.x = element_blank(),
        text = element_text(colour = "black", size=10),
        axis.text=element_text(colour="black", size=10),
        axis.text.x = element_text(angle=90, vjust=.93, hjust=1.01,size=10),
        strip.background=element_rect(colour="white"),
        strip.text.x = element_blank(),
        legend.text = element_text(size=10)) +
  labs(title = "") + facet_grid(.~FMT_Label, scales = "free",space="free")

ggsave(engraft.plot2,filename = "figures/SameSTR_Figs/ABS_Engraftment_FMT_Shared_Strains_Percentage_Plot3_barplot2.png", width = 1800,height = 400,units = "px",dpi = 200,create.dir = TRUE)

#### Combine FMT Plots into 1 Plot ####

fmt.bplot.combo<-ggarrange(a020.baseline.plot2,donor.baseline.plot2,engraft.plot2,ncol=1,nrow=3,legend="right",common.legend = TRUE)

# Annotate the figure by adding a common labels

ggsave(fmt.bplot.combo,filename = "figures/SameSTR_Figs/ABS_Engraftment_FMT_Shared_Strains_Percentage_All_Plots.png", width = 1800,height = 1200,units = "px",dpi = 200,create.dir = TRUE)

sessionInfo() # see loaded packages

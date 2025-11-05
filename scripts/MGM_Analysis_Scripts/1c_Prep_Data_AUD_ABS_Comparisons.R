#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory

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
  library(nationalparkcolors)
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
})

# DATA NOTES
## ABS patients - Flare is an ABS flare up, remission is when that patient is in remission, HHP is healthy household partner

#### Import ABS Data from HuManN - Pathways ####

# first import MetaCyc Pathway data of ABS patients
abs.paths <- as.data.frame(read.delim('data/Metagenomes/SeqProcessing/Humann3_Results/output_merged/output_merged_relab/ABS_all_pathabundance_relab.tsv', sep="\t"))
# ^ contains relative abundance of pathways found with raw, mgm reads using BioBakery 3 suite (mainly HuManN)
# let's clean this file up just a little bit..
colnames(abs.paths)[which(colnames(abs.paths) == "X..Pathway")]<-gsub("X..","", colnames(abs.paths[1])) # clean up Pathway col name
# colnames(abs.paths)[which(colnames(abs.paths) == "X.Pathway")]<-gsub("X.","", colnames(abs.paths)[which(colnames(abs.paths) == "X.Pathway")])
colnames(abs.paths)<-gsub("(.*)_S.*_L00.*_Abundance.RELAB","\\1", colnames(abs.paths)) # drop _Abundance.RELAB from patient ID names in columns
colnames(abs.paths)<-gsub("S_(A0.*)","\\1", colnames(abs.paths)) # drop _Abundance.RELAB from patient ID names in columns
colnames(abs.paths)<-gsub("(A02.*_2023)_(.*)_(.*)","\\1\\.\\2\\.\\3", colnames(abs.paths)) # drop _Abundance.RELAB from patient ID names in columns

rownames(abs.paths)<-abs.paths$Pathway
# change A020_2023_09_15 to A020_2023.09.15
abs.gen.paths<-abs.paths[!grepl("\\|", abs.paths$Pathway),] # keeping only general paths, not species-specific pathways
abs.gen.paths[1:5,1:5]
dim(abs.gen.paths)
# Total pathways in all ABS patients: 11938, general pathways that are not species specific: 471

# drop the unmapped and unintegrated pathway data 
abs.gen.paths<-abs.gen.paths[!rownames(abs.gen.paths) %in% c("UNMAPPED","UNINTEGRATED"),]

# transpose df for merging later
abs.gen.paths.t<-as.data.frame(t(abs.gen.paths[,-1]))
abs.gen.paths.t$SampleID<-rownames(abs.gen.paths.t)

#### Import ABS Data from HuManN - KOs ####

# then import KO gene family assignment data of ABS patients
abs.KOs <- as.data.frame(read.delim('data/Metagenomes/SeqProcessing/Humann3_Results/output_merged/output_merged_relab/ABS_all_ko_genefamilies_relab.tsv', sep="\t"))
# ^ contains relative abundance of pathways found with raw, mgm reads using BioBakery 3 suite (mainly HuManN)
# let's clean this file up just a little bit..
colnames(abs.KOs)[which(colnames(abs.KOs) == "X..Gene.Family")]<-gsub("X..","", colnames(abs.KOs[1])) # clean up Pathway col name
colnames(abs.KOs)<-gsub("(.*)_S.*_L00.*_Abundance.RELAB","\\1", colnames(abs.KOs)) # drop _Abundance.RELAB from patient ID names in columns
colnames(abs.KOs)<-gsub("S_(A0.*)","\\1", colnames(abs.KOs)) # drop _Abundance.RELAB from patient ID names in columns
colnames(abs.KOs)<-gsub("(A02.*_2023)_(.*)_(.*)","\\1\\.\\2\\.\\3", colnames(abs.KOs)) # drop _Abundance.RELAB from patient ID names in columns
rownames(abs.KOs)<-abs.KOs$Gene.Family

abs.gen.KOs<-abs.KOs[!grepl("\\|", abs.KOs$Gene.Family),] # keeping only general KOs, not species-specific pathways

# drop the unmapped and ungrouped KO data 
abs.gen.KOs<-abs.gen.KOs[!rownames(abs.gen.KOs) %in% c("UNMAPPED","UNGROUPED"),]

# transpose df for merging later
abs.gen.KOs.t<-as.data.frame(t(abs.gen.KOs[,-1]))
abs.gen.KOs.t$SampleID<-rownames(abs.gen.KOs.t)

#### Import ABS-HC Data from HuManN - Pathways ####

# first import MetaCyc Pathway data of ABS patients
abs.HC.paths <- as.data.frame(read.delim('data/Metagenomes/Revisions_3.3.2025/data/ABS_HC_Proj_Humann3_output_merged/ABS_HC_Proj_PRJEB32762_all_pathabundance_relab.tsv', sep="\t"))
# ^ contains relative abundance of pathways found with raw, mgm reads using BioBakery 3 suite (mainly HuManN)
# let's clean this file up just a little bit..
colnames(abs.HC.paths)[which(colnames(abs.HC.paths) == "X..Pathway")]<-gsub("X..","", colnames(abs.HC.paths[1])) # clean up Pathway col name
# colnames(abs.HC.paths)[which(colnames(abs.HC.paths) == "X.Pathway")]<-gsub("X.","", colnames(abs.HC.paths)[which(colnames(abs.HC.paths) == "X.Pathway")])
colnames(abs.HC.paths)<-gsub("(.*)_Abundance.RELAB","\\1", colnames(abs.HC.paths)) # drop _Abundance.RELAB from patient ID names in columns

rownames(abs.HC.paths)<-abs.HC.paths$Pathway

abs.HC.prep.paths<-abs.HC.paths[!grepl("\\|", abs.HC.paths$Pathway),] # keeping only preperal paths, not species-specific pathways
abs.HC.prep.paths[1:5,1:5]

# drop the unmapped and unintegrated pathway data 
abs.HC.prep.paths<-abs.HC.prep.paths[!rownames(abs.HC.prep.paths) %in% c("UNMAPPED","UNINTEGRATED"),]

# melt this df to merge with metadata later & change SampleID names!
abs.HC.prep.paths.m<-melt(abs.HC.prep.paths,by="Pathway")
colnames(abs.HC.prep.paths.m)[colnames(abs.HC.prep.paths.m)=="variable"]<-"run_accession"
colnames(abs.HC.prep.paths.m)[colnames(abs.HC.prep.paths.m)=="value"]<-"Pathway_RelAb"
colnames(abs.HC.prep.paths.m) # sanity check

#### Import ABS-HC Data from HuManN - KOs ####

# then import KO gene family assignment data of ABS patients
abs.HC.KOs <- as.data.frame(read.delim('data/Metagenomes/Revisions_3.3.2025/data/ABS_HC_Proj_Humann3_output_merged/ABS_HC_Proj_PRJEB32762_all_ko_genefamilies_relab.tsv', sep="\t"))
# ^ contains relative abundance of pathways found with raw, mgm reads using BioBakery 3 suite (mainly HuManN)
# let's clean this file up just a little bit..
colnames(abs.HC.KOs)[which(colnames(abs.HC.KOs) == "X..Gene.Family")]<-gsub("X..","", colnames(abs.HC.KOs[1])) # clean up Pathway col name
colnames(abs.HC.KOs)<-gsub("(.*)_Abundance.RELAB","\\1", colnames(abs.HC.KOs)) # drop _Abundance.RELAB from patient ID names in columns
rownames(abs.HC.KOs)<-abs.HC.KOs$Gene.Family

abs.HC.prep.KOs<-abs.HC.KOs[!grepl("\\|", abs.HC.KOs$Gene.Family),] # keeping only general KOs, not species-specific pathways
# drop the unmapped and ungrouped KO data 
abs.HC.prep.KOs<-abs.HC.prep.KOs[!rownames(abs.HC.prep.KOs) %in% c("UNMAPPED","UNGROUPED"),]

# melt this df to merge with metadata later & change SampleID names!
abs.HC.prep.KOs.m<-melt(abs.HC.prep.KOs,by="Gene.Family")
colnames(abs.HC.prep.KOs.m)[colnames(abs.HC.prep.KOs.m)=="variable"]<-"run_accession"
colnames(abs.HC.prep.KOs.m)[colnames(abs.HC.prep.KOs.m)=="value"]<-"Gene.Family_RelAb"
colnames(abs.HC.prep.KOs.m) # sanity check

# transpose df for merging later
#abs.HC.gen.KOs.t<-as.data.frame(t(abs.HC.gen.KOs[,-1]))

#### Import ABS-HC Metadata & Fix ABS-HC Data SampleID Names ####
abs.HC.metadata<-as.data.frame(read_xlsx("data/Metagenomes/Revisions_3.3.2025/data/seqIDmatched_clean_modified.xlsx", sheet="seqIDmatched_clean_modified"))
head(abs.HC.metadata)
# ^^ match these data by patient ID that starts with ERR..

# first the pathway data
## create abs.HC.gen.paths
abs.HC.path.rename<-merge(abs.HC.prep.paths.m,abs.HC.metadata,by="run_accession")
abs.HC.gen.paths<-dcast(abs.HC.path.rename,Pathway ~ SampleID, value.var = "Pathway_RelAb")
abs.HC.gen.paths[1:5,1:5]
rownames(abs.HC.gen.paths)<-abs.HC.gen.paths$Pathway

# transpose df for merging later
abs.HC.gen.paths.t<-as.data.frame(t(abs.HC.gen.paths[,-1]))
abs.HC.gen.paths.t[1:5,1:5]
abs.HC.gen.paths.t$SampleID<-rownames(abs.HC.gen.paths.t)

# then KO data
## create abs.HC.gen.KOs
abs.HC.KO.rename<-merge(abs.HC.prep.KOs.m,abs.HC.metadata,by="run_accession")
abs.HC.gen.KOs<-dcast(abs.HC.KO.rename,Gene.Family ~ SampleID, value.var = "Gene.Family_RelAb")
abs.HC.gen.KOs[1:5,1:5]
rownames(abs.HC.gen.KOs)<-abs.HC.gen.KOs$Gene.Family

# transpose df for merging later
abs.HC.gen.KOs.t<-as.data.frame(t(abs.HC.gen.KOs[,-1]))
abs.HC.gen.KOs.t[1:5,1:5]
abs.HC.gen.KOs.t$SampleID<-rownames(abs.HC.gen.KOs.t)

#### Merge ABS & ABS_HC Pathway Data ####

abs.gen.paths.t[1:5,1:5] # ABS pathway relative abundances from Humann
abs.HC.gen.paths.t[1:5,1:5] # ABS-HCs (for revisions) pathway relative abundances from Humann

# do these data frames have enough matching column names for binding?
colnames(abs.gen.paths.t) %in% colnames(abs.HC.gen.paths.t)

# rbind.fill can merge datasets of different kinds if most of their columns match by name!!!!
abs.proj.gen.paths<-rbind.fill(abs.gen.paths.t,abs.HC.gen.paths.t)
abs.proj.gen.paths[1:4,1:10] # sanity check 1
abs.proj.gen.paths[163:167,1:10] # sanity check 2
abs.proj.gen.paths[is.na(abs.proj.gen.paths)] <- 0 # convert all NAs to 0

abs.proj.gen.paths[1:5,1:5]
rownames(abs.proj.gen.paths)<-abs.proj.gen.paths$SampleID

#### Merge ABS & ABS_HC KO Data ####

abs.gen.KOs.t[1:5,1:5] # ABS KO relative abundances from Humann
abs.HC.gen.KOs.t[1:5,1:5] # ABS-HCs (for revisions) KO relative abundances from Humann

# do these data frames have enough matching column names for binding?
colnames(abs.gen.KOs.t) %in% colnames(abs.HC.gen.KOs.t)

# rbind.fill can merge datasets of different kinds if most of their columns match by name!!!!
abs.proj.gen.KOs<-rbind.fill(abs.gen.KOs.t,abs.HC.gen.KOs.t)
abs.proj.gen.KOs[1:4,1:10] # sanity check 1
abs.proj.gen.KOs[163:167,1:10] # sanity check 2
abs.proj.gen.KOs[is.na(abs.proj.gen.KOs)] <- 0 # convert all NAs to 0

rownames(abs.proj.gen.KOs)<-abs.proj.gen.KOs$SampleID

#### Import ABS Project Metadata & Combine with ABS_HC Metadata ####

abs.metadata<-as.data.frame(read_xlsx("data/SubmittedSamplesMaster_metadata_CLHclean.xlsx", sheet="Sheet1_Metadata"))
head(abs.metadata)
colnames(abs.metadata)[colnames(abs.metadata)=="Flare"]<-"Dx2"
head(abs.metadata)

# ^^ match these data by patient ID that starts with ERR..

# create group color palette to merge with abs.metadata
grp.clrs2 = as.data.frame(t(data.frame("Flare"="#8ac926","Remission"="#4d908e","HHP"="purple1")))

grp.clrs2$Dx2<-rownames(grp.clrs2)
colnames(grp.clrs2)[which(names(grp.clrs2) == "V1")] <- "Group_Color"
grp.clrs2

abs.metadata<-merge(abs.metadata,grp.clrs2,by="Dx2")
unique(abs.metadata$Group_Color)

# give ABS-HC a group color for later
head(abs.HC.metadata)
abs.HC.metadata$Group_Color<-"skyblue1"

# compare the data you're about to rbind.fill()
head(abs.metadata)
head(abs.HC.metadata)

# created combined metadata for AUD and ABS comparisons
## only including relevant information like Dx, Dx2, SampleID, Patient, Group_Color, etc
abs.all.metadata<-rbind.fill(abs.metadata[,c(1:2,4:5,7,15)],abs.HC.metadata[,c(1:3,ncol(abs.HC.metadata))])

#### Import ABS (MGM) MetaPhlan Output ####

# first import MetaCyc Pathway data of ABS patients
abs.mgm.tax.wide <- as.data.frame(read.delim('data/Metagenomes/SeqProcessing/ABS_Metaphlan_Results/Batch1to6andnewHC_metaphlan_abundance_table_March2025.csv', sep=","))
# ^ contains relative abundance of taxa found with raw, mgm reads using BioBakery 3 suite (MetaPhlan)
abs.mgm.tax<-melt(abs.mgm.tax.wide,by=vars(Species))

colnames(abs.mgm.tax)[colnames(abs.mgm.tax)=="variable"]<-"ID"
colnames(abs.mgm.tax)[colnames(abs.mgm.tax)=="value"]<-"relative_abundance"
head(abs.mgm.tax)

# pull out first column of MetaPhlan results for splitting
abs.taxa.combined<-data.frame(Species=abs.mgm.tax$Species)

# split clade info by "|" into 8 separate columns per taxonomic level, then assign column names by taxonomic level
abs.taxa.clean1<-as.data.frame(str_split_fixed(abs.taxa.combined$Species, fixed("|"), 8))
names(abs.taxa.clean1)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain")
head(abs.taxa.clean1)

# drop the leading "k__", "p__", etc from the entire data frame
abs.taxa.clean1[,1:8] <- sapply(abs.taxa.clean1[,1:8],function(x) gsub("[:a-z:]__","",as.character(x)))

# cbind the cleaned-up clade name info (8 columns of taxonomic info) to the original data you imported
## order of the dfs has not changed which is how you can use cbind here!
abs.mgm.tax[1:4,] # check that the taxa order match in aud.mgm.tax$clade_name & taxa.clean
abs.taxa.clean1[1:4,]

abs.mgm.tax.clean<-cbind(abs.taxa.clean1,abs.mgm.tax[,-1])
# ^ contains relative abundance of taxa found with raw, mgm reads using BioBakery 3 suite (mainly MetaPhlan)
head(abs.mgm.tax.clean)

#### Merge MetaPhlan RelAb and Metadata Together ####

head(abs.mgm.tax.clean)

# drop all Eukaryotic hits from merged taxa data
"Eukaryota" %in% abs.mgm.tax.clean$Kingdom
abs.mgm.tax.clean<-abs.mgm.tax.clean[!abs.mgm.tax.clean$Kingdom %in% "Eukaryota",]
"Eukaryota" %in% abs.mgm.tax.clean$Kingdom

unique(abs.all.metadata$ID) %in% unique(abs.mgm.tax.clean$ID)
length(unique(abs.all.metadata$ID))
## DATA NOTE: Humann3 data we used all data, including duplicate samples from individuals for ABS-HC data
## for MetaPhlan, Cynthia took the average of the Metaphlan results for each patient and made that one sample
## this is why there are 167 samples in the Humann3 analyses, but only 154 here

which(!unique(abs.all.metadata$ID) %in% unique(abs.mgm.tax.clean$ID))

# merge metadata and pathway relative abundance data together here
abs.tax.meta<-merge(abs.mgm.tax.clean,abs.all.metadata,by="ID") # 
head(abs.tax.meta)

abs.tax.SampID<-unique(abs.tax.meta[,names(abs.tax.meta) %in% c("Kingdom", "Phylum", "Class", "Order", "Family",
                                                                   "Genus","Species","Strain","relative_abundance","SampleID")])
# make sure that we did not lose any samples in this step
length(unique(abs.tax.meta$SampleID))
length(unique(abs.tax.meta$ID))
length(unique(abs.tax.SampID$SampleID))

head(abs.tax.SampID)
# ** will use abs.tax.SampID to create other dfs so we can keep track of samples with legible IDs rather than sequencing IDs


#### Import Gut Brain Modules Data ####

# GBMs - MetaCyc Pathways
gbm.paths<-as.data.frame(read_xlsx("data/Fxns_Of_Interest/GutBrainModules_VallesColomer/GutBrainModules_VC2019.xlsx", sheet="MetaCyc_GBM_List"))
head(gbm.paths)

# GBMs - KO IDs
gbms.kos<-as.data.frame(read_xlsx("data/Fxns_Of_Interest/GutBrainModules_VallesColomer/GutBrainModules_VC2019.xlsx", sheet="GBM_List_KO"))
head(gbms.kos)
which(table(gbms.kos$KO_ID)>1) # confirming that we have no repeat KOs assigned to multiple GBMs

#### Import Ethanol Production Data ####

eth.paths<-as.data.frame(read_xlsx("data/Fxns_Of_Interest/Ethanol_Production_Pathways_MetaCyc.xlsx", sheet="Sheet1"))
head(eth.paths)
which(table(eth.paths$Pathway)>1) # confirming that we have no repeat KOs assigned to different overall eth pathways

eth.KOs<-as.data.frame(read_xlsx("data/Fxns_Of_Interest/Fermentation_Ethanol_KOs.xlsx", sheet="Fermentation_KOs"))
head(eth.KOs)
which(table(eth.KOs$Gene.Family)>1) # confirming that we have no repeat KOs assigned to different overall eth pathways

#### Import KEGG Functions of Interest ####

LPS.kegg<-read.table("data/Fxns_Of_Interest/LPS_Biosynthesis_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
colnames(LPS.kegg)[which(colnames(LPS.kegg) == "KO_ID")]<-"Gene.Family" # rename KO_ID to Gene.Family for easy merging later

quorsens.kegg<-read.table("data/Fxns_Of_Interest/QuorumSensing_KOs_KEGG.txt", header = TRUE, quote = "", sep = "\t", dec = ".")
colnames(quorsens.kegg)[which(colnames(quorsens.kegg) == "KO_ID")]<-"Gene.Family" # rename KO_ID to Gene.Family for easy merging later

antbio.res.kegg<-as.data.frame(read_xlsx("data/Fxns_Of_Interest/AntibioticResistance_KEGG.xlsx", sheet="KOs_DrugGroup"))

#### Arcsin(SqRt()) Transformation of Combined ABS Pathway Relative Abundance Data ####

# now we are going to use the Arcsin(Square Root()) Transformation to transform the pathway relative abundances
## sources for this transformation: LLorens-Rico et al 2021
## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)
## NOTE: even if your data is proportional, it must be between 0-1 for the arcsin(sqrt()) transformation to work!

abs.proj.gen.paths[1:5,1:5] # pathway relative abundance data for all samples, by SampleID
max(abs.proj.gen.paths[,-ncol(abs.proj.gen.paths)]) # sanity check that the largest value in your dt is < 1 for this transformation

# my workaround for keeping SampleID in the first column for easier indexing later
abs.gen.paths.trans<-abs.proj.gen.paths
abs.gen.paths.trans[,-ncol(abs.gen.paths.trans)]<-asin(sqrt(abs.proj.gen.paths[,-ncol(abs.proj.gen.paths)])) # Arcsin(sqrt()) transformation of the pathway rel abs

abs.gen.paths.trans[1:5,1:5]
abs.proj.gen.paths[1:5,1:5] # sanity check

# melt the merged, arcsin(sqrt()) transformed pathway rel ab data
abs.gen.paths.trans.m<-melt(abs.gen.paths.trans,by=c("SampleID"))
head(abs.gen.paths.trans.m)
colnames(abs.gen.paths.trans.m)[which(names(abs.gen.paths.trans.m) == "variable")] <- "Pathway"
colnames(abs.gen.paths.trans.m)[which(names(abs.gen.paths.trans.m) == "value")] <- "Path_Transformed_RelAb"
head(abs.gen.paths.trans.m)

# merge the melted arcsin(sqrt()) transformed pathway relative abundance and patient metadata
abs.meta.trans.pth.all<-merge(abs.all.metadata[,-c(3)],abs.gen.paths.trans.m,by="SampleID") # contains arcsin(sqrt()) transformed pathway relative abundance and patient metadata
head(abs.meta.trans.pth.all)
abs.meta.trans.pth.all$SampleID = factor(abs.meta.trans.pth.all$SampleID, levels=unique(abs.meta.trans.pth.all$SampleID[order(abs.meta.trans.pth.all$Patient,abs.meta.trans.pth.all$Dx2)]), ordered=TRUE)

#### Subset Transformed Relative Abundance Data for Ethanol Production Only ####
abs.gen.paths.trans[1:5,1:5] # pathway relative abundance data for all samples, by SampleID

# first we will subset the wide-veresion of the df that contains the Arcsin(sqrt()) transformed relative abundances only for the ethanol production pathways
path.abs.eth.trans<-abs.gen.paths.trans[,colnames(abs.gen.paths.trans) %in% eth.paths$Pathway]
head(path.abs.eth.trans)
path.abs.eth.trans$SampleID<-rownames(path.abs.eth.trans)

# then merge arcsin(sqrt()) transformed relative abundances for only the ethanol production pathways with the patient metadata
p.abs.eth.t.meta<-merge(abs.all.metadata,path.abs.eth.trans,by="SampleID") # contains arcsin(sqrt()) transformed pathway relative abundance and patient metadata
head(p.abs.eth.t.meta)
p.abs.eth.t.meta$SampleID = factor(p.abs.eth.t.meta$SampleID, levels=unique(p.abs.eth.t.meta$SampleID[order(p.abs.eth.t.meta$Patient,p.abs.eth.t.meta$Dx2)]), ordered=TRUE)

# melt the merged, arcsin(sqrt()) transformed pathway rel ab data we just created for visualizing purposes (i.e., the boxplots and)
p.abs.eth.t.meta.m<-melt(p.abs.eth.t.meta,id.vars=c("SampleID","Patient","Dx","Dx2","Group_Color"))
head(p.abs.eth.t.meta.m)
colnames(p.abs.eth.t.meta.m)[which(names(p.abs.eth.t.meta.m) == "variable")] <- "Pathway"
colnames(p.abs.eth.t.meta.m)[which(names(p.abs.eth.t.meta.m) == "value")] <- "Path_Transformed_RelAb"
head(p.abs.eth.t.meta.m)

# a quick look at the transformed data for ethanol production pathways
ggplot(p.abs.eth.t.meta.m, aes(x=reorder(SampleID,Path_Transformed_RelAb), y=Path_Transformed_RelAb, fill=Pathway))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=12),
        axis.text.x = element_text(hjust=1,angle=45), legend.title = element_text(size=13,hjust=0),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1)) + labs(x="Sample ID", y="Relative Abundance", title="Active Alcohol Use vs. Abstinence vs. Healthy Controls",fill="Ethanol Production") +
  facet_wrap(vars(Dx2), scales = "free_x")

#### Subset Transformed Data for GBM Pathways ####
head(abs.gen.paths.trans.m) # contains the melted arcsin(sqrt()) transformed pathway rel ab data
head(gbm.paths) # contains the gut brain modules from the Valles-Colomer 2019 paper

# first we will subset the wide-version of the df that contains the Arcsin(sqrt()) transformed relative abundances only for GBM pathways
abs.p.gbm.trans<-abs.gen.paths.trans[,(colnames(abs.gen.paths.trans) %in% gbm.paths$Pathway)]
head(abs.p.gbm.trans)
abs.p.gbm.trans$SampleID<-rownames(abs.p.gbm.trans)

# merge ^^ with metadata
abs.gbm.pths.t.meta<-merge(abs.p.gbm.trans,abs.all.metadata,by.x="SampleID",by.y="SampleID")

# merge the arc-sin sqrt transformed pathway rel ab data and the gbm pathway data
abs.gbm.pths.trans.m<-merge(abs.gen.paths.trans.m,gbm.paths,by.x="Pathway",by.y="Pathway")

# merge ^^ with metadata
abs.gbm.pths.t.m.meta<-merge(abs.gbm.pths.trans.m,abs.all.metadata,by.x="SampleID",by.y="SampleID")

# create shortened names for MetaCyc pathways for visualization purposes
abs.gbm.pths.t.m.meta$PathShort<-gsub("(.*):.*","\\1",abs.gbm.pths.t.m.meta$Pathway)
## ^ for gsub() - \\1 corresponds to everything in (.*) that preceeded : in pathway name

# a quick look at the transformed data for GBM pathways
ggplot(abs.gbm.pths.t.m.meta, aes(x=reorder(SampleID,Path_Transformed_RelAb), y=Path_Transformed_RelAb, fill=GBM))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=12),
        axis.text.x = element_text(hjust=1,angle=45), legend.title = element_text(size=13,hjust=0),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1)) + labs(x="Sample ID", y="Relative Abundance") +
  facet_wrap(vars(Dx2), scales = "free_x")


#### Arcsin(SqRt()) Transformation of All KO Relative Abundance Data ####

# now we are going to use the Arcsin(Square Root()) Transformation to transform the KO relative abundances
## sources for this transformation: LLorens-Rico et al 2021
## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)
## NOTE: even if your data is proportional, it must be between 0-1 for the arcsin(sqrt()) transformation to work!

abs.proj.gen.KOs[1:5,1:5] # KO relative abundance data for all samples, by SampleID
max(abs.proj.gen.KOs[,-1]) # sanity check that the largest value in your dt is < 1 for this transformation

# my workaround for keeping SampleID in the first column for easier indexing later
abs.gen.KOs.trans<-abs.proj.gen.KOs
abs.gen.KOs.trans[,-1]<-asin(sqrt(abs.proj.gen.KOs[,-1])) # Arcsin(sqrt()) transformation of the KO rel abs

abs.gen.KOs.trans[1:5,1:5]
abs.proj.gen.KOs[1:5,1:5] # sanity check

# melt the merged, arcsin(sqrt()) transformed KO rel ab data
abs.gen.KOs.trans.m<-melt(abs.gen.KOs.trans,by=c("SampleID"))
head(abs.gen.KOs.trans.m)
colnames(abs.gen.KOs.trans.m)[which(names(abs.gen.KOs.trans.m) == "variable")] <- "Gene.Family"
colnames(abs.gen.KOs.trans.m)[which(names(abs.gen.KOs.trans.m) == "value")] <- "KO_Transformed_RelAb"
head(abs.gen.KOs.trans.m)

# merge the melted arcsin(sqrt()) transformed KO relative abundance and patient metadata
abs.meta.trans.KO.all<-merge(abs.all.metadata,abs.gen.KOs.trans.m,by="SampleID") # contains arcsin(sqrt()) transformed KO relative abundance and patient metadata
head(abs.meta.trans.KO.all)
abs.meta.trans.KO.all$SampleID = factor(abs.meta.trans.KO.all$SampleID, levels=unique(abs.meta.trans.KO.all$SampleID[order(abs.meta.trans.KO.all$Patient,abs.meta.trans.KO.all$Dx2)]), ordered=TRUE)

#### Subset Arcsin(SqRt()) Transformed Relative Abundance Data for Ethanol Production Gene.Families Only ####
abs.gen.KOs.trans[1:5,1:5] # KO relative abundance data for all samples, by SampleID

# first we will subset the wide-version of the df that contains the Arcsin(sqrt()) transformed relative abundances only for the ethanol production KOs
KO.abs.eth.trans<-abs.gen.KOs.trans[,colnames(abs.gen.KOs.trans) %in% eth.KOs$Gene.Family]
head(KO.abs.eth.trans)
KO.abs.eth.trans$SampleID<-rownames(KO.abs.eth.trans)

# then merge arcsin(sqrt()) transformed relative abundances for only the ethanol production KOs with the patient metadata
KO.abs.eth.t.meta<-merge(abs.all.metadata,KO.abs.eth.trans,by="SampleID") # contains arcsin(sqrt()) transformed KO relative abundance and patient metadata
head(KO.abs.eth.t.meta)
KO.abs.eth.t.meta$SampleID = factor(KO.abs.eth.t.meta$SampleID, levels=unique(KO.abs.eth.t.meta$SampleID[order(KO.abs.eth.t.meta$Patient,KO.abs.eth.t.meta$Dx2)]), ordered=TRUE)

# melt the merged, arcsin(sqrt()) transformed KO rel ab data we just created for visualizing purposes (i.e., the boxplots and)
KO.abs.eth.t.meta.m<-melt(KO.abs.eth.t.meta,id.vars=c("SampleID","Patient","Dx","Dx2","Group_Color"))
head(KO.abs.eth.t.meta.m)
colnames(KO.abs.eth.t.meta.m)[which(names(KO.abs.eth.t.meta.m) == "variable")] <- "Gene.Family"
colnames(KO.abs.eth.t.meta.m)[which(names(KO.abs.eth.t.meta.m) == "value")] <- "KO_Transformed_RelAb"
head(KO.abs.eth.t.meta.m)

# merge with the Fermentation KO metadata to keep track of KO functions
KO.fxn.abs.eth.t.meta.m<-merge(KO.abs.eth.t.meta.m,eth.KOs,by=c("Gene.Family"))

# a quick look at the transformed data for ethanol production KOs
ggplot(KO.fxn.abs.eth.t.meta.m, aes(x=reorder(SampleID,KO_Transformed_RelAb), y=KO_Transformed_RelAb, fill=KO_Function))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=12),
        axis.text.x = element_text(hjust=1,angle=45), legend.title = element_text(size=13,hjust=0),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1)) + labs(x="Sample ID", y="Relative Abundance", title="Active Alcohol Use vs. Abstinence vs. Healthy Controls",fill="Ethanol Production KOs") +
  facet_wrap(vars(Dx2), scales = "free_x")

#### Subset Transformed Data for GBM KOs ####
head(abs.gen.KOs.trans.m) # contains the melted arcsin(sqrt()) transformed KO rel ab data
head(gbms.kos) # contains the gut brain modules from the Valles-Colomer 2019 paper

# merge the arc-sin sqrt transformed KO rel ab data and the gbm KO data
abs.gbms.kos.trans.m<-merge(abs.gen.KOs.trans.m,gbms.kos,by.x="Gene.Family",by.y="Gene.Family")

# create long (dcast) version of abs.gbms.kos.trans.m
abs.gbms.kos.trans.table<-dcast(abs.gbms.kos.trans.m, SampleID ~ Gene.Family, value.var = "KO_Transformed_RelAb")
rownames(abs.gbms.kos.trans.table)<-abs.gbms.kos.trans.table$SampleID
abs.gbms.kos.trans.table[1:6,1:6]

# create merged transformed KO + metadata where KOs are columns
abs.gbms.kos.t.meta<-merge(abs.gbms.kos.trans.table,abs.all.metadata,by.x="SampleID",by.y="SampleID")

# merge melted transformed KO data with metadata (with only one column that holds KO IDs)
abs.gbms.kos.t.meta.m<-merge(abs.gbms.kos.trans.m,abs.all.metadata,by.x="SampleID",by.y="SampleID")

# a quick look at the transformed data for GBM KOs
ggplot(abs.gbms.kos.t.meta.m, aes(x=reorder(SampleID,KO_Transformed_RelAb), y=KO_Transformed_RelAb, fill=GBM))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=12),
        axis.text.x = element_text(hjust=1,angle=45), legend.title = element_text(size=13,hjust=0),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2)) + labs(x="Sample ID", y="Relative Abundance", title="Active Alcohol Use vs. Abstinence vs. Healthy Controls",fill="GBM") +
  facet_wrap(vars(Dx2), scales = "free_x")

#### Subset Transformed Data for Antibiotic Resistance KOs ####
head(abs.gen.KOs.trans.m) # contains the melted arcsin(sqrt()) transformed KO rel ab data
head(antbio.res.kegg) # antibiotic-resistance KOs

# merge the arc-sin sqrt transformed KO rel ab data and the gbm KO data
abs.AR.kos.trans.m<-merge(abs.gen.KOs.trans.m,antbio.res.kegg,by.x="Gene.Family",by.y="Gene.Family")

# create long (dcast) version of abs.AR.kos.trans.m
abs.AR.kos.trans.table<-dcast(abs.AR.kos.trans.m, SampleID ~ Gene.Family, value.var = "KO_Transformed_RelAb")
rownames(abs.AR.kos.trans.table)<-abs.AR.kos.trans.table$SampleID
abs.AR.kos.trans.table[1:6,]

# create merged transformed KO + metadata where KOs are columns
abs.AR.kos.t.meta<-merge(abs.AR.kos.trans.table,abs.all.metadata,by.x="SampleID",by.y="SampleID")

# merge melted transformed KO data with metadata (with only one column that holds KO IDs)
abs.AR.kos.t.meta.m<-merge(abs.AR.kos.trans.m,abs.all.metadata,by.x="SampleID",by.y="SampleID")

# a quick look at the transformed data for antibiotic resistance KOs
ggplot(abs.AR.kos.t.meta.m, aes(x=reorder(SampleID,KO_Transformed_RelAb), y=KO_Transformed_RelAb, fill=DrugGroup))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=12),
        axis.text.x = element_text(hjust=1,angle=45), legend.title = element_text(size=13,hjust=0),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1)) + labs(x="Sample ID", y="Relative Abundance", title="Active Alcohol Use vs. Abstinence vs. Healthy Controls",fill="Drug Group") +
  facet_wrap(vars(Dx2), scales = "free_x")


#### Subset Transformed Data for Quorum Sensing KOs ####
head(abs.gen.KOs.trans.m) # contains the melted arcsin(sqrt()) transformed KO rel ab data
head(quorsens.kegg) # antibiotic-resistance KOs

# merge the arc-sin sqrt transformed KO rel ab data and the gbm KO data
abs.quors.kos.trans.m<-merge(abs.gen.KOs.trans.m,quorsens.kegg,by.x="Gene.Family",by.y="Gene.Family")

# create long (dcast) version of abs.quors.kos.trans.m
abs.quors.kos.trans.table<-dcast(abs.quors.kos.trans.m, SampleID ~ Gene.Family, value.var = "KO_Transformed_RelAb")
rownames(abs.quors.kos.trans.table)<-abs.quors.kos.trans.table$SampleID
abs.quors.kos.trans.table[1:6,]

# create merged transformed KO + metadata where KOs are columns
abs.quors.kos.t.meta<-merge(abs.quors.kos.trans.table,abs.all.metadata,by.x="SampleID",by.y="SampleID")

# merge melted transformed KO data with metadata (with only one column that holds KO IDs)
abs.quors.kos.t.meta.m<-merge(abs.quors.kos.trans.m,abs.all.metadata,by.x="SampleID",by.y="SampleID")

# a quick look at the transformed data for antibiotic resistance KOs
ggplot(abs.quors.kos.t.meta.m, aes(x=reorder(SampleID,KO_Transformed_RelAb), y=KO_Transformed_RelAb,fill=KO_Function))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=12),
        axis.text.x = element_text(hjust=1,angle=45), legend.title = element_text(size=13,hjust=0),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1)) + labs(x="Sample ID", y="Relative Abundance", title="Active Alcohol Use vs. Abstinence vs. Healthy Controls",fill="KO_Function") +
  facet_wrap(vars(Dx2), scales = "free_x")



#### Save Data to Rdata File for Downstream Analyses ####

save.image("data/AUD_vs_ABS/AUD_vs_ABS_Transformed_HumannOutput.Rdata") # import into other scripts to pick up here for downstream analyses

#### AUD + ABS Combined Data & Object Descriptions ####
# # ^^ what is important in here...
# anything with *.gen.paths includes only general Humann pathways, NOT species-specific pathways!
# ^^ only considering general pathways and general KOs, not species-specific results from Humann

# Notes
# anything with *.gen.paths includes only general Humann pathways, NOT species-specific pathways!
# anything with abs. means that this includes combined AUD + ABS Humann3 Outputs!
## .eth = ethanol production pathways via fermentation
## .t. or trans = contains arcsin(square root()) transformed relative abundance data
## .m = melted data frame
## AUD = alcohol use disorder; W1 = active alcohol use, W2 = after two weeks of abstinence
## Flare = ABS flare; Remission = remission from ABS flare
## HC = healthy control (AUD project); HHP = healthy household partner (ABS project); ABS_HC = newly added HCs for ABS project
## GBM = gut-brain module
## KO.fxn.eth = includes KO functions for ethanol production via fermentation

### Pathway-Related Objects

abs.proj.gen.paths[1:4,1:4] ## MetaCyc pathway relative abundance for AUD & ABS samples (with updated SampleIDs for AUD samples)
abs.gen.paths.trans[1:4,1:4] ## arcsin(sqrrt()) transformed pathway relative abundances for AUD + ABS data

head(abs.gen.paths.trans.m) ## melted version of abs.gen.paths.trans
head(abs.meta.trans.pth.all) ## metadata + arcsin(sqrrt()) transformed Pathway relative abundances

path.abs.eth.trans[1:4,1:4] ## subsetted arcsin(sqrrt()) transformed Pathway relative abundances – only ethanol production pathways
p.abs.eth.t.meta[1:4,] ## path.eth.trans + metadata
p.abs.eth.t.meta.m[1:4,] ## melted version of p.eth.trans.meta

abs.p.gbm.trans[1:4,] ## table aka cast format of GBM-only arcsin(sqrt()) transformed pathways RelAbs
abs.gbm.pths.t.meta[1:4,] ## metadata + path.gbm.trans
abs.gbm.pths.trans.m[1:4,] ## arcsin(square root()) transformed Pathway relative abundances
abs.gbm.pths.t.m.meta[1:4,] ## metadata + gbm.pths.trans.melt

### KO-Related Objects

## all KOs
abs.proj.gen.KOs[1:5,1:5] # KO relative abundance data for all samples, by SampleID

abs.gen.KOs.trans[1:5,1:5] ## arcsin(square root()) transformed, KO relative abundances
abs.gen.KOs.trans.m[1:4,] ## melted version of abs.gen.KOs.trans
abs.meta.trans.KO.all[1:4,] ## metadata + arcsin(sqrt()) transformed KO RelAbs merged

## ethanol fermentation KOs
KO.abs.eth.trans[1:4,1:4] ## subsetted, arcsin(sqrt()) transformed KO relative abundances – only fermentation KOs

KO.abs.eth.t.meta[1:4,] ## KO.eth.trans and metadata merged into one df
KO.abs.eth.t.meta.m[1:4,] ## melted version of KO.eth.trans.meta  for merging

KO.fxn.abs.eth.t.meta.m[1:4,] ## contains arcsin(sqrt()) transformed KO rel abs + metadata + KO function information

## GBM KOs
abs.gbms.kos.trans.table[1:4,] ## GBM KOs where KOs are columns, not melted df
abs.gbms.kos.trans.m[1:4,] ## merged GBM KOs + arcsin(sqrt()) transformed KO relative abundance data

abs.gbms.kos.t.meta[1:4,] ## merged abs.gbms.kos.trans.table + metadata
# ^ columns are separate KOs
abs.gbms.kos.t.meta.m[1:4,] ## merged abs.gbms.kos.trans.m + metadata
# ^ has one column with all KO IDs (melted version)

## Antibiotic Resistance KOs
abs.AR.kos.trans.table[1:4,] ## Antibiotic Resistance KOs where KOs are columns, not melted df
abs.AR.kos.trans.m[1:4,] ## merged Antibiotic Resistance KOs + arcsin(sqrt()) transformed KO relative abundance data

abs.AR.kos.t.meta[1:4,] ## merged abs.AR.kos.trans.table + metadata
# ^ columns are separate KOs
abs.AR.kos.t.meta.m[1:4,] ## merged abs.AR.kos.trans.m + metadata
# ^ has one column with all KO IDs (melted version)

## Quorum Sensing KOs
abs.quors.kos.trans.table[1:4,] ## Quorum Sensing KOs where KOs are columns, not melted df
abs.quors.kos.trans.m[1:4,] ## merged Quorum Sensing KOs + arcsin(sqrt()) transformed KO relative abundance data

abs.quors.kos.t.meta[1:4,] ## merged abs.quors.kos.trans.table + metadata
# ^ columns are separate KOs
abs.quors.kos.t.meta.m[1:4,] ## merged abs.quors.kos.trans.m + metadata
# ^ has one column with all KO IDs (melted version)


## NOTE - OBJECTS BELOW ARE ONLY AUD PROJECT, DO NOT INCLUDE ABS DATA!
## objects that were imported from "data/AUD_HC_Pathway_RelativeAbundances_HumannOutput_BeforeTransformation.Rdata" are below

# # ^^ what is important in here...
# anything with *.gen.paths includes only general Humann pathways, NOT species-specific pathways!
# ^^ only considering general pathways and general KOs, not species-specific results from Humann

# gut-brain module-related objects:
head(gbm.paths) # Valles-Colomer 2019 paper identifying gut brain modules, connecting these to MetaCyc pathways
head(gbms.kos) # Valles-Colomer 2019 paper identifying gut brain modules, connecting these to KEGG Orthology IDs (KO)

# ethanol production via fermentation objects
head(eth.paths)
head(eth.KOs) # data frame of genes and KOs involved in various types of ethanol fermentation

# Metadata objects:
head(metadata.simple) # metadata variables of interest for AUD W1W2 patients and HC patients
head(filt.metadata) # project metadata excluding AUD samples that have only W1 or W2
head(w1.w2.samps.only) # list of patients that have both W1 and W2 measurements
head(w1.w2.samps.only.clean) # list of patients that have both W1 and W2 measurements, and their sequencing sample IDs that we use to merge

# Raw Pathway Relative Abundances (i.e., the output directly from Humann as is) objects:
head(all.gen.paths) # df of combined Humann outputs for HC patients and AUD patients that have both W1 and W2 measurements
t.all.gen.paths[1:5,1:5] # a transposed version of all.gen.paths

head(meta.path.RA.all) # combined metadata + Pathway RelAbs for patients
head(all.gen.paths.SampID) # new data frame where rows are new SampleIDs that are easier to read, and columns are pathways

head(aud.w1w2.paths) # pathway relative abundances for AUD patients with W1 and W2 measurements only

# Raw KO Relative Abundances
head(all.gen.KOs) # df of combined Humann outputs for HC patients and AUD patients that have both W1 and W2 measurements
t.all.gen.KOs[1:5,1:5] # a transposed version of all.gen.KOs

head(meta.KO.RA.all) # combined metadata + KO RelAbs for patients
head(all.gen.KOs.SampID) # new data frame where rows are new SampleIDs that are easier to read, and columns are KOss

head(aud.w1w2.KOs) # KOs relative abundances for AUD patients with W1 and W2 measurements only

# less important objects loaded but good to be aware of...
head(aud.gen.paths) # all AUD patients' Humann pathway relative abundances, regardless of having W1 or W2 or both W1 and W2 measurements
head(hc.paths) # all HC patients' Humann pathway relative abundances
head(all.gen.paths.melt) # melted version of all.gen.paths

head(aud.gen.KOs) # all AUD patients' Humann KO relative abundances, regardless of having W1 or W2 or both W1 and W2 measurements
head(hc.KOs) # all HC patients' Humann KO relative abundances
head(all.gen.KOs.melt) # melted version of all.gen.KOs

head(metadata) # all project metadata

# Notes
# anything with *.gen.paths includes only general Humann pathways, NOT species-specific pathways!
## .eth = of interest; indicates pathways or KOs relating to ethanol production via fermentation
## AUD = alcohol use disorder; W1 = active alcohol use, W2 = after two weeks of abstinence
## HC = healthy control
## w1.w2 = only patients who have both W1 and W2 measurements included
## GBM = gut-brain module
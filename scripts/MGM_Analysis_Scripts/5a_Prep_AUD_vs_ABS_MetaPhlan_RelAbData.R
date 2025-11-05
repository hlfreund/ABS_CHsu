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
})

# #### Import AUD (MGM) MetaPhlan Output ####
# 
# # first import data of MGM patients
# aud.mgm.tax <- as.data.frame(read.delim('data/Metagenomes/AUD_Proj_MetaPhlan_merged_abundance_table.txt', sep="\t"))
# head(aud.mgm.tax$SampleID)
# colnames(aud.mgm.tax)[which(names(aud.mgm.tax) == "SampleID")] <- "ID" # for merging with metadata and getting real sampleIDs later
# 
# # pull out first column of MetaPhlan results for splitting
# taxa.combined<-data.frame(clade_name=aud.mgm.tax$clade_name)
# 
# # split clade info by "|" into 8 separate columns per taxonomic level, then assign column names by taxonomic level
# taxa.clean<-as.data.frame(str_split_fixed(taxa.combined$clade_name, fixed("|"), 8))
# names(taxa.clean)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain")
# 
# # drop the leading "k__", "p__", etc from the entire data frame
# taxa.clean[,1:8] <- sapply(taxa.clean[,1:8],function(x) gsub("[:a-z:]__","",as.character(x)))
# 
# # cbind the cleaned-up clade name info (8 columns of taxonomic info) to the original data you imported
# ## order of the dfs has not changed which is how you can use cbind here!
# aud.mgm.tax[1:4,1:4] # check that the taxa order match in aud.mgm.tax$clade_name & taxa.clean
# taxa.clean[1:4,]
# 
# aud.mgm.tax.clean1<-cbind(taxa.clean,aud.mgm.tax[,-1])
# # ^ contains relative abundance of taxa found with raw, mgm reads using BioBakery 3 suite (mainly MetaPhlan)
# head(aud.mgm.tax.clean1)
# 
# aud.mgm.tax.clean<-subset(aud.mgm.tax.clean1, select=-c(NCBI_tax_id,additional_species))
# head(aud.mgm.tax.clean)
# 
# 
# #### Import AUD-HC (MGM) MetaPhlan Output ####
# 
# # first import MetaCyc Pathway data of ABS patients
# aud.hc.mgm.tax.wide <- as.data.frame(read.delim('data/240528_humann_output_merged_healthy/all_metaphlan_bugs_list_240528_healthycontrols.csv', sep=","))
# # ^ contains relative abundance of taxa found with raw, mgm reads using BioBakery 3 suite (MetaPhlan)
# aud.hc.mgm.tax<-melt(aud.hc.mgm.tax.wide,by=vars(Species))
# 
# colnames(aud.hc.mgm.tax)[colnames(aud.hc.mgm.tax)=="variable"]<-"ID"
# colnames(aud.hc.mgm.tax)[colnames(aud.hc.mgm.tax)=="value"]<-"relative_abundance"
# head(aud.hc.mgm.tax)
# 
# # pull out first column of MetaPhlan results for splitting
# aud.hc.taxa.combined<-data.frame(Species=aud.hc.mgm.tax$Species)
# 
# # split clade info by "|" into 8 separate columns per taxonomic level, then assign column names by taxonomic level
# aud.hc.taxa.clean1<-as.data.frame(str_split_fixed(aud.hc.taxa.combined$Species, fixed("|"), 8))
# names(aud.hc.taxa.clean1)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Strain")
# head(aud.hc.taxa.clean1)
# 
# # drop the leading "k__", "p__", etc from the entire data frame
# aud.hc.taxa.clean1[,1:8] <- sapply(aud.hc.taxa.clean1[,1:8],function(x) gsub("[:a-z:]__","",as.character(x)))
# 
# # cbind the cleaned-up clade name info (8 columns of taxonomic info) to the original data you imported
# ## order of the dfs has not changed which is how you can use cbind here!
# aud.hc.mgm.tax[1:4,] # check that the taxa order match in aud.mgm.tax$clade_name & taxa.clean
# aud.hc.taxa.clean1[1:4,]
# 
# aud.hc.mgm.tax.clean<-cbind(aud.hc.taxa.clean1,aud.hc.mgm.tax[,-1])
# # ^ contains relative abundance of taxa found with raw, mgm reads using BioBakery 3 suite (mainly MetaPhlan)
# head(aud.hc.mgm.tax.clean)
# 

# #### Import AUD Project Metadata ####
# 
# metadata<-as.data.frame(read_xlsx("data/PS_Metadata.xlsx", sheet="AUD_W1W2_HC_Subset"))
# head(metadata)
# # ^^ match these data by patient ID that starts with ERR...
# 
# metadata$SampleID<-as.character(interaction(metadata$Patient,metadata$Dx2,sep="."))
# 
# metadata$Dx2<-factor(metadata$Dx2,levels=c("AUD_W1","AUD_W2","HC"))
# 
# # create group color palette to merge with metadata
# grp.clrs = as.data.frame(t(data.frame("AUD_W1"="#d00000","AUD_W2"="#ffba08","HC"="#3f88c5")))
# 
# grp.clrs$Dx2<-rownames(grp.clrs)
# colnames(grp.clrs)[which(names(grp.clrs) == "V1")] <- "Group_Color"
# grp.clrs
# 
# metadata<-merge(metadata,grp.clrs,by="Dx2")
# unique(metadata$Group_Color)
# 
# metadata.simple<-subset(metadata,select=c("Dx2","Patient","Dx","ID","SampleID","Group_Color"))
# head(metadata.simple)
# 
#### Import ABS (MGM) MetaPhlan Output ####

# first import MetaCyc Pathway data of ABS patients
abs.mgm.tax.wide <- as.data.frame(read.delim('data/Metagenomes/ABS_Metaphlan_Results/Batch1to6andnewHC_metaphlan_abundance_table_March2025.csv', sep=","))
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


#### Import ABS Project Metadata & Combine with ABS-HC Metadata ####

abs.metadata<-as.data.frame(read_xlsx("data/Metagenomes/ABS_Humann3_Results/SubmittedSamplesMaster_metadata_CLHclean.xlsx", sheet="Sheet1_Metadata"))
head(abs.metadata) # using the original ABS metadata Cynthia provided me to match up with sequencing ID names used in relab table
# ^^ match these data by patient ID
colnames(abs.metadata)[colnames(abs.metadata)=="Flare"]<-"Dx2" # change Flare column name to Dx2

# create group color palette to merge with abs.metadata
grp.clrs2 = as.data.frame(t(data.frame("Flare"="#8ac926","Remission"="#52b788","HHP"="purple1")))

grp.clrs2$Dx2<-rownames(grp.clrs2)
colnames(grp.clrs2)[which(names(grp.clrs2) == "V1")] <- "Group_Color"
grp.clrs2

abs.metadata<-merge(abs.metadata,grp.clrs2,by="Dx2")
unique(abs.metadata$Group_Color)

# now import ABS-HC patient metadata
abs.HC.metadata<-as.data.frame(read_xlsx("data/Metagenomes/ABS_Humann3_Results/revisions_output_merged_relab/seqIDmatched_clean_modified.xlsx", sheet="seqIDmatched_clean_nodups"))
head(abs.HC.metadata)
# ^^ match these data with the abs.HC.metadata$Combo_ID column & metaphlan table

# give ABS-HC a group color for later
head(abs.HC.metadata)
abs.HC.metadata$Group_Color<-"skyblue1"

# compare the data you're about to rbind.fill()
head(abs.metadata)
head(abs.HC.metadata)

# created combined metadata for AUD and ABS comparisons
## only including relevant information like Dx, Dx2, SampleID, Patient, Group_Color, etc
comb.metadata<-rbind.fill(abs.metadata,abs.HC.metadata[,c(1:4,(ncol(abs.HC.metadata)-1):ncol(abs.HC.metadata))])

comb.metadata$Dx2 = factor(comb.metadata$Dx2, levels=c("Flare","Remission","HHP","ABS_HC"), ordered=TRUE)

#### Merge MetaPhlan RelAb and Metadata Together ####

# head(aud.mgm.tax.clean)
# head(aud.hc.mgm.tax.clean)
head(abs.mgm.tax.clean)

comb.mgm.tax.all<-rbind.fill(abs.mgm.tax.clean)
length(unique(comb.mgm.tax.all$ID))

# drop all Eukaryotic hits from merged taxa data
"Eukaryota" %in% comb.mgm.tax.all$Kingdom
comb.mgm.tax.clean<-comb.mgm.tax.all[!comb.mgm.tax.all$Kingdom %in% "Eukaryota",]
"Eukaryota" %in% comb.mgm.tax.clean$Kingdom

unique(comb.metadata$ID) %in% unique(comb.mgm.tax.clean$ID)
length(unique(comb.metadata$ID))
## DATA NOTE: Humann3 data we used all data, including duplicate samples from individuals for ABS-HC data
## for MetaPhlan, Cynthia took the average of the Metaphlan results for each patient and made that one sample
## this is why there are 167 samples in the Humann3 analyses, but only 154 here

which(!unique(comb.metadata$ID) %in% unique(comb.mgm.tax.clean$ID))

# merge metadata and pathway relative abundance data together here
comb.tax.meta<-merge(comb.mgm.tax.clean,comb.metadata,by="ID") # 
head(comb.tax.meta)

comb.tax.SampID<-unique(comb.tax.meta[,names(comb.tax.meta) %in% c("Kingdom", "Phylum", "Class", "Order", "Family",
                                                                "Genus","Species","Strain","relative_abundance","SampleID")])
# make sure that we did not lose any samples in this step
length(unique(comb.tax.meta$SampleID))
length(unique(comb.tax.meta$ID))
length(unique(comb.tax.SampID$SampleID))

head(comb.tax.SampID)
# ** will use comb.tax.SampID to create other dfs so we can keep track of samples with legible IDs rather than sequencing IDs

#### Create Kingdom x Relative Abundance Table ####

# Filter out Genus Relative Abundance

# drop empty Phylums AND keep empty Class, Order, Family, Genera, Species and Strain columns -- prevents repeating relative abundance values!
king.only.clean <- comb.tax.SampID[ which(comb.tax.SampID$Kingdom!="" & comb.tax.SampID$Phylum=="" & comb.tax.SampID$Class=="" & comb.tax.SampID$Order=="" & comb.tax.SampID$Family=="" & comb.tax.SampID$Genus=="" & comb.tax.SampID$Species=="" & comb.tax.SampID$Strain=="") , ]
head(king.only.clean)

# keep genera with relative abundance of >= 1%
#king.only.filt<-king.only.clean[king.only.clean$relative_abundance>=1.0,]

# create JUST class relative abundance table with newly filtered order relative abundance df
mgm.king.relab<-as.data.frame(dcast(king.only.clean,SampleID~Kingdom,value.var = "relative_abundance",fun.aggregate=sum))
mgm.king.relab[1:8,] # a preview of the data frame we just made
rownames(mgm.king.relab)<-mgm.king.relab$SampleID

# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.king.relab[,-1])
mgm.king.relab$SampleID[which(rowSums(mgm.king.relab[,-1])>100)]
mgm.king.relab$SampleID[which(rowSums(mgm.king.relab[,-1])<50)]

## MetaPhlan Relative Abundance Notes!!!
# for samples that rowSums < 100
## this indicates that the remaining reads that would up to 100% were not assigned at the genus level
# for samples that rowSums > 100
## these were samples that had duplicated sequence data, which Cynthia took the average of
## taking the averages led to the > 100 RelAbs


#### Create Phylum x Relative Abundance Table ####

# Filter out Genus Relative Abundance
phy.spec.typ<-subset(comb.tax.SampID,select=c(Phylum,Class,Order,Family,Genus,Species,Strain,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...
head(phy.spec.typ)

# drop empty Phylums AND keep empty Class, Order, Family, Genera, Species and Strain columns -- prevents repeating relative abundance values!
phy.only.clean <- phy.spec.typ[ which(phy.spec.typ$Phylum!="" & phy.spec.typ$Class=="" & phy.spec.typ$Order=="" & phy.spec.typ$Family=="" & phy.spec.typ$Genus=="" & phy.spec.typ$Species=="" & phy.spec.typ$Strain=="") , ]
head(phy.only.clean)

# keep genera with relative abundance of >= 1%
#phy.only.filt<-phy.only.clean[phy.only.clean$relative_abundance>=1.0,]

# create JUST class relative abundance table with newly filtered order relative abundance df
mgm.phy.relab<-as.data.frame(dcast(phy.only.clean,SampleID~Phylum,value.var = "relative_abundance",fun.aggregate=sum))
mgm.phy.relab[1:8,1:6] # a preview of the data frame we just made
rownames(mgm.phy.relab)<-mgm.phy.relab$SampleID

# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.phy.relab[,-1])
mgm.phy.relab$SampleID[which(rowSums(mgm.phy.relab[,-1])>100)]
mgm.phy.relab$SampleID[which(rowSums(mgm.phy.relab[,-1])<50)]

## MetaPhlan Relative Abundance Notes!!!
# for samples that rowSums < 100
## this indicates that the remaining reads that would up to 100% were not assigned at the genus level
# for samples that rowSums > 100
## these were samples that had duplicated sequence data, which Cynthia took the average of
## taking the averages led to the > 100 RelAbs

#### Create Class x Relative Abundance Table ####

# Filter out Genus Relative Abundance
cls.spec.typ<-subset(comb.tax.SampID,select=c(Class,Order,Family,Genus,Species,Strain,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...
head(cls.spec.typ)

# drop empty Classes AND keep empty Order, Family, Genera, Species and Strain columns -- prevents repeating relative abundance values!
cls.only.clean <- cls.spec.typ[ which(cls.spec.typ$Class!="" & cls.spec.typ$Order=="" & cls.spec.typ$Family=="" & cls.spec.typ$Genus=="" & cls.spec.typ$Species=="" & cls.spec.typ$Strain=="") , ]
head(cls.only.clean)

# keep genera with relative abundance of >= 1%
#cls.only.filt<-cls.only.clean[cls.only.clean$relative_abundance>=1.0,]

# create JUST class relative abundance table with newly filtered order relative abundance df
mgm.cls.relab<-as.data.frame(dcast(cls.only.clean,SampleID~Class,value.var = "relative_abundance",fun.aggregate=sum))
mgm.cls.relab[1:8,1:6] # a preview of the data frame we just made
rownames(mgm.cls.relab)<-mgm.cls.relab$SampleID

# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.cls.relab[,-1])
mgm.cls.relab$SampleID[which(rowSums(mgm.cls.relab[,-1])>100)]
mgm.cls.relab$SampleID[which(rowSums(mgm.cls.relab[,-1])<50)]

## MetaPhlan Relative Abundance Notes!!!
# for samples that rowSums < 100
## this indicates that the remaining reads that would up to 100% were not assigned at the genus level
# for samples that rowSums > 100
## these were samples that had duplicated sequence data, which Cynthia took the average of
## taking the averages led to the > 100 RelAbs


#### Create Order x Relative Abundance Table ####

# Filter out Genus Relative Abundance
ord.spec.typ<-subset(comb.tax.SampID,select=c(Order,Family,Genus,Species,Strain,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...
head(ord.spec.typ)

# drop empty Orders AND keep empty Family, Genera, Species and Strain columns -- prevents repeating relative abundance values!
ord.only.clean <- ord.spec.typ[ which(ord.spec.typ$Order!="" & ord.spec.typ$Family=="" & ord.spec.typ$Genus=="" & ord.spec.typ$Species=="" & ord.spec.typ$Strain=="") , ]
head(ord.only.clean)

# keep genera with relative abundance of >= 1%
#ord.only.filt<-ord.only.clean[ord.only.clean$relative_abundance>=1.0,]

# create JUST order relative abundance table with newly filtered order relative abundance df
mgm.ord.relab<-as.data.frame(dcast(ord.only.clean,SampleID~Order,value.var = "relative_abundance",fun.aggregate=sum))
mgm.ord.relab[1:8,1:6] # a preview of the data frame we just made
rownames(mgm.ord.relab)<-mgm.ord.relab$SampleID

# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.ord.relab[,-1])
mgm.ord.relab$SampleID[which(rowSums(mgm.ord.relab[,-1])>100)]
mgm.ord.relab$SampleID[which(rowSums(mgm.ord.relab[,-1])<50)]

## MetaPhlan Relative Abundance Notes!!!
# for samples that rowSums < 100
## this indicates that the remaining reads that would up to 100% were not assigned at the genus level
# for samples that rowSums > 100
## these were samples that had duplicated sequence data, which Cynthia took the average of
## taking the averages led to the > 100 RelAbs

#### Create Family x Relative Abundance Table ####

# Filter out Genus Relative Abundance
fam.spec.typ<-subset(comb.tax.SampID,select=c(Family,Genus,Species,Strain,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...
head(fam.spec.typ)

# drop empty Familys AND keep empty Family, Genera, Species and Strain columns -- prevents repeating relative abundance values!
fam.only.clean <- fam.spec.typ[ which(fam.spec.typ$Family!="" & fam.spec.typ$Genus=="" & fam.spec.typ$Species=="" & fam.spec.typ$Strain=="") , ]
head(fam.only.clean)

# keep genera with relative abundance of >= 1%
#fam.only.filt<-fam.only.clean[fam.only.clean$relative_abundance>=1.0,]

# create JUST order relative abundance table with newly filtered order relative abundance df
mgm.fam.relab<-as.data.frame(dcast(fam.only.clean,SampleID~Family,value.var = "relative_abundance",fun.aggregate=sum))
mgm.fam.relab[1:8,1:6] # a preview of the data frame we just made
rownames(mgm.fam.relab)<-mgm.fam.relab$SampleID

# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.fam.relab[,-1])
mgm.fam.relab$SampleID[which(rowSums(mgm.fam.relab[,-1])>100)]
mgm.fam.relab$SampleID[which(rowSums(mgm.fam.relab[,-1])<50)]

## MetaPhlan Relative Abundance Notes!!!
# for samples that rowSums < 100
## this indicates that the remaining reads that would up to 100% were not assigned at the genus level
# for samples that rowSums > 100
## these were samples that had duplicated sequence data, which Cynthia took the average of
## taking the averages led to the > 100 RelAbs

#### Create Genus x Relative Abundance Table ####

# Filter out Genus Relative Abundance
gen.spec.typ<-subset(comb.tax.SampID,select=c(Genus,Species,Strain,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...

# drop empty Genera AND keep empty Species and Strain columns -- prevents repeating relative abundance values!
gen.only.clean <- gen.spec.typ[ which(gen.spec.typ$Genus!="" & gen.spec.typ$Species=="" & gen.spec.typ$Strain=="") , ]

# keep genera with relative abundance of >= 1%
#gen.only.filt<-gen.only.clean[gen.only.clean$relative_abundance>=1.0,]

# create JUST genus relative abundance table with newly filtered genus relative abundance df
mgm.gen.relab<-as.data.frame(dcast(gen.only.clean,SampleID~Genus,value.var = "relative_abundance",fun.aggregate=sum))
mgm.gen.relab[1:8,1:6] # a preview of the data frame we just made
rownames(mgm.gen.relab)<-mgm.gen.relab$SampleID

# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.gen.relab[,-1])
mgm.gen.relab$SampleID[which(rowSums(mgm.gen.relab[,-1])>100)]

## MetaPhlan Relative Abundance Notes!!!
# for samples that rowSums < 100
## this indicates that the remaining reads that would up to 100% were not assigned at the genus level
# for samples that rowSums > 100
## these were samples that had duplicated sequence data, which Cynthia took the average of
## taking the averages led to the > 100 RelAbs

#### Create Species x Relative Abundance Table ####

# Filter out Species Relative Abundance
spec.typ<-subset(comb.tax.SampID,select=c("Species","Strain","relative_abundance","SampleID"))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...

# drop empty Species AND keep empty Strain columns -- prevents repeating relative abundance values!
spec.only.clean <- spec.typ[ which(spec.typ$Species!="" & spec.typ$Strain=="") , ]

# keep species with relative abundance of >= 1%
#spec.only.filt<-spec.only.clean[spec.only.clean$relative_abundance>=1.0,]

# create JUST specus relative abundance table with newly filtered species relative abundance df
mgm.spec.relab<-as.data.frame(dcast(spec.only.clean,SampleID~Species,value.var = "relative_abundance",fun.aggregate=sum))
mgm.spec.relab[1:4,1:4]
rownames(mgm.spec.relab)<-mgm.spec.relab$SampleID
# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.spec.relab[,-1])

#### Create Strain x Relative Abundance Table ####

# Filter out Strain Relative Abundance
strain.typ<-subset(comb.tax.SampID,select=c("Strain","relative_abundance","SampleID"))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...

# drop empty Species AND keep empty Strain columns -- prevents repeating relative abundance values!
strain.only.clean <- strain.typ[ which(strain.typ$Species!="" & strain.typ$Strain=="") , ]

# keep species with relative abundance of >= 1%
#strain.only.filt<-strain.only.clean[strain.only.clean$relative_abundance>=1.0,]

# create JUST specus relative abundance table with newly filtered species relative abundance df
mgm.strain.relab<-as.data.frame(dcast(strain.only.clean,SampleID~Species,value.var = "relative_abundance",fun.aggregate=sum))
mgm.strain.relab[1:4,1:4]
rownames(mgm.strain.relab)<-mgm.strain.relab$SampleID
# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.strain.relab[,-1])


#### Arcsin(SqRt()) Transformation of PHYLUM Relative Abundance Data ####

# now we are going to use the Arcsin(Square Root()) Transformation to transform the pathway relative abundances
## sources for this transformation: LLorens-Rico et al 2021
## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)
## NOTE: arcsin(square root()) transformation requires that the numbers be between 0 - 1

mgm.phy.relab[1:4,1:4] # phyla relative abundances
max(mgm.phy.relab[,-1]) # make sure that when we divide by 100, our max # is < 1 before transformation!
max(mgm.phy.relab[,-1]/100)

phy.trans.relab<-mgm.phy.relab #creating new data frame for transformation so that we can keep our SampleIDs & order of df
phy.trans.relab[1:4,1:4] # sanity check

# now for the transformation! 
## first we divide all values by 100 so we are working with scaled down #s < 1
## then we do the arcsin(square root()) transformation
phy.trans.relab[,-1]<-asin(sqrt(phy.trans.relab[,-1]/100)) # Arcsin(sqrt()) transformation of the taxa rel abs

phy.trans.relab[1:10,1:5]

# converting any NaN to 0
phy.trans.relab[is.na(phy.trans.relab)] # are there NAs?
#phy.trans.relab[is.na(phy.trans.relab)]<-0



#### Arcsin(SqRt()) Transformation of CLASS Relative Abundance Data ####

# now we are going to use the Arcsin(Square Root()) Transformation to transform the pathway relative abundances
## sources for this transformation: LLorens-Rico et al 2021
## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)
## NOTE: arcsin(square root()) transformation requires that the numbers be between 0 - 1

mgm.cls.relab[1:4,1:4] # class relative abundances
max(mgm.cls.relab[,-1]) # make sure that when we divide by 100, our max # is < 1 before transformation!
max(mgm.cls.relab[,-1]/100)

cls.trans.relab<-mgm.cls.relab #creating new data frame for transformation so that we can keep our SampleIDs & order of df
cls.trans.relab[1:4,1:4] # sanity check

# now for the transformation! 
## first we divide all values by 100 so we are working with scaled down #s < 1
## then we do the arcsin(square root()) transformation
cls.trans.relab[,-1]<-asin(sqrt(cls.trans.relab[,-1]/100)) # Arcsin(sqrt()) transformation of the taxa rel abs

cls.trans.relab[1:10,1:5]

# converting any NaN to 0
cls.trans.relab[is.na(cls.trans.relab)] # are there NAs?
#cls.trans.relab[is.na(cls.trans.relab)]<-0


#### Arcsin(SqRt()) Transformation of ORDER Relative Abundance Data ####

# now we are going to use the Arcsin(Square Root()) Transformation to transform the pathway relative abundances
## sources for this transformation: LLorens-Rico et al 2021
## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)
## NOTE: arcsin(square root()) transformation requires that the numbers be between 0 - 1

mgm.ord.relab[1:4,1:4] # order relative abundances
max(mgm.ord.relab[,-1]) # make sure that when we divide by 100, our max # is < 1 before transformation!
max(mgm.ord.relab[,-1]/100)

ord.trans.relab<-mgm.ord.relab #creating new data frame for transformation so that we can keep our SampleIDs & order of df
ord.trans.relab[1:4,1:4] # sanity check

# now for the transformation! 
## first we divide all values by 100 so we are working with scaled down #s < 1
## then we do the arcsin(square root()) transformation
ord.trans.relab[,-1]<-asin(sqrt(ord.trans.relab[,-1]/100)) # Arcsin(sqrt()) transformation of the taxa rel abs

ord.trans.relab[1:10,1:5]

# converting any NaN to 0
ord.trans.relab[is.na(ord.trans.relab)] # are there NAs?
#ord.trans.relab[is.na(ord.trans.relab)]<-0

#### Arcsin(SqRt()) Transformation of FAMILY Relative Abundance Data ####

# now we are going to use the Arcsin(Square Root()) Transformation to transform the pathway relative abundances
## sources for this transformation: LLorens-Rico et al 2021
## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)
## NOTE: arcsin(square root()) transformation requires that the numbers be between 0 - 1

mgm.fam.relab[1:4,1:4] # family relative abundances
max(mgm.fam.relab[,-1]) # make sure that when we divide by 100, our max # is < 1 before transformation!
max(mgm.fam.relab[,-1]/100)

fam.trans.relab<-mgm.fam.relab #creating new data frame for transformation so that we can keep our SampleIDs & order of df
fam.trans.relab[1:4,1:4] # sanity check

# now for the transformation! 
## first we divide all values by 100 so we are working with scaled down #s < 1
## then we do the arcsin(square root()) transformation
fam.trans.relab[,-1]<-asin(sqrt(fam.trans.relab[,-1]/100)) # Arcsin(sqrt()) transformation of the taxa rel abs

fam.trans.relab[1:10,1:5]

# converting any NaN to 0
fam.trans.relab[is.na(fam.trans.relab)] # are there NAs?
#fam.trans.relab[is.na(fam.trans.relab)]<-0


#### Arcsin(SqRt()) Transformation of GENUS Relative Abundance Data ####

# now we are going to use the Arcsin(Square Root()) Transformation to transform the pathway relative abundances
## sources for this transformation: LLorens-Rico et al 2021
## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)
## NOTE: arcsin(square root()) transformation requires that the numbers be between 0 - 1

mgm.gen.relab[1:4,1:4] # genera relative abundances
max(mgm.gen.relab[,-1]) # make sure that when we divide by 100, our max # is < 1 before transformation!
max(mgm.gen.relab[,-1]/100)

gen.trans.relab<-mgm.gen.relab #creating new data frame for transformation so that we can keep our SampleIDs & order of df
gen.trans.relab[1:4,1:4] # sanity check

# now for the transformation! 
## first we divide all values by 100 so we are working with scaled down #s < 1
## then we do the arcsin(square root()) transformation
gen.trans.relab[,-1]<-asin(sqrt(gen.trans.relab[,-1]/100)) # Arcsin(sqrt()) transformation of the taxa rel abs

gen.trans.relab[1:10,1:5]

# converting any NaN to 0
gen.trans.relab[is.na(gen.trans.relab)] # are there NAs?
#gen.trans.relab[is.na(gen.trans.relab)]<-0

#### Arcsin(SqRt()) Transformation of SPECIES Relative Abundance Data ####

# now we are going to use the Arcsin(Square Root()) Transformation to transform the pathway relative abundances
## sources for this transformation: LLorens-Rico et al 2021
## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)

mgm.spec.relab[1:4,1:4] # species relative abundances
max(mgm.spec.relab[,-1]) # make sure that when we divide by 100, our max # is < 1 before transformation!
max(mgm.spec.relab[,-1]/100)

spec.trans.relab<-mgm.spec.relab #creating new data frame for transformation so that we can keep our SampleIDs & order of df
spec.trans.relab[1:4,1:4] # sanity check

# now for the transformation! 
## first we divide all values by 100 so we are working with scaled down #s < 1
## then we do the arcsin(square root()) transformation
spec.trans.relab[,-1]<-asin(sqrt(spec.trans.relab[,-1]/100)) # Arcsin(sqrt()) transformation of the taxa rel abs

spec.trans.relab[1:5,1:5]

# converting any NaN to 0
spec.trans.relab[is.na(spec.trans.relab)] # are there NAs?
#spec.trans.relab[is.na(spec.trans.relab)]<-0


#### Visualize Kingdom (Raw RelAb) by Dx2 ####

# first merge raw kingdom relative abundances + metadata
mgm.king.relab[1:5,]
mgm.king.relab.m<-melt(mgm.king.relab,by=vars("SampleID"))
head(mgm.king.relab.m)
colnames(mgm.king.relab.m)[colnames(mgm.king.relab.m)=="variable"]<-"Kingdom"
colnames(mgm.king.relab.m)[colnames(mgm.king.relab.m)=="value"]<-"RelAb"

head(comb.metadata)

comb.king.RA.meta<-merge(mgm.king.relab.m[mgm.king.relab.m$RelAb>0,],comb.metadata,by="SampleID")
comb.king.RA.meta[1:3,]
comb.king.RA.meta <- comb.king.RA.meta[order(comb.king.RA.meta$RelAb, decreasing = TRUE),]

# visualize the relab data
rawking.bplt1<-ggplot(comb.king.RA.meta, aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Kingdom)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Taxa by Kingdom in Metagenomes", x="", y="Relative Abundance",fill="Kingdom",
       subtitle="Eukaryotic Hits were Removed")

ggsave(rawking.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Kingdom/AUD_vs_ABS_MGM_Metaphlan_Phyla_RawRA_barplot.png", width = 7000,height = 3000,units = "px",dpi = 200,create.dir = TRUE)

#### Visualize Phyla by Dx2 ####

# first merge raw phyla relative abundances + metadata
mgm.phy.relab[1:5,1:5]
mgm.phy.relab.m<-melt(mgm.phy.relab,by=vars("SampleID"))
head(mgm.phy.relab.m)
colnames(mgm.phy.relab.m)[colnames(mgm.phy.relab.m)=="variable"]<-"Phylum"
colnames(mgm.phy.relab.m)[colnames(mgm.phy.relab.m)=="value"]<-"RelAb"

head(comb.metadata)

comb.phy.RA.meta<-merge(mgm.phy.relab.m[mgm.phy.relab.m$RelAb>0,],comb.metadata,by="SampleID")
comb.phy.RA.meta[1:3,]
comb.phy.RA.meta <- comb.phy.RA.meta[order(comb.phy.RA.meta$RelAb, decreasing = TRUE),]

# visualize the relab data
rawphy.bplt1<-ggplot(comb.phy.RA.meta, aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Phylum)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Phyla in Metagenomes", x="", y="Relative Abundance",fill="Phylum")

ggsave(rawphy.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Phylum/AUD_vs_ABS_MGM_Metaphlan_Phyla_RawRA_barplot.png", width = 7000,height = 3000,units = "px",dpi = 200,create.dir = TRUE)

rawphy.bplt2<-ggplot(comb.phy.RA.meta, aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Phylum)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Phyla in Metagenomes", x="", y="Relative Abundance",fill="Phylum")

ggsave(rawphy.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Phylum/AUD_vs_ABS_MGM_Metaphlan_Phyla_RawRA_barplot_b.png", width = 7000,height = 3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
rawphy.ts1<-ggplot(comb.phy.RA.meta, aes(Phylum, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(comb.phy.RA.meta$Group_Color[order(comb.phy.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Phylum", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Phyla")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawphy.ts1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Phylum/AUD_vs_ABS_MGM_Metaphlan_Phylum_RawRA_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, phyla relative abundances + metadata
phy.trans.relab[1:5,1:5]
phy.trans.relab.m<-melt(phy.trans.relab,by=vars("SampleID"))
head(phy.trans.relab.m)
colnames(phy.trans.relab.m)[colnames(phy.trans.relab.m)=="variable"]<-"Phylum"
colnames(phy.trans.relab.m)[colnames(phy.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(comb.metadata)

comb.phy.tRA.meta<-merge(phy.trans.relab.m[phy.trans.relab.m$Transformed_RelAb>0,],comb.metadata,by="SampleID")
comb.phy.tRA.meta[1:3,]
comb.phy.tRA.meta <- comb.phy.tRA.meta[order(comb.phy.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what is the max of the transformed relab?
max(comb.phy.tRA.meta$Transformed_RelAb)

# visualize transformed, relab data
trphy.bplt1<-ggplot(comb.phy.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Phylum)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Phyla in Metagenomes", x="", y="Transformed Relative Abundance",fill="Phylum")

ggsave(trphy.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Phylum/AUD_vs_ABS_MGM_Metaphlan_Phyla_AsinSqrtRA_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trphy.bplt2<-ggplot(comb.phy.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Phylum)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Phyla in Metagenomes", x="", y="Transformed Relative Abundance",fill="Phylum")

ggsave(trphy.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Phylum/AUD_vs_ABS_MGM_Metaphlan_Phyla_AsinSqrtRA_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trphy.bplt3<-ggplot(comb.phy.tRA.meta[comb.phy.tRA.meta$Transformed_RelAb>0.5,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Phylum)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Phyla in Metagenomes", x="", y="Transformed Relative Abundance",fill="Phylum",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.5")

ggsave(trphy.bplt3,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Phylum/AUD_vs_ABS_MGM_Metaphlan_Phyla_AsinSqrtRA_Over0.5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)


#### Visualize Classes by Dx2 ####

# first merge raw class relative abundances + metadata
mgm.cls.relab[1:5,1:5]
mgm.cls.relab.m<-melt(mgm.cls.relab,by=vars("SampleID"))
head(mgm.cls.relab.m)
colnames(mgm.cls.relab.m)[colnames(mgm.cls.relab.m)=="variable"]<-"Class"
colnames(mgm.cls.relab.m)[colnames(mgm.cls.relab.m)=="value"]<-"RelAb"

# convert Humann3 specific labels to Other
mgm.cls.relab.m$Class<-gsub("CFGB.*","Other",mgm.cls.relab.m$Class)
# convert any unclassified taxa as just Unclassified
mgm.cls.relab.m$Class<-gsub(".*_unclassified","Unclassified",mgm.cls.relab.m$Class)

head(comb.metadata)

comb.cls.RA.meta<-merge(mgm.cls.relab.m[mgm.cls.relab.m$RelAb>0,],comb.metadata,by="SampleID")
comb.cls.RA.meta[1:3,]
comb.cls.RA.meta <- comb.cls.RA.meta[order(comb.cls.RA.meta$RelAb, decreasing = TRUE),]

# what are the min, max of the raw RelAbs?
min(comb.cls.RA.meta$RelAb)

# how many colors do you need in the palette?
unique(comb.cls.RA.meta$Class[comb.cls.RA.meta$RelAb>5])

# visualize the relab data
rawcls.bplt1<-ggplot(comb.cls.RA.meta[comb.cls.RA.meta$RelAb>5,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Class)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Classes in Metagenomes", x="", y="Relative Abundance",fill="Class",
       subtitle="Includes Taxa with Relative Abundance > 5%")

ggsave(rawcls.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Class/AUD_vs_ABS_MGM_Metaphlan_Class_RawRA_Over5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# how many colors do you need in the palette?
unique(comb.cls.RA.meta$Class[comb.cls.RA.meta$RelAb>20])

rawcls.bplt2<-ggplot(comb.cls.RA.meta[comb.cls.RA.meta$RelAb>20,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Class)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Classes in Metagenomes", x="", y="Relative Abundance",fill="Class",
       subtitle="Includes Taxa with Relative Abundance > 20%")

ggsave(rawcls.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Class/AUD_vs_ABS_MGM_Metaphlan_Class_RawRA_Over20_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
rawcls.ts1<-ggplot(comb.cls.RA.meta[comb.cls.RA.meta$RelAb>5,], aes(Class, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(comb.cls.RA.meta$Group_Color[order(comb.cls.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle = element_text(size=35)) +
  labs(x="Class", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Classes",subtitle="Includes taxa with Relative Abundance > 5%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawcls.ts1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Class/AUD_vs_ABS_MGM_Metaphlan_Class_RawRA_Over5_taxo.sum.png", width = 5000,height=7000,units = "px",dpi = 200,create.dir = TRUE)

rawcls.ts2<-ggplot(comb.cls.RA.meta[comb.cls.RA.meta$RelAb>10,], aes(Class, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(comb.cls.RA.meta$Group_Color[order(comb.cls.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),plot.subtitle = element_text(size=35)) +
  labs(x="Class", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Classes",subtitle="Includes taxa with Relative Abundance > 10%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawcls.ts2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Class/AUD_vs_ABS_MGM_Metaphlan_Class_RawRA_Over10_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, class relative abundances + metadata
cls.trans.relab[1:5,1:5]
cls.trans.relab.m<-melt(cls.trans.relab,by=vars("SampleID"))
head(cls.trans.relab.m)
colnames(cls.trans.relab.m)[colnames(cls.trans.relab.m)=="variable"]<-"Class"
colnames(cls.trans.relab.m)[colnames(cls.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(comb.metadata)

comb.cls.tRA.meta<-merge(cls.trans.relab.m[cls.trans.relab.m$Transformed_RelAb>0,],comb.metadata,by="SampleID")
comb.cls.tRA.meta[1:3,]
comb.cls.tRA.meta <- comb.cls.tRA.meta[order(comb.cls.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what are the mins, maxs of transformed rel ab data?
max(comb.cls.tRA.meta$Transformed_RelAb)
min(comb.cls.tRA.meta$Transformed_RelAb)

# how many colors do you need in the palette?
unique(comb.cls.tRA.meta$Class[comb.cls.tRA.meta$Transformed_RelAb>0.3])

# visualize transformed, relab data
ggplot(comb.cls.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Classes in Metagenomes", x="", y="Transformed Relative Abundance",fill="Class")

# how many colors do you need in the palette?
unique(comb.cls.tRA.meta$Class[comb.cls.tRA.meta$Transformed_RelAb>0.3])

trcls.bplt1<-ggplot(comb.cls.tRA.meta[comb.cls.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Classes in Metagenomes", x="", y="Transformed Relative Abundance",fill="Class",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trcls.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Class/AUD_vs_ABS_MGM_Metaphlan_Class_AsinSqrtRA_Over0.3_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trcls.bplt2<-ggplot(comb.cls.tRA.meta[comb.cls.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Classes in Metagenomes", x="", y="Transformed Relative Abundance",fill="Class",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trcls.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Class/AUD_vs_ABS_MGM_Metaphlan_Class_AsinSqrtRA_Over0.3_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# how many colors do you need in the palette?
unique(comb.cls.tRA.meta$Class[comb.cls.tRA.meta$Transformed_RelAb>0.6])

trcls.bplt3<-ggplot(comb.cls.tRA.meta[comb.cls.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Classes in Metagenomes", x="", y="Transformed Relative Abundance",fill="Class",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trcls.bplt3,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Class/AUD_vs_ABS_MGM_Metaphlan_Class_AsinSqrtRA_Over0.6_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trcls.bplt4<-ggplot(comb.cls.tRA.meta[comb.cls.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Classes in Metagenomes", x="", y="Transformed Relative Abundance",fill="Class",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trcls.bplt4,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Class/AUD_vs_ABS_MGM_Metaphlan_Class_AsinSqrtRA_Over0.6_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

#### Visualize Orders by Dx2 ####

# first merge raw class relative abundances + metadata
mgm.ord.relab[1:5,1:5]
mgm.ord.relab.m<-melt(mgm.ord.relab,by=vars("SampleID"))
head(mgm.ord.relab.m)
colnames(mgm.ord.relab.m)[colnames(mgm.ord.relab.m)=="variable"]<-"Order"
colnames(mgm.ord.relab.m)[colnames(mgm.ord.relab.m)=="value"]<-"RelAb"

# convert Humann3 specific labels to Other
unique(mgm.ord.relab.m$Order) # to see what to change
mgm.ord.relab.m$Order<-gsub("OFGB.*","Other",mgm.ord.relab.m$Order)
# convert any unclassified taxa as just Unclassified
mgm.ord.relab.m$Order<-gsub(".*_unclassified","Unclassified",mgm.ord.relab.m$Order)

head(comb.metadata)

comb.ord.RA.meta<-merge(mgm.ord.relab.m[mgm.ord.relab.m$RelAb>0,],comb.metadata,by="SampleID")
comb.ord.RA.meta[1:3,]
comb.ord.RA.meta <- comb.ord.RA.meta[order(comb.ord.RA.meta$RelAb, decreasing = TRUE),]

# what are the min, max of the raw RelAbs?
min(comb.ord.RA.meta$RelAb)
max(comb.ord.RA.meta$RelAb)

# how many colors for palette?
unique(comb.ord.RA.meta$Order)

# visualize the relab data
raword.bplt1<-ggplot(comb.ord.RA.meta[comb.ord.RA.meta$RelAb>5,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Order)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Orders in Metagenomes", x="", y="Relative Abundance",fill="Order",
       subtitle="Includes Taxa with Relative Abundance > 5%")

ggsave(raword.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Order/AUD_vs_ABS_MGM_Metaphlan_Order_RawRA_Over5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

raword.bplt2<-ggplot(comb.ord.RA.meta[comb.ord.RA.meta$RelAb>10,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Order)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Orders in Metagenomes", x="", y="Relative Abundance",fill="Order",
       subtitle="Includes Taxa with Relative Abundance > 10%")

ggsave(raword.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Order/AUD_vs_ABS_MGM_Metaphlan_Order_RawRA_Over10_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

raword.bplt3<-ggplot(comb.ord.RA.meta[comb.ord.RA.meta$RelAb>20,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Order)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Orders in Metagenomes", x="", y="Relative Abundance",fill="Order",
       subtitle="Includes Taxa with Relative Abundance > 20%")

ggsave(raword.bplt3,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Order/AUD_vs_ABS_MGM_Metaphlan_Order_RawRA_Over20_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
raword.ts1<-ggplot(comb.ord.RA.meta[comb.ord.RA.meta$RelAb>5,], aes(Order, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(comb.ord.RA.meta$Group_Color[order(comb.ord.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, legend.title = element_text(size=40),
        legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle = element_text(size=35)) +
  labs(x="Order", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Orders",subtitle="Includes taxa with Relative Abundance > 5%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(raword.ts1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Order/AUD_vs_ABS_MGM_Metaphlan_Order_RawRA_Over5_taxo.sum.png", width = 5000,height=7000,units = "px",dpi = 200,create.dir = TRUE)

raword.ts2<-ggplot(comb.ord.RA.meta[comb.ord.RA.meta$RelAb>10,], aes(Order, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(comb.ord.RA.meta$Group_Color[order(comb.ord.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, legend.title = element_text(size=40),
        legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle = element_text(size=35)) +
  labs(x="Order", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Orders",subtitle="Includes taxa with Relative Abundance > 10%")+
  coord_flip()

ggsave(raword.ts2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Order/AUD_vs_ABS_MGM_Metaphlan_Order_RawRA_Over10_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, class relative abundances + metadata
ord.trans.relab[1:5,1:5]
ord.trans.relab.m<-melt(ord.trans.relab,by=vars("SampleID"))
head(ord.trans.relab.m)
colnames(ord.trans.relab.m)[colnames(ord.trans.relab.m)=="variable"]<-"Order"
colnames(ord.trans.relab.m)[colnames(ord.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(comb.metadata)

comb.ord.tRA.meta<-merge(ord.trans.relab.m[ord.trans.relab.m$Transformed_RelAb>0,],comb.metadata,by="SampleID")
comb.ord.tRA.meta[1:3,]
comb.ord.tRA.meta <- comb.ord.tRA.meta[order(comb.ord.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what are the min, max of transformed order relabs?
max(comb.ord.tRA.meta$Transformed_RelAb)
min(comb.ord.tRA.meta$Transformed_RelAb)

# visualize transformed, relab data
ggplot(comb.ord.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Orders in Metagenomes", x="", y="Transformed Relative Abundance",fill="Order")

trord.bplt1<-ggplot(comb.ord.tRA.meta[comb.ord.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Orders in Metagenomes", x="", y="Transformed Relative Abundance",fill="Order",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trord.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Order/AUD_vs_ABS_MGM_Metaphlan_Order_AsinSqrtRA_Over0.3_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trord.bplt2<-ggplot(comb.ord.tRA.meta[comb.ord.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Orders in Metagenomes", x="", y="Transformed Relative Abundance",fill="Order",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trord.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Order/AUD_vs_ABS_MGM_Metaphlan_Order_AsinSqrtRA_Over0.3_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trord.bplt3<-ggplot(comb.ord.tRA.meta[comb.ord.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Orders in Metagenomes", x="", y="Transformed Relative Abundance",fill="Order",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trord.bplt3,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Order/AUD_vs_ABS_MGM_Metaphlan_Order_AsinSqrtRA_Over0.6_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trord.bplt4<-ggplot(comb.ord.tRA.meta[comb.ord.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Orders in Metagenomes", x="", y="Transformed Relative Abundance",fill="Order",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trord.bplt4,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Order/AUD_vs_ABS_MGM_Metaphlan_Order_AsinSqrtRA_Over0.6_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

#### Visualize Families by Dx2 ####

# first merge raw class relative abundances + metadata
mgm.fam.relab[1:5,1:5]
mgm.fam.relab.m<-melt(mgm.fam.relab,by=vars("SampleID"))
head(mgm.fam.relab.m)
colnames(mgm.fam.relab.m)[colnames(mgm.fam.relab.m)=="variable"]<-"Family"
colnames(mgm.fam.relab.m)[colnames(mgm.fam.relab.m)=="value"]<-"RelAb"

# convert Humann3 specific labels to Other
unique(mgm.fam.relab.m$Family) # to see what to change
mgm.fam.relab.m$Family<-gsub("FGB.*","Other",mgm.fam.relab.m$Family)
# convert any unclassified taxa as just Unclassified
mgm.fam.relab.m$Family<-gsub(".*_unclassified","Unclassified",mgm.fam.relab.m$Family)

head(comb.metadata)

comb.fam.RA.meta<-merge(mgm.fam.relab.m[mgm.fam.relab.m$RelAb>0,],comb.metadata,by="SampleID")
comb.fam.RA.meta[1:3,]
comb.fam.RA.meta <- comb.fam.RA.meta[order(comb.fam.RA.meta$RelAb, decreasing = TRUE),]

# what are the min, max of the raw RelAbs?
min(comb.fam.RA.meta$RelAb)
max(comb.fam.RA.meta$RelAb)

# how many colors for palette?
unique(comb.fam.RA.meta$Family)

# visualize the relab data
rawfam.bplt1<-ggplot(comb.fam.RA.meta[comb.fam.RA.meta$RelAb>5,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Family)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Families in Metagenomes", x="", y="Relative Abundance",fill="Family",
       subtitle="Includes Taxa with Relative Abundance > 5%")

ggsave(rawfam.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Family/AUD_vs_ABS_MGM_Metaphlan_Family_RawRA_Over5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawfam.bplt2<-ggplot(comb.fam.RA.meta[comb.fam.RA.meta$RelAb>10,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Family)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Families in Metagenomes", x="", y="Relative Abundance",fill="Family",
       subtitle="Includes Taxa with Relative Abundance > 10%")

ggsave(rawfam.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Family/AUD_vs_ABS_MGM_Metaphlan_Family_RawRA_Over10_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawfam.bplt3<-ggplot(comb.fam.RA.meta[comb.fam.RA.meta$RelAb>20,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Family)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Families in Metagenomes", x="", y="Relative Abundance",fill="Family",
       subtitle="Includes Taxa with Relative Abundance > 20%")

ggsave(rawfam.bplt3,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Family/AUD_vs_ABS_MGM_Metaphlan_Family_RawRA_Over20_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
rawfam.ts1<-ggplot(comb.fam.RA.meta[comb.fam.RA.meta$RelAb>5,], aes(Family, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(comb.fam.RA.meta$Group_Color[order(comb.fam.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40,),
        axis.text = element_text(size=35),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Family", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Families",subtitle="Includes taxa with Relative Abundance > 5%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawfam.ts1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Family/AUD_vs_ABS_MGM_Metaphlan_Family_RawRA_Over10_taxo.sum.png", width = 5000,height=7000,units = "px",dpi = 200,create.dir = TRUE)

rawfam.ts2<-ggplot(comb.fam.RA.meta[comb.fam.RA.meta$RelAb>10,], aes(Family, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(comb.fam.RA.meta$Group_Color[order(comb.fam.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle = element_text(size=35)) +
  labs(x="Family", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Families",subtitle="Includes taxa with Relative Abundance > 10%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawfam.ts2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Family/AUD_vs_ABS_MGM_Metaphlan_Family_RawRA_Over10_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, class relative abundances + metadata
fam.trans.relab[1:5,1:5]
fam.trans.relab.m<-melt(fam.trans.relab,by=vars("SampleID"))
head(fam.trans.relab.m)
colnames(fam.trans.relab.m)[colnames(fam.trans.relab.m)=="variable"]<-"Family"
colnames(fam.trans.relab.m)[colnames(fam.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(comb.metadata)

comb.fam.tRA.meta<-merge(fam.trans.relab.m[fam.trans.relab.m$Transformed_RelAb>0,],comb.metadata,by="SampleID")
comb.fam.tRA.meta[1:3,]
comb.fam.tRA.meta <- comb.fam.tRA.meta[order(comb.fam.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what are the min, max of transformed order relabs?
max(comb.fam.tRA.meta$Transformed_RelAb)
min(comb.fam.tRA.meta$Transformed_RelAb)

# visualize transformed, relab data
ggplot(comb.fam.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Families in Metagenomes", x="", y="Transformed Relative Abundance",fill="Family")

trfam.bplt1<-ggplot(comb.fam.tRA.meta[comb.fam.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Families in Metagenomes", x="", y="Transformed Relative Abundance",fill="Family",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trfam.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Family/AUD_vs_ABS_MGM_Metaphlan_Family_AsinSqrtRA_Over0.3_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trfam.bplt2<-ggplot(comb.fam.tRA.meta[comb.fam.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Families in Metagenomes", x="", y="Transformed Relative Abundance",fill="Family",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trfam.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Family/AUD_vs_ABS_MGM_Metaphlan_Family_AsinSqrtRA_Over0.3_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trfam.bplt3<-ggplot(comb.fam.tRA.meta[comb.fam.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Families in Metagenomes", x="", y="Transformed Relative Abundance",fill="Family",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trfam.bplt3,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Family/AUD_vs_ABS_MGM_Metaphlan_Family_AsinSqrtRA_Over0.6_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trfam.bplt4<-ggplot(comb.fam.tRA.meta[comb.fam.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Families in Metagenomes", x="", y="Transformed Relative Abundance",fill="Family",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trfam.bplt4,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Family/AUD_vs_ABS_MGM_Metaphlan_Family_AsinSqrtRA_Over0.6_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)


#### Visualize Genera by Dx2 ####

# first merge raw class relative abundances + metadata
mgm.gen.relab[1:5,1:5]
mgm.gen.relab.m<-melt(mgm.gen.relab,by=vars("SampleID"))
head(mgm.gen.relab.m)
colnames(mgm.gen.relab.m)[colnames(mgm.gen.relab.m)=="variable"]<-"Genus"
colnames(mgm.gen.relab.m)[colnames(mgm.gen.relab.m)=="value"]<-"RelAb"

# convert Humann3 specific labels to Other
unique(mgm.gen.relab.m$Genus) # to see what to change
mgm.gen.relab.m$Genus<-gsub("GGB.*","Other",mgm.gen.relab.m$Genus)
# convert any unclassified taxa as just Unclassified
mgm.gen.relab.m$Genus<-gsub(".*_unclassified","Unclassified",mgm.gen.relab.m$Genus)

head(comb.metadata)

comb.gen.RA.meta<-merge(mgm.gen.relab.m[mgm.gen.relab.m$RelAb>0,],comb.metadata,by="SampleID")
comb.gen.RA.meta[1:3,]
comb.gen.RA.meta <- comb.gen.RA.meta[order(comb.gen.RA.meta$RelAb, decreasing = TRUE),]

# what are the min, max of the raw RelAbs?
min(comb.gen.RA.meta$RelAb)
max(comb.gen.RA.meta$RelAb)

# how many colors for palette?
unique(comb.gen.RA.meta$Genus)

# visualize the relab data
rawgen.bplt1<-ggplot(comb.gen.RA.meta[comb.gen.RA.meta$RelAb>5,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Genus)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  brightness(scale_fill_discrete_qualitative(palette="Dark 3"),0.95)+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Relative Abundance",fill="Genus",
       subtitle="Includes Taxa with Relative Abundance > 5%") +
  guides(fill=guide_legend(ncol=4))

ggsave(rawgen.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_RawRA_Over5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawgen.bplt2<-ggplot(comb.gen.RA.meta[comb.gen.RA.meta$RelAb>10,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Genus)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Relative Abundance",fill="Genus",
       subtitle="Includes Taxa with Relative Abundance > 10%")

ggsave(rawgen.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_RawRA_Over10_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawgen.bplt3<-ggplot(comb.gen.RA.meta[comb.gen.RA.meta$RelAb>20,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Genus)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Relative Abundance",fill="Genus",
       subtitle="Includes Taxa with Relative Abundance > 20%")

ggsave(rawgen.bplt3,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_RawRA_Over20_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawgen.bplt3<-ggplot(comb.gen.RA.meta[comb.gen.RA.meta$RelAb>30,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Genus)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance (%)",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Relative Abundance",fill="Genus",
       subtitle="Includes Taxa with Relative Abundance > 30%")

ggsave(rawgen.bplt3,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_RawRA_Over30_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
rawgen.ts1<-ggplot(comb.gen.RA.meta[comb.gen.RA.meta$RelAb>5,], aes(Genus, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(comb.gen.RA.meta$Group_Color[order(comb.gen.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle=element_text(size=35)) +
  labs(x="Genera", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Genera",subtitle="Includes taxa with Relative Abundance > 5%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawgen.ts1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_RawRA_Over5_taxo.sum.png", width = 5000,height=7000,units = "px",dpi = 200,create.dir = TRUE)

rawgen.ts2<-ggplot(comb.gen.RA.meta[comb.gen.RA.meta$RelAb>10,], aes(Genus, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(comb.gen.RA.meta$Group_Color[order(comb.gen.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle=element_text(size=35)) +
  labs(x="Genera", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Genera",subtitle="Includes taxa with Relative Abundance > 10%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawgen.ts2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_RawRA_Over10_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, class relative abundances + metadata
gen.trans.relab[1:5,1:5]
gen.trans.relab.m<-melt(gen.trans.relab,by=vars("SampleID"))
head(gen.trans.relab.m)
colnames(gen.trans.relab.m)[colnames(gen.trans.relab.m)=="variable"]<-"Genus"
colnames(gen.trans.relab.m)[colnames(gen.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(comb.metadata)

comb.gen.tRA.meta<-merge(gen.trans.relab.m[gen.trans.relab.m$Transformed_RelAb>0,],comb.metadata,by="SampleID")
comb.gen.tRA.meta[1:3,]
comb.gen.tRA.meta <- comb.gen.tRA.meta[order(comb.gen.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what are the min, max of transformed order relabs?
max(comb.gen.tRA.meta$Transformed_RelAb)
min(comb.gen.tRA.meta$Transformed_RelAb)

# visualize transformed, relab data
ggplot(comb.gen.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Genus")

trgen.bplt1<-ggplot(comb.gen.tRA.meta[comb.gen.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Genus",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trgen.bplt1,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_AsinSqrtRA_Over0.3_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trgen.bplt2<-ggplot(comb.gen.tRA.meta[comb.gen.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Genus",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trgen.bplt2,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_AsinSqrtRA_Over0.3_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trgen.bplt3<-ggplot(comb.gen.tRA.meta[comb.gen.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
  geom_col(position = "stack")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Genus",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trgen.bplt3,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_AsinSqrtRA_Over0.6_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trgen.bplt4<-ggplot(comb.gen.tRA.meta[comb.gen.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
  geom_col(position = "fill")+
  facet_grid(~Dx2, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Transformed Relative Abundance",expand=c(0,0))+
  theme(plot.title = element_text(colour = "black", size=40),
        plot.subtitle = element_text(colour = "black", size=35),
        axis.title.x = element_blank(),
        text = element_text(colour = "black", size=35),
        axis.text=element_text(colour="black", size=35),
        axis.text.x = element_blank(),
        strip.text.x = element_text(size = 21, face = "bold"), 
        strip.background=element_rect(fill="ghostwhite")) +
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Genus",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trgen.bplt4,filename = "figures/AUD_vs_ABS/Taxa_RelAbs/Genus/AUD_vs_ABS_MGM_Metaphlan_Genus_AsinSqrtRA_Over0.6_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)



#### Save Global Envionrment as RData Object ####
save.image("data/AUD_vs_ABS/AUD_vs_ABS_MetaPhlan_BacteriaRelAb_ArcSinTransformed_Ready.Rdata")
# ^ includes all data combined in our global environment

# to load this global environment so that you don't have to run all this code again
## aka quickly upload every object made in this code
# load("data/MetaPhlan_BacteriaRelAb_Filtered_ArcSinTransformed_Ready.Rdata") # load Rdata to global env

# Version Information
sessionInfo()

#### List of Objects & Descriptions ####

# More important data frames to be used for analysis
gen.trans.relab[1:5,1:5] # relative abundance of genera with Sample IDs as rows, genera as columns
gen.trans.relab[1:5,1:5] # arcsine(sqr-root()) transformed relative abundance (genera)

mgm.spec.relab[1:5,1:5] # relative abundance of species with Sample IDs as rows, species as columns
spec.trans.relab[1:5,1:5] # arcsine(sqr-root()) transformed relative abundance (species)

# Less important data frames that were used to reformat the data
head(aud.mgm.tax) # aud.mgm.tax is the original import of the MetaPhlan data
head(aud.mgm.tax.clean) # parsed taxonomy version, where each taxonomic level has its own column now
head(aud.mgm.tax.clean1) # same as aud.mgm.tax.clean1 but without NCBI and additional_species columns

head(gen.spec.typ) # subsetted version of aud.mgm.tax.clean1 but only includes Genus, Species, Strain, relative_abundance, and SampleID columns
head(spec.typ) # subsetted version of aud.mgm.tax.clean1 but only includes Species, Strain, relative_abundance, and SampleID columns (no Genus!)

gen.only.clean[1:5,1:5] # bacterial species relative abundance

spec.only.clean[1:5,] # bacterial species relative abundance

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
# names(taxa.clean)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Clade")
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
# names(aud.hc.taxa.clean1)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Clade")
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
names(abs.taxa.clean1)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Species","Clade")
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

abs.metadata<-as.data.frame(read_xlsx("data/ABS_confirmedonly_metadata.xlsx", sheet="ABS_confirmedonly_metadata"))
head(abs.metadata) # using the original ABS metadata Cynthia provided me to match up with sequencing ID names used in relab table
# ^^ match these data by patient ID
colnames(abs.metadata)[colnames(abs.metadata)=="Flare"]<-"Dx2" # change Flare column name to Dx2

# create group color palette to merge with abs.metadata
grp.clrs2 = as.data.frame(t(data.frame("HC"="#40557B","Flare"="#ffa600","Remission"="#ff6361","HHP"="#007561")))

grp.clrs2$Dx2<-rownames(grp.clrs2)
colnames(grp.clrs2)[which(names(grp.clrs2) == "V1")] <- "Group_Color"
grp.clrs2

abs.metadata<-merge(abs.metadata,grp.clrs2,by="Dx2")
unique(abs.metadata$Group_Color)

# # now import ABS-HC patient metadata
# abs.HC.metadata<-as.data.frame(read_xlsx("data/Metagenomes/ABS_Humann3_Results/revisions_output_merged_relab/seqIDmatched_clean_modified.xlsx", sheet="seqIDmatched_clean_nodups"))
# head(abs.HC.metadata)
# # ^^ match these data with the abs.HC.metadata$Combo_ID column & metaphlan table
# 
# # give ABS-HC a group color for later
# head(abs.HC.metadata)
# abs.HC.metadata$Group_Color<-"#40557B"
# 
# # compare the data you're about to rbind.fill()
# head(abs.metadata)
# head(abs.HC.metadata)
# 
# # created combined metadata ABS + HCs
# ## only including relevant information like Dx, Dx2, SampleID, Patient, Group_Color, etc
# abs.metadata<-rbind.fill(abs.metadata,abs.HC.metadata[,c(1:4,(ncol(abs.HC.metadata)-1):ncol(abs.HC.metadata))])
# 
# abs.metadata$Dx2 = factor(abs.metadata$Dx2, levels=c("Flare","Remission","HHP","ABS_HC"), ordered=TRUE)

#### Merge MetaPhlan RelAb and Metadata Together ####

# head(aud.mgm.tax.clean)
# head(aud.hc.mgm.tax.clean)
head(abs.mgm.tax.clean)

unique(abs.metadata$ID) %in% unique(abs.mgm.tax.clean$ID)
unique(abs.mgm.tax.clean$ID) %in% unique(abs.metadata$ID)

length(unique(abs.metadata$ID) %in% unique(abs.mgm.tax.clean$ID))

# drop all Eukaryotic hits from merged taxa data
"Eukaryota" %in% abs.mgm.tax.clean$Kingdom
abs.mgm.tax.clean<-abs.mgm.tax.clean[!abs.mgm.tax.clean$Kingdom %in% "Eukaryota",]
"Eukaryota" %in% abs.mgm.tax.clean$Kingdom

unique(abs.metadata$ID) %in% unique(abs.mgm.tax.clean$ID) # only Donor not included in Metaphlan results
unique(abs.metadata$ID)[which(!unique(abs.metadata$ID) %in% unique(abs.mgm.tax.clean$ID))]

unique(abs.mgm.tax.clean$ID)[which(!unique(abs.mgm.tax.clean$ID) %in% unique(abs.metadata$ID))]

length(unique(abs.metadata$ID))
## DATA NOTE: Humann3 data we used all data, including duplicate samples from individuals for ABS-HC data
## for MetaPhlan, Cynthia took the average of the Metaphlan results for each patient and made that one sample
## this is why there are 167 samples in the Humann3 analyses, but only 154 here

which(!unique(abs.metadata$ID) %in% unique(abs.mgm.tax.clean$ID))

# merge metadata and pathway relative abundance data together here
abs.tax.meta<-merge(abs.mgm.tax.clean,abs.metadata,by="ID") # 
head(abs.tax.meta)
length(unique(abs.tax.meta$ID)) # should be same as metadata minus the Donor sample

abs.tax.SampID<-unique(abs.tax.meta[,names(abs.tax.meta) %in% c("Kingdom", "Phylum", "Class", "Order", "Family",
                                                                "Genus","Species","Clade","relative_abundance","SampleID")])
# make sure that we did not lose any samples in this step
length(unique(abs.tax.meta$SampleID))
length(unique(abs.tax.meta$ID))
length(unique(abs.tax.SampID$SampleID))

head(abs.tax.SampID)
# ** will use abs.tax.SampID to create other dfs so we can keep track of samples with legible IDs rather than sequencing IDs

#### Create Kingdom x Relative Abundance Table ####

# Filter out Genus Relative Abundance

# drop empty Phylums AND keep empty Class, Order, Family, Genera, Species and Clade columns -- prevents repeating relative abundance values!
king.only.clean <- abs.tax.SampID[ which(abs.tax.SampID$Kingdom!="" & abs.tax.SampID$Phylum=="" & abs.tax.SampID$Class=="" & abs.tax.SampID$Order=="" & abs.tax.SampID$Family=="" & abs.tax.SampID$Genus=="" & abs.tax.SampID$Species=="" & abs.tax.SampID$Clade=="") , ]
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
phy.spec.typ<-subset(abs.tax.SampID,select=c(Phylum,Class,Order,Family,Genus,Species,Clade,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...
head(phy.spec.typ)

# drop empty Phylums AND keep empty Class, Order, Family, Genera, Species and Clade columns -- prevents repeating relative abundance values!
phy.only.clean <- phy.spec.typ[ which(phy.spec.typ$Phylum!="" & phy.spec.typ$Class=="" & phy.spec.typ$Order=="" & phy.spec.typ$Family=="" & phy.spec.typ$Genus=="" & phy.spec.typ$Species=="" & phy.spec.typ$Clade=="") , ]
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
cls.spec.typ<-subset(abs.tax.SampID,select=c(Class,Order,Family,Genus,Species,Clade,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...
head(cls.spec.typ)

# drop empty Classes AND keep empty Order, Family, Genera, Species and Clade columns -- prevents repeating relative abundance values!
cls.only.clean <- cls.spec.typ[ which(cls.spec.typ$Class!="" & cls.spec.typ$Order=="" & cls.spec.typ$Family=="" & cls.spec.typ$Genus=="" & cls.spec.typ$Species=="" & cls.spec.typ$Clade=="") , ]
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
ord.spec.typ<-subset(abs.tax.SampID,select=c(Order,Family,Genus,Species,Clade,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...
head(ord.spec.typ)

# drop empty Orders AND keep empty Family, Genera, Species and Clade columns -- prevents repeating relative abundance values!
ord.only.clean <- ord.spec.typ[ which(ord.spec.typ$Order!="" & ord.spec.typ$Family=="" & ord.spec.typ$Genus=="" & ord.spec.typ$Species=="" & ord.spec.typ$Clade=="") , ]
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
fam.spec.typ<-subset(abs.tax.SampID,select=c(Family,Genus,Species,Clade,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...
head(fam.spec.typ)

# drop empty Familys AND keep empty Family, Genera, Species and Clade columns -- prevents repeating relative abundance values!
fam.only.clean <- fam.spec.typ[ which(fam.spec.typ$Family!="" & fam.spec.typ$Genus=="" & fam.spec.typ$Species=="" & fam.spec.typ$Clade=="") , ]
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
gen.spec.typ<-subset(abs.tax.SampID,select=c(Genus,Species,Clade,relative_abundance,SampleID))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...

# drop empty Genera AND keep empty Species and Clade columns -- prevents repeating relative abundance values!
gen.only.clean <- gen.spec.typ[ which(gen.spec.typ$Genus!="" & gen.spec.typ$Species=="" & gen.spec.typ$Clade=="") , ]

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
spec.typ<-subset(abs.tax.SampID,select=c("Species","Clade","relative_abundance","SampleID"))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...

# drop empty Species AND keep empty Clade columns -- prevents repeating relative abundance values!
spec.only.clean <- spec.typ[ which(spec.typ$Species!="" & spec.typ$Clade=="") , ]

# keep species with relative abundance of >= 1%
#spec.only.filt<-spec.only.clean[spec.only.clean$relative_abundance>=1.0,]

# create JUST specus relative abundance table with newly filtered species relative abundance df
mgm.spec.relab<-as.data.frame(dcast(spec.only.clean,SampleID~Species,value.var = "relative_abundance",fun.aggregate=sum))
mgm.spec.relab[1:4,1:4]
rownames(mgm.spec.relab)<-mgm.spec.relab$SampleID
# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.spec.relab[,-1])

write.table(mgm.spec.relab,file="data/Metagenomes/ABS_Metaphlan_Results/Batch1to6andnewHC_metaphlan_abundance_table_March2025_SpeciesOnly.tsv",row.names=FALSE,col.names=TRUE,sep="\t")

#### Create Clade x Relative Abundance Table ####

# Subset out Species+Clade Relative Abundance
spec.strain.typ<-subset(abs.tax.SampID,select=c("Species","Clade","relative_abundance","SampleID"))
spec.strain.typ <- spec.strain.typ[ which(spec.strain.typ$Species!="" & spec.strain.typ$Clade!="") , ]
# ^ helps us identify which strain goes with which species

# create clade list (species + "strain") for later
clade.list<-unique(spec.strain.typ[,1:2])
clade.list$SpeciesClean<-gsub("_SGB[0-9].*","",clade.list$Species)
clade.list$FullClade<-interaction(clade.list$SpeciesClean,clade.list$Clade,sep="_")

# drop empty Species AND keep empty Clade columns -- prevents repeating relative abundance values!
#spec.only.clean <- spec.strain.typ[ which(spec.strain.typ$Species!="" & spec.strain.typ$Clade=="") , ]

# Filter out Clade Relative Abundance
strain.typ<-subset(abs.tax.SampID,select=c("Clade","relative_abundance","SampleID"))
## ^^ have to do this because relative abundance is calculated by clade, and mulitple values for rel.ab. appear if the genus is repeated...

# drop empty Species AND keep empty Clade columns -- prevents repeating relative abundance values!
strain.only.clean <- strain.typ[ which(strain.typ$Clade!="") , ]

# keep species with relative abundance of > 0
strain.only.filt<-strain.only.clean[strain.only.clean$relative_abundance>0,]

# create JUST species relative abundance table with newly filtered species relative abundance df
mgm.strain.relab<-as.data.frame(dcast(strain.only.clean,SampleID~Clade,value.var = "relative_abundance",fun.aggregate=sum))
mgm.strain.relab[1:4,1:4]
rownames(mgm.strain.relab)<-mgm.strain.relab$SampleID
# sanity check that the relative abundances add up to 100, and they do! or are close!!
rowSums(mgm.strain.relab[,-1])
dim(mgm.strain.relab[,-1])

write.table(mgm.strain.relab,file="data/Metagenomes/ABS_Metaphlan_Results/Batch1to6andnewHC_metaphlan_abundance_table_March2025_CladesOnly.tsv",row.names=FALSE,col.names=TRUE,sep="\t")


#### PCoA with Species Data ####
mgm.spec.relab[1:5,1:5]
rowSums(mgm.spec.relab[,-1])

# PCOA w/ Jaccard distance matrix (of binary data)
#pcoa(vegdist(abs.cooccur.table.sym,na.rm = TRUE,binary=TRUE,method="jaccard"),correction="lingoes") # 
abs.jac.pcoa1 <- pcoa(dist(mgm.spec.relab[,-1],method="binary"),correction="lingoes") # 

# The proportion of variances explained is in its element values$Relative_eig
abs.jac.pcoa1$values

# extract principal coordinates
abs.jac.pcoa1.vectors<-data.frame(abs.jac.pcoa1$vectors)
abs.jac.pcoa1.vectors$SampleID<-rownames(abs.jac.pcoa1$vectors)

# merge pcoa coordinates w/ abs.metadata
abs.jac.pcoa1.meta<-merge(abs.jac.pcoa1.vectors, abs.metadata, by.x="SampleID", by.y="SampleID")

head(abs.jac.pcoa1.meta)
rownames(abs.jac.pcoa1.meta)<-abs.jac.pcoa1.meta$newID

head(abs.jac.pcoa1$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 4.74%, PC2 = 2.15%

#### Homogeneity of Variance & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(abs.metadata)
head(mgm.spec.relab[,-1])

# are both dfs in the same order? Need this for stats below
abs.metadata$SampleID %in% rownames(mgm.spec.relab)

no.donor.metadata<-abs.metadata[abs.metadata$SampleID %in% rownames(mgm.spec.relab),]
no.donor.metadata<-unique(no.donor.metadata[,-c(2)])
rownames(no.donor.metadata)<-no.donor.metadata$SampleID

rownames(no.donor.metadata) %in% rownames(mgm.spec.relab)
no.donor.metadata=no.donor.metadata[rownames(mgm.spec.relab),] ## reorder metadata to match order of shared strain-clade data

# first by compare dispersions by sampling date
abs.disper1<-betadisper((dist(mgm.spec.relab[,-1],method="binary")), no.donor.metadata$Dx2)
abs.disper1

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(abs.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#           Flare      HC     HHP Remission
# Flare             0.54500 0.57100     0.566
# HC        0.54373         0.99000     0.964
# HHP       0.58338 0.98895             0.954
# Remission 0.56528 0.96370 0.95908          

anova(abs.disper1) # p = 0.9001 --> reject the Null H, spatial medians (a measure of dispersion) are significantly difference across conditions
# ANOVA adjusted p-value
aov.beta.p1<-anova(abs.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(abs.disper1) # tells us which category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj
# HC-Flare         0.0165996598 -0.05616096 0.08936028 0.9331992
# HHP-Flare        0.0169563240 -0.05192828 0.08584093 0.9178917
# Remission-Flare  0.0156107884 -0.05047733 0.08169890 0.9265585
# HHP-HC           0.0003566642 -0.07299442 0.07370775 0.9999993
# Remission-HC    -0.0009888714 -0.07172028 0.06974253 0.9999824
# Remission-HHP   -0.0013455356 -0.06808318 0.06539211 0.9999473

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself
## adonis2 will fail if there are samples with rowSums = 0, let's drop them
#mgm.spec.relab[is.na(mgm.spec.relab)]<-0
mgm.spec.relab2<-subset(mgm.spec.relab[,-1],rowSums(mgm.spec.relab[,-1])!=0)
no.donor.metadata2<-subset(no.donor.metadata,rownames(no.donor.metadata) %in% rownames(mgm.spec.relab2))
# went from 129 --> 97 samples, meaning 32 samples didn't share strains with any other samples

mgm.spec.PA.table<-ifelse(mgm.spec.relab[,-1]>0,1,0)
mgm.spec.PA.table[1:5,1:5]
mgm.spec.relab[1:5,1:6]

# use the pvalue below for PCOA that excludes Donor
pnova1<-adonis2(mgm.spec.PA.table ~ Dx2,data=no.donor.metadata,by="terms",method="jaccard",permutations= 10000)
pnova1 # p-value = 9.999e-05
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 0.00029997

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

clade.jac.dist = dist(mgm.spec.relab[,-1],method="binary") #Jaccard distance matrix
pair.mod1<-pairwise.adonis(clade.jac.dist,no.donor.metadata$Dx2, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1
#        pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1        Flare vs HC  1  1.557065 3.966016 0.09681234   1e-04     0.0006  **
# 2       Flare vs HHP  1  1.254836 2.935974 0.06124782   2e-04     0.0012   *
# 3 Flare vs Remission  1  1.622496 3.781513 0.07302824   1e-04     0.0006  **
# 4          HC vs HHP  1  1.500746 3.899418 0.09306616   1e-04     0.0006  **
# 5    HC vs Remission  1  1.346660 3.453688 0.07769181   1e-04     0.0006  **
# 6   HHP vs Remission  1  1.240358 2.937042 0.05655005   1e-04     0.0006  **

# Visualize dispersions
par(mar=c(1,1,1,1))
png('figures/SameSTR_Figs/ABS_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(abs.disper1,main = "Centroids and Dispersion based on Jaccard Distance")
dev.off()

par(mar=c(1,1,1,1))
png('figures/SameSTR_Figs/ABS_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(abs.disper1,xlab="By GroupName", main = "Distance to Centroid by Category", sub="Based on Jaccard Distance", col=grp.clrs$GroupName_Color)
dev.off()


#### Visualize PCoA ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints
head(abs.jac.pcoa1$values)
# create PCoA ggplot fig

pcoa1<-ggplot(abs.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Dx2)) +geom_point(size=4)+theme_bw()+
  labs(color="Flare",title="PCoA of Presence/Absence of Shared Clades",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
        axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  scale_fill_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  xlab("PC1 [11.70%]") + ylab("PC2 [8.97%]")+ annotate("text",x=0.25,y=0.5,label="PERMANOVA, p.adj = 0.0003")

ggsave(pcoa1,filename = "figures/RelativeAbundance/PCOA/ABS_MGMs_Species_Clade_Binary_PCOA1.png", 
       width = 2500,height = 1500,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa2<-ggplot(abs.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Dx2)) +geom_point(size=4)+theme_bw()+
  labs(color="Flare",title="PCoA of Presence/Absence of Shared Clades",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
        axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  scale_fill_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  xlab("PC1 [11.70%]") + ylab("PC2 [8.97%]")+ stat_ellipse(alpha = 0.5) + annotate("text",x=0.80,y=0.55,label="PERMANOVA, p.adj = 0.0003")

ggsave(pcoa2,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedStrains_Binary_Ellipses_PCOA2.png", 
       width = 2500,height = 1500,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa3<-ggplot(abs.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Dx2)) +geom_text(label=abs.jac.pcoa1.meta$SampleID,position="identity")+theme_bw()+
  labs(color="Flare",title="PCoA of Presence/Absence of Shared Clades",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
        axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  scale_fill_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  xlab("PC1 [11.70%]") + ylab("PC2 [8.97%]")+ annotate("text",x=0.40,y=0.5,label="PERMANOVA, p.adj = 0.0003")

ggsave(pcoa3,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedStrains_Binary_Labeled_by_SampleID_PCOA3.png", 
       width = 2500,height = 1500,units = "px",
       dpi = 200,create.dir = TRUE)

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


#### Arcsin(SqRt()) Transformation of CLADE Relative Abundance Data ####

# now we are going to use the Arcsin(Square Root()) Transformation to transform the pathway relative abundances
## sources for this transformation: LLorens-Rico et al 2021
## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)

mgm.strain.relab[1:4,1:4] # strains relative abundances
max(mgm.strain.relab[,-1]) # make sure that when we divide by 100, our max # is < 1 before transformation!
max(mgm.strain.relab[,-1]/100)

strain.trans.relab<-mgm.strain.relab #creating new data frame for transformation so that we can keep our SampleIDs & order of df
strain.trans.relab[1:4,1:4] # sanity check

# now for the transformation! 
## first we divide all values by 100 so we are working with scaled down #s < 1
## then we do the arcsin(square root()) transformation
strain.trans.relab[,-1]<-asin(sqrt(strain.trans.relab[,-1]/100)) # Arcsin(sqrt()) transformation of the taxa rel abs

strain.trans.relab[1:5,1:5]

# converting any NaN to 0
strain.trans.relab[is.na(strain.trans.relab)] # are there NAs?
#strain.trans.relab[is.na(strain.trans.relab)]<-0




#### Visualize Kingdom (Raw RelAb) by Dx2 ####

# first merge raw kingdom relative abundances + metadata
mgm.king.relab[1:5,]
mgm.king.relab.m<-melt(mgm.king.relab,by=vars("SampleID"))
head(mgm.king.relab.m)
colnames(mgm.king.relab.m)[colnames(mgm.king.relab.m)=="variable"]<-"Kingdom"
colnames(mgm.king.relab.m)[colnames(mgm.king.relab.m)=="value"]<-"RelAb"

head(abs.metadata)

abs.king.RA.meta<-merge(mgm.king.relab.m[mgm.king.relab.m$RelAb>0,],abs.metadata,by="SampleID")
abs.king.RA.meta[1:3,]
abs.king.RA.meta <- abs.king.RA.meta[order(abs.king.RA.meta$RelAb, decreasing = TRUE),]

# visualize the relab data
rawking.bplt1<-ggplot(abs.king.RA.meta, aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Kingdom)) +
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

ggsave(rawking.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Kingdom/RelativeAbundance_MGM_Metaphlan_Phyla_RawRA_barplot.png", width = 7000,height = 3000,units = "px",dpi = 200,create.dir = TRUE)

#### Visualize Phyla by Dx2 ####

# first merge raw phyla relative abundances + metadata
mgm.phy.relab[1:5,1:5]
mgm.phy.relab.m<-melt(mgm.phy.relab,by=vars("SampleID"))
head(mgm.phy.relab.m)
colnames(mgm.phy.relab.m)[colnames(mgm.phy.relab.m)=="variable"]<-"Phylum"
colnames(mgm.phy.relab.m)[colnames(mgm.phy.relab.m)=="value"]<-"RelAb"

head(abs.metadata)

abs.phy.RA.meta<-merge(mgm.phy.relab.m[mgm.phy.relab.m$RelAb>0,],abs.metadata,by="SampleID")
abs.phy.RA.meta[1:3,]
abs.phy.RA.meta <- abs.phy.RA.meta[order(abs.phy.RA.meta$RelAb, decreasing = TRUE),]

# visualize the relab data
rawphy.bplt1<-ggplot(abs.phy.RA.meta, aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Phylum)) +
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

ggsave(rawphy.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Phylum/RelativeAbundance_MGM_Metaphlan_Phyla_RawRA_barplot.png", width = 7000,height = 3000,units = "px",dpi = 200,create.dir = TRUE)

rawphy.bplt2<-ggplot(abs.phy.RA.meta, aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Phylum)) +
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

ggsave(rawphy.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Phylum/RelativeAbundance_MGM_Metaphlan_Phyla_RawRA_barplot_b.png", width = 7000,height = 3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
rawphy.ts1<-ggplot(abs.phy.RA.meta, aes(Phylum, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.phy.RA.meta$Group_Color[order(abs.phy.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Phylum", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Phyla")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawphy.ts1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Phylum/RelativeAbundance_MGM_Metaphlan_Phylum_RawRA_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, phyla relative abundances + metadata
phy.trans.relab[1:5,1:5]
phy.trans.relab.m<-melt(phy.trans.relab,by=vars("SampleID"))
head(phy.trans.relab.m)
colnames(phy.trans.relab.m)[colnames(phy.trans.relab.m)=="variable"]<-"Phylum"
colnames(phy.trans.relab.m)[colnames(phy.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(abs.metadata)

abs.phy.tRA.meta<-merge(phy.trans.relab.m[phy.trans.relab.m$Transformed_RelAb>0,],abs.metadata,by="SampleID")
abs.phy.tRA.meta[1:3,]
abs.phy.tRA.meta <- abs.phy.tRA.meta[order(abs.phy.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what is the max of the transformed relab?
max(abs.phy.tRA.meta$Transformed_RelAb)

# visualize transformed, relab data
trphy.bplt1<-ggplot(abs.phy.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Phylum)) +
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

ggsave(trphy.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Phylum/RelativeAbundance_MGM_Metaphlan_Phyla_AsinSqrtRA_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trphy.bplt2<-ggplot(abs.phy.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Phylum)) +
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

ggsave(trphy.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Phylum/RelativeAbundance_MGM_Metaphlan_Phyla_AsinSqrtRA_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trphy.bplt3<-ggplot(abs.phy.tRA.meta[abs.phy.tRA.meta$Transformed_RelAb>0.5,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Phylum)) +
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

ggsave(trphy.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Phylum/RelativeAbundance_MGM_Metaphlan_Phyla_AsinSqrtRA_Over0.5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)


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

head(abs.metadata)

abs.cls.RA.meta<-merge(mgm.cls.relab.m[mgm.cls.relab.m$RelAb>0,],abs.metadata,by="SampleID")
abs.cls.RA.meta[1:3,]
abs.cls.RA.meta <- abs.cls.RA.meta[order(abs.cls.RA.meta$RelAb, decreasing = TRUE),]

# what are the min, max of the raw RelAbs?
min(abs.cls.RA.meta$RelAb)

# how many colors do you need in the palette?
unique(abs.cls.RA.meta$Class[abs.cls.RA.meta$RelAb>5])

# visualize the relab data
rawcls.bplt1<-ggplot(abs.cls.RA.meta[abs.cls.RA.meta$RelAb>5,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Class)) +
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

ggsave(rawcls.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Class/RelativeAbundance_MGM_Metaphlan_Class_RawRA_Over5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# how many colors do you need in the palette?
unique(abs.cls.RA.meta$Class[abs.cls.RA.meta$RelAb>20])

rawcls.bplt2<-ggplot(abs.cls.RA.meta[abs.cls.RA.meta$RelAb>20,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Class)) +
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

ggsave(rawcls.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Class/RelativeAbundance_MGM_Metaphlan_Class_RawRA_Over20_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
rawcls.ts1<-ggplot(abs.cls.RA.meta[abs.cls.RA.meta$RelAb>5,], aes(Class, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.cls.RA.meta$Group_Color[order(abs.cls.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle = element_text(size=35)) +
  labs(x="Class", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Classes",subtitle="Includes taxa with Relative Abundance > 5%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawcls.ts1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Class/RelativeAbundance_MGM_Metaphlan_Class_RawRA_Over5_taxo.sum.png", width = 5000,height=7000,units = "px",dpi = 200,create.dir = TRUE)

rawcls.ts2<-ggplot(abs.cls.RA.meta[abs.cls.RA.meta$RelAb>10,], aes(Class, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.cls.RA.meta$Group_Color[order(abs.cls.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),plot.subtitle = element_text(size=35)) +
  labs(x="Class", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Classes",subtitle="Includes taxa with Relative Abundance > 10%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawcls.ts2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Class/RelativeAbundance_MGM_Metaphlan_Class_RawRA_Over10_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, class relative abundances + metadata
cls.trans.relab[1:5,1:5]
cls.trans.relab.m<-melt(cls.trans.relab,by=vars("SampleID"))
head(cls.trans.relab.m)
colnames(cls.trans.relab.m)[colnames(cls.trans.relab.m)=="variable"]<-"Class"
colnames(cls.trans.relab.m)[colnames(cls.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(abs.metadata)

abs.cls.tRA.meta<-merge(cls.trans.relab.m[cls.trans.relab.m$Transformed_RelAb>0,],abs.metadata,by="SampleID")
abs.cls.tRA.meta[1:3,]
abs.cls.tRA.meta <- abs.cls.tRA.meta[order(abs.cls.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what are the mins, maxs of transformed rel ab data?
max(abs.cls.tRA.meta$Transformed_RelAb)
min(abs.cls.tRA.meta$Transformed_RelAb)

# how many colors do you need in the palette?
unique(abs.cls.tRA.meta$Class[abs.cls.tRA.meta$Transformed_RelAb>0.3])

# visualize transformed, relab data
ggplot(abs.cls.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
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
unique(abs.cls.tRA.meta$Class[abs.cls.tRA.meta$Transformed_RelAb>0.3])

trcls.bplt1<-ggplot(abs.cls.tRA.meta[abs.cls.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
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

ggsave(trcls.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Class/RelativeAbundance_MGM_Metaphlan_Class_AsinSqrtRA_Over0.3_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trcls.bplt2<-ggplot(abs.cls.tRA.meta[abs.cls.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
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

ggsave(trcls.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Class/RelativeAbundance_MGM_Metaphlan_Class_AsinSqrtRA_Over0.3_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# how many colors do you need in the palette?
unique(abs.cls.tRA.meta$Class[abs.cls.tRA.meta$Transformed_RelAb>0.6])

trcls.bplt3<-ggplot(abs.cls.tRA.meta[abs.cls.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
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

ggsave(trcls.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Class/RelativeAbundance_MGM_Metaphlan_Class_AsinSqrtRA_Over0.6_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trcls.bplt4<-ggplot(abs.cls.tRA.meta[abs.cls.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Class)) +
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

ggsave(trcls.bplt4,filename = "figures/RelativeAbundance/Taxa_RelAbs/Class/RelativeAbundance_MGM_Metaphlan_Class_AsinSqrtRA_Over0.6_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

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

head(abs.metadata)

abs.ord.RA.meta<-merge(mgm.ord.relab.m[mgm.ord.relab.m$RelAb>0,],abs.metadata,by="SampleID")
abs.ord.RA.meta[1:3,]
abs.ord.RA.meta <- abs.ord.RA.meta[order(abs.ord.RA.meta$RelAb, decreasing = TRUE),]

# what are the min, max of the raw RelAbs?
min(abs.ord.RA.meta$RelAb)
max(abs.ord.RA.meta$RelAb)

# how many colors for palette?
unique(abs.ord.RA.meta$Order)

# visualize the relab data
raword.bplt1<-ggplot(abs.ord.RA.meta[abs.ord.RA.meta$RelAb>5,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Order)) +
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

ggsave(raword.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Order/RelativeAbundance_MGM_Metaphlan_Order_RawRA_Over5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

raword.bplt2<-ggplot(abs.ord.RA.meta[abs.ord.RA.meta$RelAb>10,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Order)) +
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

ggsave(raword.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Order/RelativeAbundance_MGM_Metaphlan_Order_RawRA_Over10_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

raword.bplt3<-ggplot(abs.ord.RA.meta[abs.ord.RA.meta$RelAb>20,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Order)) +
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

ggsave(raword.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Order/RelativeAbundance_MGM_Metaphlan_Order_RawRA_Over20_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
raword.ts1<-ggplot(abs.ord.RA.meta[abs.ord.RA.meta$RelAb>5,], aes(Order, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.ord.RA.meta$Group_Color[order(abs.ord.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, legend.title = element_text(size=40),
        legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle = element_text(size=35)) +
  labs(x="Order", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Orders",subtitle="Includes taxa with Relative Abundance > 5%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(raword.ts1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Order/RelativeAbundance_MGM_Metaphlan_Order_RawRA_Over5_taxo.sum.png", width = 5000,height=7000,units = "px",dpi = 200,create.dir = TRUE)

raword.ts2<-ggplot(abs.ord.RA.meta[abs.ord.RA.meta$RelAb>10,], aes(Order, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.ord.RA.meta$Group_Color[order(abs.ord.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, legend.title = element_text(size=40),
        legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle = element_text(size=35)) +
  labs(x="Order", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Orders",subtitle="Includes taxa with Relative Abundance > 10%")+
  coord_flip()

ggsave(raword.ts2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Order/RelativeAbundance_MGM_Metaphlan_Order_RawRA_Over10_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, class relative abundances + metadata
ord.trans.relab[1:5,1:5]
ord.trans.relab.m<-melt(ord.trans.relab,by=vars("SampleID"))
head(ord.trans.relab.m)
colnames(ord.trans.relab.m)[colnames(ord.trans.relab.m)=="variable"]<-"Order"
colnames(ord.trans.relab.m)[colnames(ord.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(abs.metadata)

abs.ord.tRA.meta<-merge(ord.trans.relab.m[ord.trans.relab.m$Transformed_RelAb>0,],abs.metadata,by="SampleID")
abs.ord.tRA.meta[1:3,]
abs.ord.tRA.meta <- abs.ord.tRA.meta[order(abs.ord.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what are the min, max of transformed order relabs?
max(abs.ord.tRA.meta$Transformed_RelAb)
min(abs.ord.tRA.meta$Transformed_RelAb)

# visualize transformed, relab data
ggplot(abs.ord.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
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

trord.bplt1<-ggplot(abs.ord.tRA.meta[abs.ord.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
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

ggsave(trord.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Order/RelativeAbundance_MGM_Metaphlan_Order_AsinSqrtRA_Over0.3_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trord.bplt2<-ggplot(abs.ord.tRA.meta[abs.ord.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
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

ggsave(trord.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Order/RelativeAbundance_MGM_Metaphlan_Order_AsinSqrtRA_Over0.3_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trord.bplt3<-ggplot(abs.ord.tRA.meta[abs.ord.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
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

ggsave(trord.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Order/RelativeAbundance_MGM_Metaphlan_Order_AsinSqrtRA_Over0.6_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trord.bplt4<-ggplot(abs.ord.tRA.meta[abs.ord.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Order)) +
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

ggsave(trord.bplt4,filename = "figures/RelativeAbundance/Taxa_RelAbs/Order/RelativeAbundance_MGM_Metaphlan_Order_AsinSqrtRA_Over0.6_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

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

head(abs.metadata)

abs.fam.RA.meta<-merge(mgm.fam.relab.m[mgm.fam.relab.m$RelAb>0,],abs.metadata,by="SampleID")
abs.fam.RA.meta[1:3,]
abs.fam.RA.meta <- abs.fam.RA.meta[order(abs.fam.RA.meta$RelAb, decreasing = TRUE),]

# what are the min, max of the raw RelAbs?
min(abs.fam.RA.meta$RelAb)
max(abs.fam.RA.meta$RelAb)

# how many colors for palette?
unique(abs.fam.RA.meta$Family)

# visualize the relab data
rawfam.bplt1<-ggplot(abs.fam.RA.meta[abs.fam.RA.meta$RelAb>5,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Family)) +
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

ggsave(rawfam.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Family/RelativeAbundance_MGM_Metaphlan_Family_RawRA_Over5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawfam.bplt2<-ggplot(abs.fam.RA.meta[abs.fam.RA.meta$RelAb>10,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Family)) +
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

ggsave(rawfam.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Family/RelativeAbundance_MGM_Metaphlan_Family_RawRA_Over10_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawfam.bplt3<-ggplot(abs.fam.RA.meta[abs.fam.RA.meta$RelAb>20,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Family)) +
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

ggsave(rawfam.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Family/RelativeAbundance_MGM_Metaphlan_Family_RawRA_Over20_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
rawfam.ts1<-ggplot(abs.fam.RA.meta[abs.fam.RA.meta$RelAb>5,], aes(Family, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.fam.RA.meta$Group_Color[order(abs.fam.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40,),
        axis.text = element_text(size=35),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Family", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Families",subtitle="Includes taxa with Relative Abundance > 5%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawfam.ts1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Family/RelativeAbundance_MGM_Metaphlan_Family_RawRA_Over10_taxo.sum.png", width = 5000,height=7000,units = "px",dpi = 200,create.dir = TRUE)

rawfam.ts2<-ggplot(abs.fam.RA.meta[abs.fam.RA.meta$RelAb>10,], aes(Family, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.fam.RA.meta$Group_Color[order(abs.fam.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle = element_text(size=35)) +
  labs(x="Family", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Families",subtitle="Includes taxa with Relative Abundance > 10%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawfam.ts2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Family/RelativeAbundance_MGM_Metaphlan_Family_RawRA_Over10_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, class relative abundances + metadata
fam.trans.relab[1:5,1:5]
fam.trans.relab.m<-melt(fam.trans.relab,by=vars("SampleID"))
head(fam.trans.relab.m)
colnames(fam.trans.relab.m)[colnames(fam.trans.relab.m)=="variable"]<-"Family"
colnames(fam.trans.relab.m)[colnames(fam.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(abs.metadata)

abs.fam.tRA.meta<-merge(fam.trans.relab.m[fam.trans.relab.m$Transformed_RelAb>0,],abs.metadata,by="SampleID")
abs.fam.tRA.meta[1:3,]
abs.fam.tRA.meta <- abs.fam.tRA.meta[order(abs.fam.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what are the min, max of transformed order relabs?
max(abs.fam.tRA.meta$Transformed_RelAb)
min(abs.fam.tRA.meta$Transformed_RelAb)

# visualize transformed, relab data
ggplot(abs.fam.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
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

trfam.bplt1<-ggplot(abs.fam.tRA.meta[abs.fam.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
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

ggsave(trfam.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Family/RelativeAbundance_MGM_Metaphlan_Family_AsinSqrtRA_Over0.3_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trfam.bplt2<-ggplot(abs.fam.tRA.meta[abs.fam.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
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

ggsave(trfam.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Family/RelativeAbundance_MGM_Metaphlan_Family_AsinSqrtRA_Over0.3_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trfam.bplt3<-ggplot(abs.fam.tRA.meta[abs.fam.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
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

ggsave(trfam.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Family/RelativeAbundance_MGM_Metaphlan_Family_AsinSqrtRA_Over0.6_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trfam.bplt4<-ggplot(abs.fam.tRA.meta[abs.fam.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Family)) +
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

ggsave(trfam.bplt4,filename = "figures/RelativeAbundance/Taxa_RelAbs/Family/RelativeAbundance_MGM_Metaphlan_Family_AsinSqrtRA_Over0.6_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)


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

head(abs.metadata)

abs.gen.RA.meta<-merge(mgm.gen.relab.m[mgm.gen.relab.m$RelAb>0,],abs.metadata,by="SampleID")
abs.gen.RA.meta[1:3,]
abs.gen.RA.meta <- abs.gen.RA.meta[order(abs.gen.RA.meta$RelAb, decreasing = TRUE),]

# what are the min, max of the raw RelAbs?
min(abs.gen.RA.meta$RelAb)
max(abs.gen.RA.meta$RelAb)

# how many colors for palette?
unique(abs.gen.RA.meta$Genus)

# visualize the relab data
rawgen.bplt1<-ggplot(abs.gen.RA.meta[abs.gen.RA.meta$RelAb>5,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Genus)) +
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

ggsave(rawgen.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_RawRA_Over5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawgen.bplt2<-ggplot(abs.gen.RA.meta[abs.gen.RA.meta$RelAb>10,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Genus)) +
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

ggsave(rawgen.bplt2,filename = "figures/ABS/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_RawRA_Over10_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawgen.bplt3<-ggplot(abs.gen.RA.meta[abs.gen.RA.meta$RelAb>20,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Genus)) +
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

ggsave(rawgen.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_RawRA_Over20_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawgen.bplt3<-ggplot(abs.gen.RA.meta[abs.gen.RA.meta$RelAb>30,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Genus)) +
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

ggsave(rawgen.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_RawRA_Over30_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
rawgen.ts1<-ggplot(abs.gen.RA.meta[abs.gen.RA.meta$RelAb>5,], aes(Genus, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.gen.RA.meta$Group_Color[order(abs.gen.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle=element_text(size=35)) +
  labs(x="Genera", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Genera",subtitle="Includes taxa with Relative Abundance > 5%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawgen.ts1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_RawRA_Over5_taxo.sum.png", width = 5000,height=7000,units = "px",dpi = 200,create.dir = TRUE)

rawgen.ts2<-ggplot(abs.gen.RA.meta[abs.gen.RA.meta$RelAb>10,], aes(Genus, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.gen.RA.meta$Group_Color[order(abs.gen.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle=element_text(size=35)) +
  labs(x="Genera", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Genera",subtitle="Includes taxa with Relative Abundance > 10%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawgen.ts2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_RawRA_Over10_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, class relative abundances + metadata
gen.trans.relab[1:5,1:5]
gen.trans.relab.m<-melt(gen.trans.relab,by=vars("SampleID"))
head(gen.trans.relab.m)
colnames(gen.trans.relab.m)[colnames(gen.trans.relab.m)=="variable"]<-"Genus"
colnames(gen.trans.relab.m)[colnames(gen.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(abs.metadata)

abs.gen.tRA.meta<-merge(gen.trans.relab.m[gen.trans.relab.m$Transformed_RelAb>0,],abs.metadata,by="SampleID")
abs.gen.tRA.meta[1:3,]
abs.gen.tRA.meta <- abs.gen.tRA.meta[order(abs.gen.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what are the min, max of transformed order relabs?
max(abs.gen.tRA.meta$Transformed_RelAb)
min(abs.gen.tRA.meta$Transformed_RelAb)

# visualize transformed, relab data
ggplot(abs.gen.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
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

trgen.bplt1<-ggplot(abs.gen.tRA.meta[abs.gen.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
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

ggsave(trgen.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_AsinSqrtRA_Over0.3_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trgen.bplt2<-ggplot(abs.gen.tRA.meta[abs.gen.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
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

ggsave(trgen.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_AsinSqrtRA_Over0.3_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trgen.bplt3<-ggplot(abs.gen.tRA.meta[abs.gen.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
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

ggsave(trgen.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_AsinSqrtRA_Over0.6_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trgen.bplt4<-ggplot(abs.gen.tRA.meta[abs.gen.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Genus)) +
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

ggsave(trgen.bplt4,filename = "figures/RelativeAbundance/Taxa_RelAbs/Genus/RelativeAbundance_MGM_Metaphlan_Genus_AsinSqrtRA_Over0.6_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)



#### Visualize Clade by Dx2 ####

# first merge raw class relative abundances + metadata
mgm.strain.relab[1:5,1:5]
mgm.strain.relab.m<-melt(mgm.strain.relab,by=vars("SampleID"))
head(mgm.strain.relab.m)
colnames(mgm.strain.relab.m)[colnames(mgm.strain.relab.m)=="variable"]<-"Clade"
colnames(mgm.strain.relab.m)[colnames(mgm.strain.relab.m)=="value"]<-"RelAb"

# convert Humann3 specific labels to Other
unique(mgm.strain.relab.m$Clade) # to see what to change
#mgm.strain.relab.m$Clade<-gsub("GGB.*","Other",mgm.strain.relab.m$Clade)
# convert any unclassified taxa as just Unclassified
mgm.strain.relab.m$Clade<-gsub(".*_unclassified","Unclassified",mgm.strain.relab.m$Clade)

head(abs.metadata)

abs.strain.RA.meta<-merge(mgm.strain.relab.m[mgm.strain.relab.m$RelAb>0,],abs.metadata,by="SampleID")
abs.strain.RA.meta[1:3,]
abs.strain.RA.meta <- abs.strain.RA.meta[order(abs.strain.RA.meta$RelAb, decreasing = TRUE),]

# what are the min, max of the raw RelAbs?
min(abs.strain.RA.meta$RelAb)
max(abs.strain.RA.meta$RelAb)

# how many colors for palette?
unique(abs.strain.RA.meta$Clade)

# visualize the relab data
rawstrain.bplt1<-ggplot(abs.strain.RA.meta[abs.strain.RA.meta$RelAb>5,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Clade)) +
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
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Relative Abundance",fill="Clade",
       subtitle="Includes Taxa with Relative Abundance > 5%") +
  guides(fill=guide_legend(ncol=4))

ggsave(rawstrain.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_RawRA_Over5_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawstrain.bplt2<-ggplot(abs.strain.RA.meta[abs.strain.RA.meta$RelAb>10,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Clade)) +
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
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Relative Abundance",fill="Clade",
       subtitle="Includes Taxa with Relative Abundance > 10%")

ggsave(rawstrain.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_RawRA_Over10_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawstrain.bplt3<-ggplot(abs.strain.RA.meta[abs.strain.RA.meta$RelAb>20,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Clade)) +
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
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Relative Abundance",fill="Clade",
       subtitle="Includes Taxa with Relative Abundance > 20%")

ggsave(rawstrain.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_RawRA_Over20_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

rawstrain.bplt3<-ggplot(abs.strain.RA.meta[abs.strain.RA.meta$RelAb>30,], aes(x = reorder(SampleID,RelAb), y = RelAb, fill = Clade)) +
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
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Relative Abundance",fill="Clade",
       subtitle="Includes Taxa with Relative Abundance > 30%")

ggsave(rawstrain.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_RawRA_Over30_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

# taxonomic summary
rawstrain.ts1<-ggplot(abs.strain.RA.meta[abs.strain.RA.meta$RelAb>5,], aes(Clade, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.strain.RA.meta$Group_Color[order(abs.strain.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle=element_text(size=35)) +
  labs(x="Genera", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Genera",subtitle="Includes taxa with Relative Abundance > 5%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawstrain.ts1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_RawRA_Over5_taxo.sum.png", width = 5000,height=7000,units = "px",dpi = 200,create.dir = TRUE)

rawstrain.ts2<-ggplot(abs.strain.RA.meta[abs.strain.RA.meta$RelAb>10,], aes(Clade, RelAb)) +
  geom_jitter(aes(color=factor(Dx2)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Dx2",values=unique(abs.strain.RA.meta$Group_Color[order(abs.strain.RA.meta$Dx2)])) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=40),axis.title.y = element_text(size=40),
        axis.text = element_text(size=35),axis.text.y = element_text(face="italic"),legend.title.align=0.5, 
        legend.title = element_text(size=40),legend.text = element_text(size=35),plot.title = element_text(size=32),plot.margin = unit(c(1, 1, 1, 2),"cm"),
        plot.subtitle=element_text(size=35)) +
  labs(x="Genera", y="Relative Abundance", title="AUD vs ABS - MetaPhlan Bacterial Genera",subtitle="Includes taxa with Relative Abundance > 10%")+
  coord_flip() + scale_x_discrete(limits=rev)

ggsave(rawstrain.ts2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_RawRA_Over10_taxo.sum.png", width = 5000,height=6000,units = "px",dpi = 200,create.dir = TRUE)

# then merge transformed, class relative abundances + metadata
strain.trans.relab[1:5,1:5]
strain.trans.relab.m<-melt(strain.trans.relab,by=vars("SampleID"))
head(strain.trans.relab.m)
colnames(strain.trans.relab.m)[colnames(strain.trans.relab.m)=="variable"]<-"Clade"
colnames(strain.trans.relab.m)[colnames(strain.trans.relab.m)=="value"]<-"Transformed_RelAb"

head(abs.metadata)

abs.strain.tRA.meta<-merge(strain.trans.relab.m[strain.trans.relab.m$Transformed_RelAb>0,],abs.metadata,by="SampleID")
abs.strain.tRA.meta[1:3,]
abs.strain.tRA.meta <- abs.strain.tRA.meta[order(abs.strain.tRA.meta$Transformed_RelAb, decreasing = TRUE),]

# what are the min, max of transformed order relabs?
max(abs.strain.tRA.meta$Transformed_RelAb)
min(abs.strain.tRA.meta$Transformed_RelAb)

# visualize transformed, relab data
ggplot(abs.strain.tRA.meta, aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Clade)) +
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
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Clade")

trstrain.bplt1<-ggplot(abs.strain.tRA.meta[abs.strain.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Clade)) +
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
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Clade",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trstrain.bplt1,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_AsinSqrtRA_Over0.3_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trstrain.bplt2<-ggplot(abs.strain.tRA.meta[abs.strain.tRA.meta$Transformed_RelAb>0.3,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Clade)) +
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
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Clade",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.3")

ggsave(trstrain.bplt2,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_AsinSqrtRA_Over0.3_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trstrain.bplt3<-ggplot(abs.strain.tRA.meta[abs.strain.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Clade)) +
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
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Clade",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trstrain.bplt3,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_AsinSqrtRA_Over0.6_barplot.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)

trstrain.bplt4<-ggplot(abs.strain.tRA.meta[abs.strain.tRA.meta$Transformed_RelAb>0.6,], aes(x = reorder(SampleID,Transformed_RelAb), y = Transformed_RelAb, fill = Clade)) +
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
  labs(title = "AUD vs ABS - Bacterial Genera in Metagenomes", x="", y="Transformed Relative Abundance",fill="Clade",
       subtitle="Includes Taxa with Transformed Relative Abundance > 0.6")

ggsave(trstrain.bplt4,filename = "figures/RelativeAbundance/Taxa_RelAbs/Clade/RelativeAbundance_MGM_Metaphlan_Clade_AsinSqrtRA_Over0.6_barplot_b.png", width = 7000,height=3000,units = "px",dpi = 200,create.dir = TRUE)




#### PCoA with Clade Data ####
head(clade.list) # will use FullClade column for work below
mgm.strain.relab[1:5,1:5]
mgm.strain.relab.m[1:5,]

clade.relab.m<-merge(mgm.strain.relab.m,clade.list[,-1],by="Clade")
head(clade.relab.m)

clade.relab.table<-dcast(clade.relab.m,SampleID~FullClade,value.var = "RelAb")
clade.relab.table[1:5,1:5] # sanity check
rownames(clade.relab.table)<-clade.relab.table$SampleID
dim(clade.relab.table[,-1])

write.table(clade.relab.table,file="data/Metagenomes/ABS_Metaphlan_Results/Batch1to6andnewHC_metaphlan_abundance_table_March2025_SpeciesClades_Only.tsv",row.names=FALSE,col.names=TRUE,sep="\t")

rowSums(clade.relab.table[,-1])
clade.relab.table.2<-clade.relab.table[which(rowSums(clade.relab.table[,-1])>0),]

# PCOA w/ Jaccard distance matrix (of binary data)
#pcoa(vegdist(abs.cooccur.table.sym,na.rm = TRUE,binary=TRUE,method="jaccard"),correction="lingoes") # 
abs.jac.pcoa1 <- pcoa(dist(clade.relab.table2[,-1],method="binary"),correction="lingoes") # 

# The proportion of variances explained is in its element values$Relative_eig
abs.jac.pcoa1$values

# extract principal coordinates
abs.jac.pcoa1.vectors<-data.frame(abs.jac.pcoa1$vectors)
abs.jac.pcoa1.vectors$SampleID<-rownames(abs.jac.pcoa1$vectors)

# merge pcoa coordinates w/ abs.metadata
abs.jac.pcoa1.meta<-merge(abs.jac.pcoa1.vectors, abs.metadata, by.x="SampleID", by.y="SampleID")

head(abs.jac.pcoa1.meta)
rownames(abs.jac.pcoa1.meta)<-abs.jac.pcoa1.meta$newID

head(abs.jac.pcoa1$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 4.74%, PC2 = 2.15%

#### Homogeneity of Variance & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(abs.metadata)
head(clade.relab.table[,-1])

# are both dfs in the same order? Need this for stats below
rownames(abs.metadata) %in% rownames(clade.relab.table[,-1])
abs.metadata$SampleID %in% rownames(clade.relab.table)

no.donor.metadata<-abs.metadata[abs.metadata$SampleID %in% rownames(clade.relab.table),]
no.donor.metadata<-unique(no.donor.metadata[,-c(2)])
rownames(no.donor.metadata)<-no.donor.metadata$SampleID

rownames(no.donor.metadata) %in% rownames(clade.relab.table)
no.donor.metadata=no.donor.metadata[rownames(clade.relab.table),] ## reorder metadata to match order of shared strain-clade data

# first by compare dispersions by sampling date
abs.disper1<-betadisper((dist(clade.relab.table[,-1],method="binary")), no.donor.metadata$Dx2)
abs.disper1

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(abs.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#           Flare      HC     HHP Remission
# Flare                1.0000e-03 2.9400e-01     0.189
# HC        1.1629e-06            6.0000e-03     0.007
# HHP       2.8925e-01 2.2179e-03                0.852
# Remission 2.0058e-01 4.8861e-03 8.5302e-01          

anova(abs.disper1) # p = 0.001649 --> reject the Null H, spatial medians (a measure of dispersion) are significantly difference across conditions
# ANOVA adjusted p-value
aov.beta.p1<-anova(abs.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(abs.disper1) # tells us which category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj
# HC-Flare         0.42677887  0.1443847  0.70917306 0.0008204
# HHP-Flare        0.11640212 -0.1509487  0.38375296 0.6677853
# Remission-Flare  0.13839286 -0.1181044  0.39489012 0.4968675
# HHP-HC          -0.31037675 -0.5950626 -0.02569087 0.0269026
# Remission-HC    -0.28838601 -0.5629045 -0.01386748 0.0355608
# Remission-HHP    0.02199074 -0.2370274  0.28100892 0.9961350

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself
## adonis2 will fail if there are samples with rowSums = 0, let's drop them
#clade.relab.table[is.na(clade.relab.table)]<-0
clade.relab.table2<-subset(clade.relab.table[,-1],rowSums(clade.relab.table[,-1])!=0)
no.donor.metadata2<-subset(no.donor.metadata,rownames(no.donor.metadata) %in% rownames(clade.relab.table2))
# went from 129 --> 97 samples, meaning 32 samples didn't share strains with any other samples

clade.PA.final.table<-ifelse(clade.relab.table[,-1]>0,1,0)
# use the pvalue below for PCOA that excludes Donor
pnova1<-adonis2(dist(clade.relab.table2,method="binary") ~ Dx2,data=no.donor.metadata2,by="terms",permutations= 10000)
pnova1 # p-value = 9.999e-05
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 0.00029997

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

clade.jac.dist = dist(clade.relab.table[,-1],method="binary") #Jaccard distance matrix
pair.mod1<-pairwise.adonis(clade.jac.dist,no.donor.metadata$Dx2, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1
#        pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1        Flare vs HC  1  1.557065 3.966016 0.09681234   1e-04     0.0006  **
# 2       Flare vs HHP  1  1.254836 2.935974 0.06124782   2e-04     0.0012   *
# 3 Flare vs Remission  1  1.622496 3.781513 0.07302824   1e-04     0.0006  **
# 4          HC vs HHP  1  1.500746 3.899418 0.09306616   1e-04     0.0006  **
# 5    HC vs Remission  1  1.346660 3.453688 0.07769181   1e-04     0.0006  **
# 6   HHP vs Remission  1  1.240358 2.937042 0.05655005   1e-04     0.0006  **

# Visualize dispersions
par(mar=c(1,1,1,1))
png('figures/SameSTR_Figs/ABS_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(abs.disper1,main = "Centroids and Dispersion based on Jaccard Distance")
dev.off()

par(mar=c(1,1,1,1))
png('figures/SameSTR_Figs/ABS_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(abs.disper1,xlab="By GroupName", main = "Distance to Centroid by Category", sub="Based on Jaccard Distance", col=grp.clrs$GroupName_Color)
dev.off()


#### Visualize PCoA ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints
head(abs.jac.pcoa1$values)
# create PCoA ggplot fig

pcoa1<-ggplot(abs.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Dx2)) +geom_point(size=4)+theme_bw()+
  labs(color="Flare",title="PCoA of Presence/Absence of Shared Clades",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
        axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  scale_fill_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  xlab("PC1 [11.70%]") + ylab("PC2 [8.97%]")+ annotate("text",x=0.25,y=0.5,label="PERMANOVA, p.adj = 0.0003")

ggsave(pcoa1,filename = "figures/RelativeAbundance/PCOA/ABS_MGMs_Species_Clade_Binary_PCOA1.png", 
       width = 2500,height = 1500,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa2<-ggplot(abs.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Dx2)) +geom_point(size=4)+theme_bw()+
  labs(color="Flare",title="PCoA of Presence/Absence of Shared Clades",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
        axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  scale_fill_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  xlab("PC1 [11.70%]") + ylab("PC2 [8.97%]")+ stat_ellipse(alpha = 0.5) + annotate("text",x=0.80,y=0.55,label="PERMANOVA, p.adj = 0.0003")

ggsave(pcoa2,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedClades_Binary_Ellipses_PCOA2.png", 
       width = 2500,height = 1500,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa3<-ggplot(abs.jac.pcoa1.meta, aes(x=Axis.1, y=Axis.2,color=Dx2)) +geom_text(label=abs.jac.pcoa1.meta$SampleID,position="identity")+theme_bw()+
  labs(color="Flare",title="PCoA of Presence/Absence of Shared Clades",subtitle="Using Jaccard Distance (Presence/Absence Data)")+theme_classic()+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),axis.text = element_text(size=17),
        axis.text.x = element_text(vjust=1),legend.text = element_text(size=15),legend.title=element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  scale_fill_manual(name ="Flare",values=unique(abs.jac.pcoa1.meta$Group_Color[order(abs.jac.pcoa1.meta$Dx2)])) +
  xlab("PC1 [11.70%]") + ylab("PC2 [8.97%]")+ annotate("text",x=0.40,y=0.5,label="PERMANOVA, p.adj = 0.0003")

ggsave(pcoa3,filename = "figures/SameSTR_Figs/ABS_MGMs_SameSTR_SharedClades_Binary_Labeled_by_SampleID_PCOA3.png", 
       width = 2500,height = 1500,units = "px",
       dpi = 200,create.dir = TRUE)

#### Save Global Envionrment as RData Object ####
save.image("data/RelativeAbundance/RelativeAbundance_MetaPhlan_BacteriaRelAb_ArcSinTransformed_Ready.Rdata")
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

head(gen.spec.typ) # subsetted version of aud.mgm.tax.clean1 but only includes Genus, Species, Clade, relative_abundance, and SampleID columns
head(spec.typ) # subsetted version of aud.mgm.tax.clean1 but only includes Species, Clade, relative_abundance, and SampleID columns (no Genus!)

gen.only.clean[1:5,1:5] # bacterial species relative abundance

spec.only.clean[1:5,] # bacterial species relative abundance


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
})

#### Import and Prepare Data for Analyses ####

## Import ALL env plate fungal ASV count data
its2.ASV_counts<-data.frame(readRDS("data/ITS2_DADA2/ITS2_ASVs_Counts_dada2_Updated_Robject.rds", refhook = NULL))
dim(its2.ASV_counts)
its2.ASV_counts[1:4,1:4]

## Import ASV taxonomic data
its2.ASV_tax<-data.frame(readRDS("data/ITS2_DADA2/ITS2_ASVs_Taxonomy_dada2_Robject.rds", refhook = NULL))
head(its2.ASV_tax)

# Turn NAs into Unknowns
its2.ASV_tax[is.na(its2.ASV_tax)]<- "Unknown" # turn all NAs into "Unknowns"
its2.ASV_tax$Species<-gsub("Unknown", "unknown", its2.ASV_tax$Species) # change uppercase Unkonwn to lowercase unknown for unknown species classification
head(its2.ASV_tax)

# drop leading characters in front of taxonomic IDs
## keeping the index [,1:7] keeps the dataframe a dataframe and does not remove the rownames, which are our ASV IDs!
its2.ASV_tax[,1:7] <- sapply(its2.ASV_tax[,1:7],function(x) gsub("[:a-z:]__","",as.character(x)))
head(its2.ASV_tax)
its2.ASV_tax[,1:7] <- sapply(its2.ASV_tax[,1:7],function(x) gsub("_"," ",as.character(x)))

its2.ASV_tax[30:40,] # sanity check that this worked
#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata")

#### Import metadata ####

# Import all metadata
metadata<-as.data.frame(read.delim('data/ITS2fileName_metadata_merge.csv', sep=",",quote = ""))
# ^ for csv: as.data.frame(read.csv("file.csv",header = TRUE, sep = ",", quote = "",))
head(metadata)

SampNames<-metadata[,1:2]
# will change sample names from the Novogene sample names (ITS2_1,ITS2_2,etc) to SampleID names

#### Update Sample Names ####

t.its2_asvs<-as.data.frame(t(its2.ASV_counts))
t.its2_asvs[1:5,1:5]
t.its2_asvs$NewSampleName<-rownames(t.its2_asvs)

t.its2_asvs_update<-merge(SampNames,t.its2_asvs,by.x="NewSampleName",by.y="NewSampleName")
t.its2_asvs_update[1:5,1:5]

rownames(t.its2_asvs_update)<-t.its2_asvs_update$SampleID
t.its2_asvs_update[1:5,1:5]

# recreate ITS2 ASV count table with the actual sample names, not the ones Novogene gave...
its2_counts_newnames<-as.data.frame(t(t.its2_asvs_update[,-c(1:2)]))
its2_counts_newnames[1:5,1:5] # samples (with original IDs) as columns, ASVs as rows

# save total pre-contamination counts
TotalSampleCounts_withContams<-data.frame(SampleID=colnames(its2_counts_newnames[,!(colnames(its2_counts_newnames) %in% "ASV_ID")]),
                                          TotalCounts=colSums(its2_counts_newnames[,!(colnames(its2_counts_newnames) %in% "ASV_ID")]))
TotalSampleCounts_withContams

write.table(TotalSampleCounts_withContams, "data/ITS2_DADA2/ABS_ITS2_ASVs_TotalCounts_withContams.tsv",
            sep="\t", quote=F, col.names=T,row.names=F)


## Remove unwanted samples aka old Salton Seawater samples sequenced by Zymo
# remove_samples<-c("SS.OV.10m.seawater.0621", "SS.OV.2m.seawater.0621", "SS.OV.5m.seawater.0621")
# its2.ASV_counts<-its2.ASV_counts[,!(colnames(its2.ASV_counts) %in% remove_samples)]
# colnames(its2.ASV_counts)
# dim(its2.ASV_counts)

# #### Identify & Remove Contaminants ####
# ControlDF<-metadata[metadata$Flare=="Control",] # pull out samples that are controls
# 
# vector_for_decontam<-metadata$Sample_or_Control # use for decontam package
# # ^ tells us which are controls aka TRUE vs which are not aka FALSE
# 
# its2.ASV_counts[,-length(its2.ASV_counts)] <- as.data.frame(sapply(its2.ASV_counts[,-length(its2.ASV_counts)], as.numeric)) #convert data frame to numeric
# its2.ASV_c2<-t(its2.ASV_counts[,-length(its2.ASV_counts)]) # transpose so that rows are Samples and columns are ASVs
# contam_df <- isContaminant(its2.ASV_c2, neg=vector_for_decontam) # use isContaminant to identify ASVs that are contaminants
# ## ^ ASVs that are more prominant in negative controls are identified as contaminants because there is less actual DNA to compete with, allowing these contaminants to be amplifeid in neg. controls
# 
# table(contam_df$contaminant) # identify contaminants aka TRUE: 2156
# 
# contam_asvs <- (contam_df[contam_df$contaminant == TRUE, ]) # pull out ASV IDs for contaminating ASVs
# dim(contam_asvs)
# 
# its2.ASV_tax[row.names(its2.ASV_tax) %in% row.names(contam_asvs),] # see which taxa are contaminants
# 
# Contam.Taxa<-its2.ASV_tax[row.names(its2.ASV_tax) %in% row.names(contam_asvs),] # see which taxa are contaminants
# Contam.Taxa[1:4,]
# 
# write.table(Contam.Taxa, "data/ITS2_DADA2/ABS_ITS2_ContaminationASVs_List.tsv",
#             sep="\t", quote=F, col.names=T,row.names=F)
# 
# ## Create new files that EXCLUDE contaminants!!!
# 
# # making new fasta file
# #contam_indices <- which(asv_fasta %in% paste0(">", contam_asvs))
# #dont_want <- sort(c(contam_indices, contam_indices + 1))
# #asv_fasta_no_contam <- asv_fasta[- dont_want]
# 
# # making new count table
# its2.ASV_counts_no.contam <- its2.ASV_counts[!row.names(its2.ASV_counts) %in% row.names(contam_asvs), ] # drop ASVs found in contam_asvs
# head(its2.ASV_counts_no.contam)
# 
# # making new taxonomy table
# its2.ASV_tax.no.contam <- its2.ASV_tax[!row.names(its2.ASV_tax) %in% row.names(contam_asvs), ] # drop ASVs found in contam_asvs
# head(its2.ASV_tax.no.contam)
# 
# # Remove ASVs found in Controls from samples (in addition to contaminants previously ID'd by decontam)
# 
# Control_counts<-its2.ASV_counts_no.contam[,colnames(its2.ASV_counts_no.contam) %in% ControlDF$SampleID] # see which taxa are contaminants
# Control_counts
# Control_counts<-Control_counts[which(rowSums(Control_counts) > 0),] # drop ASVs that don't appear in Controls
# dim(Control_counts)
# head(Control_counts)
# 
# its2_counts_newnames<-its2.ASV_counts_no.contam[!its2.ASV_counts_no.contam$ASV_ID %in% row.names(Control_counts),!colnames(its2.ASV_counts_no.contam) %in% colnames(Control_counts)]
# its2.ASV_tax<-its2.ASV_tax.no.contam[!its2.ASV_tax.no.contam$ASV_ID %in% row.names(Control_counts),]
# 
# # sanity check
# colnames(its2_counts_newnames) # check for control sample IDs
# 
# ## and now writing them out to files
# #write(asv_fasta_no_contam, "ASVs-no-contam.fa")
# write.table(its2_counts_newnames, "data/ITS2_DADA2/ABS_ITS2_ASVs_Counts_NoContam.tsv",
#             sep="\t", quote=F, col.names=NA)
# saveRDS(its2_counts_newnames, file = "data/ITS2_DADA2/ABS_ITS2_ASVs_Counts_NoContam_Robject.rds", ascii = FALSE, version = NULL,
#         compress = TRUE, refhook = NULL)
# 
# write.table(its2.ASV_tax, "data/ITS2_DADA2/ABS_ITS2_ASVs_Taxa_NoContam.tsv",
#             sep="\t", quote=F, col.names=NA)
# saveRDS(its2.ASV_tax, file = "data/ITS2_DADA2/ABS_ITS2_ASVs_Taxa_NoContam_Robject.rds", ascii = FALSE, version = NULL,
#         compress = TRUE, refhook = NULL)

#### Update Metadata ####
# create color variable(s) to identify variables by colors
## color for Flare type
unique(metadata$Flare)
metadata$Flare<-factor(metadata$Flare, levels=c("HHP","Remission","Flare","Donor"))

colorset1 = reshape2::melt(c(HHP="#007561",Remission="#ff6361",Flare="#ffa600",Donor="grey"))

colorset1$Flare<-rownames(colorset1)
colnames(colorset1)[which(names(colorset1) == "value")] <- "Flare_Color"
colorset1

metadata<-merge(metadata, colorset1, by="Flare")
head(metadata)
metadata$Flare_Color <- as.character(metadata$Flare_Color)
rownames(metadata)<-metadata$SampleID
head(metadata)
dim(metadata)

#### Transform Data ####
# first we merge the ASV count object and the ASV taxonomy object together by column called "ASV_ID"
## then we need to melt the separate ASV & taxonomy so that we can rbind multiple data sets

# Drop singletons & zero count ASVs
dim(its2_counts_newnames) # 5424 x 63
its2_counts_newnames<-its2_counts_newnames[which(rowSums(its2_counts_newnames[,-length(its2_counts_newnames)]) > 0),]
dim(its2_counts_newnames) # 5375 x 63

# create ASV_ID column in count data and tax data for merging
its2_counts_newnames[1:5,1:5]
its2_counts_newnames$ASV_ID<-rownames(its2_counts_newnames)

head(its2.ASV_tax)
its2.ASV_tax$ASV_ID<-rownames(its2.ASV_tax)

# merge CLEAN aka singletons removed count & taxa tables
its2.ASV_all.tax<-merge(its2_counts_newnames,its2.ASV_tax, by="ASV_ID")
head(its2.ASV_all.tax)
dim(its2.ASV_all.tax)
#its2.ASV_all.tax<-its2.ASV_all.tax[, !duplicated(colnames(its2.ASV_all.tax))] # remove col duplicates
dim(its2.ASV_all.tax)
#its2.ASV_dat<-left_join(its2.ASV_tax.no.contam,its2.ASV_dat2, by=c("ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species"),all=T) # all=T to keep all rows from the dataframes, not drop rows that are missing from one DF or the other

its2.dat<-reshape2::melt(its2.ASV_all.tax,.by=vars("ASV_ID","Kingdom","Phylum","Class","Order","Family","Genus","Species"))
head(its2.dat)
colnames(its2.dat)[which(names(its2.dat) == "variable")] <- "SampleID"
colnames(its2.dat)[which(names(its2.dat) == "value")] <- "Count"

head(its2.dat)

# Create updated ASV table using correct SampleIDs (& with no singletons)
## ASV_table has SampleIDs as rows, ASVs as columns! this is for vegan package purposes
its2.ASV_table<-base::as.data.frame(dcast(its2.dat, SampleID~ASV_ID, value.var="Count", fun.aggregate=sum)) ###
its2.ASV_table[1:5,1:5]
its2.ASV_table[duplicated(rownames(its2.ASV_table))] # any samples duplicated here?
rownames(its2.ASV_table)<-its2.ASV_table$SampleID
#its2.ASV_table<-subset(its2.ASV_table, select=-c(SampleID))
rownames(its2.ASV_table)

# double check dimensions of metadata and ASV table
dim(metadata)
dim(its2.ASV_table)

# NOTE 12/3/24: Sample ITS2_8 was removed due to unsuccessful sequencing, so that is why there is one extra sample in metadata
# double check that the rownames exist + match
rownames(metadata)
rownames(its2.ASV_table)

# Find rows in metadata that are not in combined fungal asv tables
setdiff(rownames(metadata), rownames(its2.ASV_table)) # check rows in metadata not in its2.ASV_table
# ^ includes the sample that failed sequencing aka A020_10.12.20
setdiff(rownames(its2.ASV_table), rownames(metadata)) # check rows in its2.ASV_table not in metadata

# reorder metadata based off of ASV table - this will drop rows in metadata that are not in ASV table!
# here we are reordering our metadata by rows, using the rownames from our ASV table as a guide
# this indexing method will only work if the two dfs have the same row names!
metadata=metadata[rownames(its2.ASV_table),] ## will drop rows that are not shared by both dataframes!

dim(metadata) # confirm that we dropped the low quality sample from metadata
dim(its2.ASV_table)

# merge ASV table and metadata together just in case
its2.ASV_meta<-merge(its2.ASV_table,metadata, by.x="SampleID", by.y="SampleID") # drop missing sample
dim(its2.ASV_meta)
dim(metadata) # the difference between dimensions is that the original metadata contains info for the controls
dim(its2.ASV_table)
rownames(its2.ASV_meta)<-its2.ASV_meta$SampleID

# # the following section is performed if there are samples that we do not want to include in downstream analyses but were included in the sampling run!
# # recreate ASV taxonomy table so that ASVs that were removed are not included in new taxa table
# its2.tax.clean<-as.data.frame(unique(subset(its2.dat, select=-c(SampleID,Count))))
# head(its2.tax.clean)
# 
# # recreate ASV table (excluding samples we don't have metadata or counts for)
# colnames(its2.ASV_meta); rownames(its2.ASV_meta)
# head(its2.ASV_meta)
# its2.ASV_table<-subset(its2.ASV_meta, select=-c(Flare,TubeID,OriginalSampleID,SorC1,Sample_or_Control,Flare_Color))
# rownames(its2.ASV_table)<-its2.ASV_table$SampleID
# head(its2.ASV_table)
# 
# rowSums(its2.ASV_table[,-1])

# # save total post-contamination counts
# TotalSampleCounts_NoContams<-data.frame(SampleID=rownames(its2.ASV_table),
#                                           TotalCounts=rowSums(its2.ASV_table[,-1]))
# write.table(TotalSampleCounts_NoContams, "data/ITS2_DADA2/ABS_ITS2_ASVs_TotalCounts_NOContams.tsv",
#             sep="\t", quote=F, col.names=T,row.names=F)

# triple check dimensions of metadata and ASV table
dim(metadata)
dim(its2.ASV_table)

# Find rows in metadata that are not in combined fungal asv tables
setdiff(rownames(metadata), rownames(its2.ASV_table)) # check rows from metadata not in its2.ASV_table
# the difference between dimensions is that the original metadata contains info for the controls
setdiff(rownames(its2.ASV_table), rownames(metadata)) # check its2.ASV_table from metadata not in rows
# ^ should be no differences now!

# Create taxa table that doesn't have ASV_ID column
# head(its2.ASV_tax)
# its2.tax<-its2.ASV_tax[,-length(its2.ASV_tax)] # removes last column aka ASV_ID table from df
# head(its2.tax)

# create super metadata + taxa + counts data frame
its2.dat.meta<-merge(its2.dat,metadata,by.x="SampleID",by.y="SampleID")

#### Prep Dataframe for Relative Abundance ####
head(t(its2.ASV_table[,-1]))

its2.clean.counts<-as.data.frame(t(its2.ASV_table[,-1])) # all singletons, zeros, controls, contaminants dropped
head(its2.clean.counts)
its2.clean.counts$ASV_ID<-rownames(its2.clean.counts)

#its2.tax<-subset(its2.ASV_tax.no.contam, select=c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
its2.all<-merge(its2.clean.counts, its2.ASV_tax, by="ASV_ID")
head(its2.all)
its2_melt<-reshape2::melt(its2.all, id.vars = c("ASV_ID", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"))
head(its2_melt)
names(its2_melt)[which(names(its2_melt) == "variable")] <- "SampleID"
names(its2_melt)[which(names(its2_melt) == "value")] <- "Counts"
head(its2_melt)

head(metadata)

all_its2<-merge(its2_melt, metadata, by = "SampleID")
head(all_its2) # contains metadata, ASV counts, and taxonomic IDs for ASVs

#### Save Files ####

# save updated taxa file
saveRDS(its2.tax.clean, file = "data/ITS2_DADA2/ABS_ITS2_ASVs_Taxonomy_dada2_Clean_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#### Pull out ASV Sequences of Interest for MSA & Phylogenetic Tree ####
## pulling sequences from DADA2 results only for ASVs that we have in our cleaned up its2.ASV_table
## use this new set of ASVs of interest to construct phylogenetic tree, import the tree with ape package, and calculate Unifrac distance with rbiom
# giving our seq headers more manageable names (ASV_1, ASV_2...)
# load(file = "data/ITS2_DADA2/mydada_16S.V3V4.Rdata") # load your DADA2 output for this project
# 
# asv_seqs <- colnames(seqtab.nochim) # pull out ASV sequences
# asv_headers <- vector(dim(seqtab.nochim)[2], mode="character") # pull out ASV headers
# 
# # add ASV in front of each # of each sequence
# for (i in 1:dim(seqtab.nochim)[2]) {
#   asv_headers[i] <- paste("ASV", i, sep="_")
# }
# 
# # create table of ASVs and their respective sequences so that we can filter out ASVs from samples we aren't looking at
# asv_all.tax <- as.data.frame(cbind(asv_headers, asv_seqs))
# length(unique(asv_all.tax$asv_headers))
# 
# # only pull out ASVs and their sequences for ASVs in our taxa table for only samples of interest
# new.asv_all.tax<-asv_all.tax[asv_all.tax$asv_headers %in% colnames(its2.ASV_table[,-1]),]
# 
# # insert ">" before ASV IDs before we make the new FASTA file
# new.asv_all.tax$asv_headers<-gsub("ASV_",">ASV_",new.asv_all.tax$asv_headers)
# new.asv_all.tax$asv_headers[1:5]
# 
# # create and save the new fasta file - will use this to generate Multiple Sequence Alignment & build phylogenetic tree
# new.asv_fasta<-c(rbind(new.asv_all.tax$asv_headers,new.asv_all.tax$asv_seqs))
# write(new.asv_fasta, "SSD_16SV3V4_ASVs_Updated.fa")

#### Save Global Env for Import into Other Scripts ####

save.image("data/ITS2_DADA2/ABS_ITS2_DataReady.Rdata") # save global env to Rdata file

## ^ this will be loaded into other scripts for downstream analyses

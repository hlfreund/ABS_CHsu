#### Load Packages & Set WD ####
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
  library(data.table)
  library(forcats)
})

getwd()
#setwd("/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/data/Metagenomes/Revisions_3.3.2025/adhE_Gene_Search")
list.files() # sanity check

#### Import ABS & HC Combined Metadata ####

# this metadata is a file I created using the SRA submission file + added propensity-matched HCs to spreadsheet
abs.all.basic.meta<-as.data.frame(read_xlsx("data/Metagenomes/SRA_Metagenomes_metadata_basic.xlsx", sheet="AllSamples_List"))
head(abs.all.basic.meta)

abs.all.basic.meta$SampleID<-gsub("_L00[0-9]_R1_001.fastq.gz","",abs.all.basic.meta$filename)

# create group color palette to merge with abs.all.basic.meta
grp.clrs = as.data.frame(t(data.frame("HC"="lightgray","Flare"="#ffa600","Remission"="#ff6361","HHP"="#007561")))
grp.clrs

grp.clrs$Condition<-rownames(grp.clrs)
colnames(grp.clrs)[which(names(grp.clrs) == "V1")] <- "Group_Color"
grp.clrs

abs.all.basic.meta<-merge(abs.all.basic.meta,grp.clrs,by="Condition")
unique(abs.all.basic.meta$Group_Color)

# #### Import the adhE Read Depth Data ####
# # first the adhE depth in ABS mgms
# adhE.depth<-fread(file = 'data/Metagenomes/Revisions_3.3.2025/adhE_Gene_Search/ABS_adhE_bowtie2_Depth_Results.tsv', sep='\t',header = TRUE)
# head(adhE.depth)
# adhE.depth.totals<-dcast(as.data.table(adhE.depth), SampleID ~ Chrom, value.var = "ReadCount", sum)
# rownames(adhE.depth.totals)<-adhE.depth.totals$SampleID
# head(adhE.depth.totals)
# 
# adhE.depth.tots.m<-melt(adhE.depth.totals,by="SampleID")
# names(adhE.depth.tots.m)
# colnames(adhE.depth.tots.m)[which(names(adhE.depth.tots.m) == "variable")] <- "gene_id"
# colnames(adhE.depth.tots.m)[which(names(adhE.depth.tots.m) == "value")] <- "total_read_count"
# head(adhE.depth.tots.m)
# 
# # then the adhE depth in propensity-matched HCs mgms
# adhE.HC.depth<-fread(file = 'data/Metagenomes/Revisions_3.3.2025/adhE_Gene_Search/ABS_HCs_adhE_bowtie2_Depth_Results.tsv', sep='\t',header = TRUE)
# head(adhE.HC.depth)
# adhE.HC.depth.totals<-dcast(as.data.table(adhE.HC.depth), SampleID ~ Chrom, value.var = "ReadCount", sum)
# rownames(adhE.HC.depth.totals)<-adhE.HC.depth.totals$SampleID
# head(adhE.HC.depth.totals)
# 
# adhE.HC.depth.tots.m<-melt(adhE.HC.depth.totals,by="SampleID")
# names(adhE.HC.depth.tots.m)
# colnames(adhE.HC.depth.tots.m)[which(names(adhE.HC.depth.tots.m) == "variable")] <- "gene_id"
# colnames(adhE.HC.depth.tots.m)[which(names(adhE.HC.depth.tots.m) == "value")] <- "total_read_count"
# head(adhE.HC.depth.tots.m)
# 
# 
# adhE.depth.all<-rbind(adhE.depth.tots.m,adhE.HC.depth.tots.m)
# 
#### Import the adhE Read Coverage Data ####
## note: changed first column header manually in tsv before doing this step!
## ^ otherwise the column headers are not read in properly by read.table()

# first the adhE coverage in ABS mgms
adhE.cov<-read.table("data/Metagenomes/Revisions_3.3.2025/adhE_Gene_Search/ABS_adhE_bowtie2_Coverage_Results.tsv",header=TRUE,sep="\t")
names(adhE.cov) # look at column headers to ensure this imported correctly

# remove some sequencing info from sample ID names
adhE.cov$SampleID<-gsub("_L00[0-9]","",adhE.cov$SampleID)

# drop Seq18_adhE (from pathogenic Ecoli; it's nearly identical to non-pathogenic Ecoli adhE)
adhE.cov<-adhE.cov[adhE.cov$gene_id!="Seq18_adhE",]

# import metadata about adhE genes and the genomes they came from
adhE.meta<-read.table("data/Metagenomes/Revisions_3.3.2025/adhE_Gene_Search/adhE_gene_genome_metadata.txt",sep="\t")
names(adhE.meta) # look at column headers
colnames(adhE.meta)[which(names(adhE.meta) == "V1")] <- "gene_id"
colnames(adhE.meta)[which(names(adhE.meta) == "V2")] <- "genus"
colnames(adhE.meta)[which(names(adhE.meta) == "V3")] <- "species"

head(adhE.cov)

# then import the adhE coverage in propensity-matched HC mgms
adhE.HC.cov<-read.table("data/Metagenomes/Revisions_3.3.2025/adhE_Gene_Search/ABS_HCs_adhE_bowtie2_Coverage_Results.tsv",header=TRUE,sep="\t")
names(adhE.HC.cov) # look at column headers to ensure this imported correctly
head(adhE.HC.cov)

# drop Seq18_adhE (from pathogenic Ecoli)
adhE.HC.cov<-adhE.HC.cov[adhE.HC.cov$gene_id!="Seq18_adhE",]

# use rbind to combine the ABS and HC adhE coverage data together
adhE.covs.all<-rbind(adhE.cov,adhE.HC.cov)

# use merge to combine this massive adhE coverage df with the adhE metadata
adhE.all<-merge(adhE.meta,adhE.covs.all,by="gene_id")
head(adhE.all)

# order the combined data by SampleID, makes for easier review of the data in table form
adhE.all<-adhE.all[order(adhE.all$SampleID),]

# subset out samples with coverage over 50%, meandepth > 5
adhE.subset<-adhE.all[ which(adhE.all$coverage > 50 & adhE.all$meandepth > 5) , ]
# note: seem to lose HC samples when we do this filtering becuase they had such low depth of coverage for adhE genes...

# save this table!
# write.table(adhE.subset,file="data/Metagenomes/Revisions_3.3.2025/adhE_Gene_Search/ABS_and_Matched_HCs_adhE_bowtie2_coverage_with_metadata_updated.tsv",sep="\t",col.names=TRUE,row.names=FALSE)
## NOTE: saving subsetted table but also want to visualize with all adhE results; be mindful of the dfs you're referring to below

# more info on breath vs depth of coverage here:
## https://www.metagenomics.wiki/pdf/qc/coverage-read-depth
## https://www.metagenomics.wiki/tools/samtools/breadth-of-coverage

#### Prep adhE for Visualization ####

# create long (dcast) version of adhE.all
adhE.table<-reshape2::dcast(adhE.all, SampleID ~ gene_id, value.var = "meandepth",sum)
rownames(adhE.table)<-adhE.table$SampleID

# melt down long version of adhE data for visualizing
adhE.melt<-reshape2::melt(adhE.table,by="SampleID")
colnames(adhE.melt)[which(names(adhE.melt) == "variable")] <- "gene_id"
colnames(adhE.melt)[which(names(adhE.melt) == "value")] <- "meandepth"

# merge melted data with adhE data
adhE.melt2<-merge(adhE.meta,adhE.melt,by="gene_id")
head(adhE.melt2)

# merge combined ABS metadata with adhE data
head(abs.all.basic.meta)

adhE.melt3<-merge(abs.all.basic.meta[,c(1,4:6)],adhE.melt2,by="SampleID")
head(adhE.melt3)
adhE.melt3$Condition<-factor(adhE.melt3$Condition,levels=c("HC","Flare","Remission","HHP"))

#### Visualize adhE by Taxa & Sample ####

# visualize with stacked barplot
ggplot(adhE.melt3, aes(x=SampleID, y=meandepth, fill=genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
        legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
        plot.title = element_text(size=30),plot.subtitle = element_text(size=26))+
  guides(fill=guide_legend(ncol=1)) + 
  labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus")+coord_flip()

#ggsave(adhe1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_barplot.png", width = 6000,height = 7000,units = "px",dpi = 200,create.dir = TRUE)

# barplot wrapped by condition
adhe1a<-ggplot(adhE.melt3, aes(x=SampleID, y=meandepth, fill=genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
        axis.text.x = element_text(), legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
        plot.title = element_text(size=30),plot.subtitle = element_text(size=26),strip.text.y = element_text(size=26, face = "bold"),
        strip.background=element_rect(fill="ghostwhite"))+
  guides(fill=guide_legend(ncol=1)) + facet_grid(Condition~., scales = "free_y",space="free_y") +
  labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus")+coord_flip()

ggsave(adhe1a,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_Condition_barplot.png", width = 6500,height = 8500,units = "px",dpi = 200,create.dir = TRUE)

# barplot with genus & species, faceted by condition
adhe1b<-ggplot(adhE.melt3, aes(x=SampleID, y=meandepth, fill=interaction(genus,species,sep=" ")))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
        legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
        plot.title = element_text(size=30),plot.subtitle = element_text(size=26),strip.text.y = element_text(size=26, face = "bold"),
        strip.background=element_rect(fill="ghostwhite"))+
  guides(fill=guide_legend(ncol=1)) + facet_grid(Condition~., scales = "free_y",space="free_y") +
  labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus")+coord_flip()

ggsave(adhe1b,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_species_Condition_barplot.png", width = 6500,height = 8500,units = "px",dpi = 200,create.dir = TRUE)

# visualize with heatmap
ggplot(adhE.melt3, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=100,breaks=c(0,50,100,150,200),
                       labels=c("0","50","100","150","200"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

#ggsave(adhe1c,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_heatmap.png", width = 6000,height = 7500,units = "px",dpi = 200,create.dir = TRUE)

# heatmap facet wrapped by condition
adhe1d<-ggplot(adhE.melt3, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=100,breaks=c(0,50,100,150,200),
                       labels=c("0","50","100","150","200"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01,face="italic"),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") + facet_grid(Condition~., scales = "free_y",space="free_y") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(adhe1d,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_Condition_heatmap.png", width = 6500,height = 8500,units = "px",dpi = 200,create.dir = TRUE)

# heatmap with genus & species, faceted by condition
adhe1e<-ggplot(adhE.melt3, aes(SampleID, interaction(genus,species,sep=" "), fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=100,breaks=c(0,50,100,150,200),
                       labels=c("0","50","100","150","200"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01,face="italic"),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") + facet_grid(Condition~., scales = "free_y",space="free_y") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(adhe1e,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_species_Condition_heatmap.png", width = 7000,height = 9500,units = "px",dpi = 200,create.dir = TRUE)

# #### Visualize Higher Coverage adhE Results ####
# 
# # subset out samples with coverage over 50%, meandepth > 5
# adhE.subset<-adhE.all[ which( adhE.all$coverage > 50 & adhE.all$meandepth > 5) , ]
# # note: seem to lose HC samples when we do this filtering becuase they had such low depth of coverage for adhE genes...
# min(adhE.subset$meandepth) # sanity check that our filtering work...
# min(adhE.subset$coverage) # sanity check that our filtering work...
# 
# # create long (dcast) version of adhE.all
# adhE.subset.table<-reshape2::dcast(adhE.subset, SampleID ~ gene_id, value.var = "meandepth",sum)
# rownames(adhE.subset.table)<-adhE.subset.table$SampleID
# 
# # melt down long version of adhE data for visualizing
# adhE.subset.melt<-reshape2::melt(adhE.subset.table,by="SampleID")
# head(adhE.subset.melt)
# # rename column names for merging
# colnames(adhE.subset.melt)[which(names(adhE.subset.melt) == "variable")] <- "gene_id"
# colnames(adhE.subset.melt)[which(names(adhE.subset.melt) == "value")] <- "meandepth"
# head(adhE.subset.melt)
# 
# # drop the zeros that came from creating subset table
# adhE.subset.melt<-adhE.subset.melt[adhE.subset.melt$meandepth!=0,]
# head(adhE.subset.melt)
# 
# # merge subsetted adhE cov data with adhE metadata
# adhE.subset.melt2<-merge(adhE.meta,adhE.subset.melt,by="gene_id")
# head(adhE.subset.melt2)
# 
# # merge combined ABS metadata with adhE data
# adhE.subset.melt3<-merge(abs.all.basic.meta[,c(1,4:5)],adhE.subset.melt2,by="SampleID")
# head(adhE.subset.melt3)
# adhE.subset.melt3$Condition<-factor(adhE.subset.melt3$Condition,levels=c("HC","Flare","Remission","HHP"))
# 
# # stacked barplot
# ggplot(adhE.subset.melt3, aes(x=SampleID, y=meandepth, fill=genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
#         legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
#         plot.title = element_text(size=30),plot.subtitle = element_text(size=26))+
#   guides(fill=guide_legend(ncol=1)) + 
#   labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus",subtitle="Contains adhE genes with > 50% breadth of coverage and > 5 mean depth of coverage")+coord_flip()
# 
# #ggsave(adhe2,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_High_coverage_bowtie2_subset_barplot.png", width = 3000,height = 4000,units = "px",dpi = 200,create.dir = TRUE)
# 
# # barplot wrapped by condition
# adhe2a<-ggplot(adhE.subset.melt3, aes(x=SampleID, y=meandepth, fill=genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
#         legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
#         plot.title = element_text(size=30),plot.subtitle = element_text(size=26),strip.text.y = element_text(size=26, face = "bold"),
#         strip.background=element_rect(fill="ghostwhite"))+
#   guides(fill=guide_legend(ncol=1)) + facet_grid(Condition~., scales = "free_y",space="free_y") +
#   labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus")+coord_flip()
# 
# ggsave(adhe2a,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_High_coverage_bowtie2_subset_genus_Condition_barplot.png", width = 6000,height = 7000,units = "px",dpi = 200,create.dir = TRUE)
# 
# # barplot with genus & species, faceted by condition
# adhe2b<-ggplot(adhE.subset.melt3, aes(x=SampleID, y=meandepth, fill=interaction(genus,species,sep=" ")))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
#         legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
#         plot.title = element_text(size=30),plot.subtitle = element_text(size=26),strip.text.y = element_text(size=26, face = "bold"),
#         strip.background=element_rect(fill="ghostwhite"))+
#   guides(fill=guide_legend(ncol=1)) + facet_grid(Condition~., scales = "free_y",space="free_y") +
#   labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus species")+coord_flip()
# 
# ggsave(adhe2b,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_High_coverage_bowtie2_subset_genus_species_Condition_barplot.png", width = 6000,height = 8000,units = "px",dpi = 200,create.dir = TRUE)
# 
# # heatmap
# ggplot(adhE.subset.melt3, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
#   scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=100,breaks=c(0,50,100,150,200),
#                        labels=c("0","50","100","150","200"))+
#   theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
#                         axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
#                         legend.title = element_text(hjust=0.5,size=30),
#                         legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
#                         strip.text.y = element_text(size=26, face = "bold"), plot.subtitle = element_text(size=26), strip.background=element_rect(fill="ghostwhite")) +
#   labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains adhE genes with > 50% breadth of coverage and > 5 mean depth of coverage") +
#   guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()
# 
# #ggsave(adhe2b,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_High_coverage_bowtie2_subset_heatmap.png", width = 4000,height = 5000,units = "px",dpi = 200,create.dir = TRUE)
# 
# # heatmap facet wrapped by condition
# adhe2c<-ggplot(adhE.subset.melt3, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
#   scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=100,breaks=c(0,50,100,150,200),
#                        labels=c("0","50","100","150","200"))+
#   theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
#                         axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01,face="italic"),
#                         legend.title = element_text(hjust=0.5,size=30),
#                         legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
#                         strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
#   labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") + facet_grid(Condition~., scales = "free_y",space="free_y") +
#   guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()
# 
# ggsave(adhe2c,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_High_coverage_bowtie2_subset_genus_Condition_heatmap.png", width = 6000,height = 7500,units = "px",dpi = 200,create.dir = TRUE)
# 
# # heatmap with genus & species, faceted by condition
# adhe2d<-ggplot(adhE.subset.melt3, aes(SampleID, interaction(genus,species,sep=" "), fill= meandepth)) +geom_tile(color="white")+
#   scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=100,breaks=c(0,50,100,150,200),
#                        labels=c("0","50","100","150","200"))+
#   theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
#                         axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01,face="italic"),
#                         legend.title = element_text(hjust=0.5,size=30),
#                         legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
#                         strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
#   labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") + 
#   facet_grid(Condition~., scales = "free_y",space="free_y") +
#   guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()
# 
# ggsave(adhe2d,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_High_coverage_bowtie2_subset_genus_species_Condition_heatmap.png", width = 6000,height = 8000,units = "px",dpi = 200,create.dir = TRUE)
# 
#### Visualize Ecoli/Kpneumoniae adhE by Taxa - All Samples ####

# subset out Ecoli and Kpneumoniae adhE mean depth of coverage only
head(adhE.melt3)

ecoli.kpneum.adhE<-subset(adhE.melt3,species==c("coli","pneumoniae"))
ecoli.kpneum.adhE$gen_spec<-as.character(interaction(ecoli.kpneum.adhE$genus,ecoli.kpneum.adhE$species,sep=" "))
head(ecoli.kpneum.adhE)

# visualize with stacked barplot

# barplot with genus & species, faceted by condition
ec.kp.brp1<-ggplot(ecoli.kpneum.adhE, aes(x=SampleID, y=meandepth, fill=gen_spec))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
        legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
        plot.title = element_text(size=30),plot.subtitle = element_text(size=26),strip.text.y = element_text(size=26, face = "bold"),
        strip.background=element_rect(fill="ghostwhite"))+
  guides(fill=guide_legend(ncol=1)) + facet_grid(Condition~., scales = "free_y",space="free_y") +
  labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus species")+coord_flip()

ggsave(ec.kp.brp1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_Ecoli_Kpneumoniae_Condition_barplot.png", width = 5500,height = 7000,units = "px",dpi = 200,create.dir = TRUE)

# boxplot time

# create comparisons for the boxplots
mycomp <- list(c("HC", "Flare"), c("Flare", "Remission"), c("HC", "Remission"), c("HC", "HHP"),
               c("HHP", "Remission"),c("HHP", "Flare"))
mycomp1 <- list(c("Flare", "Remission"))

# don't run next line, just a reference for colors
# "HC"="lightgray","Flare"="#ffa600","Remission"="#ff6361","HHP"="#007561"

ec.kp.bxp1<-ggboxplot(
  ecoli.kpneum.adhE, x="Condition", y="meandepth", 
  add = c("jitter"), 
  facet.by = "gen_spec", nrow = 1,
  fill = "Condition", palette = c("lightgray","#ffa600", "#ff6361", "#007561"),size=1)+
  stat_compare_means(aes(group = Condition),label.y = 65,size=6) +
  stat_compare_means(aes(group = Condition),comparisons = mycomp,method="wilcox")+
  scale_y_continuous(name = "Mean Depth of Coverage",limits=c(0,65))+
  rremove("x.ticks")+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(colour = "black", size=22),
        axis.text=element_text(colour="black", size=22),
        axis.text.x = element_text(),
        strip.text.x = element_text(size=24, face = "italic"),strip.background=element_rect(fill="ghostwhite"))

ggsave(ec.kp.bxp1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_Ecoli_Kpneumoniae_Condition_boxplot.png",width = 2400,height = 2000,units = "px",dpi = 200,create.dir = TRUE)

# visualize with heatmap
ggplot(ecoli.kpneum.adhE, aes(SampleID, gen_spec, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=30,breaks=c(0,10,20,30,40,50,60),
                       labels=c("0","10","20","30","40","50","60"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_grid(Condition~., scales = "free_y",space="free_y")+coord_flip()

#ggsave(adhe1c,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_Ecoli_Kpneumoniae_heatmap.png", width = 6000,height = 7500,units = "px",dpi = 200,create.dir = TRUE)

# heatmap with genus & species, faceted by condition
ec.kp.hm1<-ggplot(ecoli.kpneum.adhE, aes(SampleID, interaction(genus,species,sep=" "), fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",,midpoint=30,breaks=c(0,10,20,30,40,50,60),
                       labels=c("0","10","20","30","40","50","60"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(face="italic"),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") + facet_grid(Condition~., scales = "free_y",space="free_y") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(ec.kp.hm1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_Ecoli_Kpneumoniae_Condition_heatmap.png",width = 5500,height = 7000,units = "px",dpi = 200,create.dir = TRUE)


#### Import FMT Metadata ####

fmt.metadata<-as.data.frame(read_xlsx("/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/data/ABS_FMT_Metadata_CLH.xlsx", sheet="Sheet1"))
fmt.metadata$SampleID<-gsub("_L00[0-9]_R1_001.fastq.gz","",fmt.metadata$filename)
unique(fmt.metadata$SampleID)

sample.order<-unique(fmt.metadata$SampleID)

# curiosity: which FMT patient/HHP/Donor samples are found in the adhE data that had High coverage?
fmt.metadata$SampleID[which(unique(fmt.metadata$SampleID) %in% unique(adhE.subset$SampleID))]

# view the df containing all adhE cov data used for visualization above...
head(adhE.melt3)

# subset out only FMT/FMT-HHP/Donor adhE cov data
adhE.fmt<-adhE.melt3[which(adhE.melt3$SampleID %in% fmt.metadata$SampleID),]
head(adhE.fmt)
unique(adhE.fmt$SampleID) %in% unique(fmt.metadata$SampleID) # sanity check, should all be TRUE

# drop rows with meandepth=0
adhE.fmt<-adhE.fmt[adhE.fmt$meandepth>0,]
min(adhE.fmt$meandepth) # sanity check

# merge FMT/FMT-HHP/Donor adhE cov data with FMT metadata
adhE.fmt.meta<-merge(adhE.fmt,fmt.metadata,by=c("SampleID","Condition","Patient"))
head(adhE.fmt.meta)

# create combined genus species column, easier for plotting later
adhE.fmt.meta$gen_spec<-interaction(adhE.fmt.meta$genus,adhE.fmt.meta$species,sep=" ")

# turn adhE.fmt.meta$FMT_Label & Condition columns into factors for ordering
unique(adhE.fmt.meta$FMT_Label)
adhE.fmt.meta$FMT_Label<-factor(adhE.fmt.meta$FMT_Label,
                                levels=c("Pre-FMT","Abx-1",
                                         "Post-FMT1","Abx-2","Post-FMT2","Donor"))

unique(adhE.fmt.meta$Condition)
adhE.fmt.meta$Condition<-factor(adhE.fmt.meta$Condition,levels=c("HC","Flare","Remission","HHP"))

adhE.fmt.meta$Patient<-factor(adhE.fmt.meta$SampleID,levels=c("A020","A021","Donor"))

adhE.fmt.meta$SampleID<-factor(adhE.fmt.meta$SampleID,levels=c(sample.order))

# order by sampleID and Patient since they are in chronological order
adhE.fmt.meta<-adhE.fmt.meta[order(adhE.fmt.meta$SampleID,adhE.fmt.meta$Patient),]

#### Visualize FMT adhE data ####

# stacked barplot, first wrapped by FMT label
adhe.fmt1<-ggplot(adhE.fmt.meta, aes(x=forcats::fct_rev(SampleID), y=meandepth, fill=gen_spec))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
        legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
        plot.title = element_text(size=30),plot.subtitle = element_text(size=26),strip.text.y = element_text(size=20, face = "bold"),
        strip.background=element_rect(fill="ghostwhite"))+
  guides(fill=guide_legend(ncol=1)) + 
  labs(x="Sample ID", y="Mean Depth of Coverage", title="Contains only FMT Patient + FMT HHP + Donor Data",fill="Genus")+
  facet_grid(FMT_Label~., scales = "free_y",space="free_y")+coord_flip()

ggsave(adhe.fmt1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_genus_FMTlabel_barplot.png", width = 4000,height = 5500,units = "px",dpi = 200,create.dir = TRUE)

# barplot wrapped by condition
adhe.fmt1a<-ggplot(adhE.fmt.meta, aes(x=forcats::fct_rev(SampleID), y=meandepth, fill=gen_spec))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
        legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
        plot.title = element_text(size=30),plot.subtitle = element_text(size=26),
        strip.background=element_rect(fill="ghostwhite"),strip.text.y = element_text(size=26, face = "bold"))+
  guides(fill=guide_legend(ncol=1)) + facet_grid(Condition~., scales = "free_y",space="free_y") +
  labs(x="Sample ID", y="Mean Depth of Coverage", title="Contains only FMT Patient + FMT HHP + Donor Data",fill="Genus")+
  facet_grid(Condition~., scales = "free_y",space="free_y")+coord_flip()

ggsave(adhe.fmt1a,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_genus_Condition_barplot.png", width = 3000,height = 4000,units = "px",dpi = 200,create.dir = TRUE)

# heatmap
adhe.fmt2<-ggplot(adhE.fmt.meta, aes(SampleID, gen_spec, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
                       labels=c("0","50","100","150"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=30),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01,face="italic"),
                        legend.title = element_text(hjust=0.5,size=28),
                        legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=30, face = "bold"), plot.subtitle = element_text(size=30), strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_grid(Condition~., scales = "free_y",space="free_y")+coord_flip()

ggsave(adhe.fmt2,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_heatmap.png", width = 7700,height = 6000,units = "px",dpi = 200,create.dir = TRUE)

# # heatmap 2
# adhe.fmt2b<-ggplot(adhE.fmt.meta, aes(SampleID, gen_spec, fill= meandepth)) +geom_tile(color="white")+
#   scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
#                        labels=c("0","50","100","150"))+
#   theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=30),
#                         axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
#                         legend.title = element_text(hjust=0.5,size=28),
#                         legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
#                         strip.text.y = element_text(size=30, face = "bold"), plot.subtitle = element_text(size=30), strip.background=element_rect(fill="ghostwhite")) +
#   labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
#   guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
#   facet_grid(.~FMT_Label, scales = "free_y",space="free_y")+coord_flip()
# 
# ggsave(adhe.fmt2b,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_heatmap2.png", width = 5000,height = 6000,units = "px",dpi = 200,create.dir = TRUE)

# heatmap 3
# adhe.fmt2c<-ggplot(adhE.fmt.meta, aes(SampleID, gen_spec, fill= meandepth)) +geom_tile(color="white")+
#   scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
#                        labels=c("0","50","100","150"))+
#   theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=30),
#                         axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),axis.text.y = element_text(face="italic"),
#                         legend.title = element_text(hjust=0.5,size=28),
#                         legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
#                         strip.text.y = element_text(size=30, face = "bold"), plot.subtitle = element_text(size=30), strip.background=element_rect(fill="ghostwhite")) +
#   labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
#   guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
#   facet_wrap(Condition~FMT_Label, scales = "free")
# 
# ggsave(adhe.fmt2c,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_heatmap3.png", width = 8000,height = 6000,units = "px",dpi = 200,create.dir = TRUE)

# heatmap 4
adhe.fmt2d<-ggplot(adhE.fmt.meta, aes(SampleID, gen_spec, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
                       labels=c("0","50","100","150"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=30),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),axis.text.y = element_text(face="italic"),
                        legend.title = element_text(hjust=0.5,size=28),
                        legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=30, face = "bold"), plot.subtitle = element_text(size=30), strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_grid(Condition~FMT_Label, scales = "free",space="free")

ggsave(adhe.fmt2d,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_heatmap4.png", width = 9000,height = 9500,units = "px",dpi = 200,create.dir = TRUE)

# heatmap without HHP
# heatmap 5
adhe.fmt2e<-ggplot(adhE.fmt.meta[adhE.fmt.meta$Condition!="HHP",], aes(SampleID, gen_spec, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
                       labels=c("0","50","100","150"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=30),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),axis.text.y = element_text(face="italic"),
                        legend.title = element_text(hjust=0.5,size=28),
                        legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=30, face = "bold"), plot.subtitle = element_text(size=30), strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_grid(Condition~FMT_Label, scales = "free",space="free")

ggsave(adhe.fmt2e,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_NoHHP_heatmap5.png", width = 6000,height = 7000,units = "px",dpi = 200,create.dir = TRUE)

#### Visualize FMT Only & Ecoli/Kpneumoniae adhE ####

head(adhE.fmt.meta)
adhE.fmt.meta.eckp<-adhE.fmt.meta[adhE.fmt.meta$species %in% c("coli","pneumoniae"),]
head(adhE.fmt.meta.eckp)

# visualize with stacked barplot

# barplot with genus & species, faceted by condition
fmt.ec.kp.brp1<-ggplot(adhE.fmt.meta.eckp, aes(x=SampleID, y=meandepth, fill=gen_spec))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=28),axis.title.y = element_text(size=28),axis.text = element_text(size=26),
        legend.title = element_text(size=26,hjust=0),legend.text = element_text(size=26,face="italic"),
        plot.title = element_text(size=30),plot.subtitle = element_text(size=26),strip.text.y = element_text(size=26, face = "bold"),
        strip.background=element_rect(fill="ghostwhite"))+
  guides(fill=guide_legend(ncol=1)) + facet_grid(Condition~., scales = "free_y",space="free_y") +
  labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus species")+coord_flip()

ggsave(fmt.ec.kp.brp1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_Ecoli_Kpneumoniae_Condition_barplot.png", width = 5500,height = 7000,units = "px",dpi = 200,create.dir = TRUE)

# boxplot time

# create comparisons for the boxplots
mycomp <- list(c("HC", "Flare"), c("Flare", "Remission"), c("HC", "Remission"), c("HC", "HHP"),
               c("HHP", "Remission"),c("HHP", "Flare"))
mycomp1 <- list(c("Flare", "Remission"))

# don't run next line, just a reference for colors
# "HC"="lightgray","Flare"="#ffa600","Remission"="#ff6361","HHP"="#007561"

fmt.ec.kp.bxp1<-ggboxplot(
  adhE.fmt.meta.eckp, x="Condition", y="meandepth", 
  add = c("jitter"), 
  facet.by = "gen_spec", nrow = 1,
  fill = "Condition", palette = c("lightgray","#ffa600", "#ff6361", "#007561"),size=1)+
  stat_compare_means(aes(group = Condition),label.y = 65,size=6) +
  stat_compare_means(aes(group = Condition),comparisons = mycomp,method="wilcox")+
  scale_y_continuous(name = "Mean Depth of Coverage",limits=c(0,65))+
  rremove("x.ticks")+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(colour = "black", size=22),
        axis.text=element_text(colour="black", size=22),
        axis.text.x = element_text(colour="black", size=21),
        strip.text.x = element_text(size=24, face = "italic"),strip.background=element_rect(fill="ghostwhite"))

ggsave(fmt.ec.kp.bxp1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_Ecoli_Kpneumoniae_Condition_boxplot.png", width = 2400,height = 2000,units = "px",dpi = 200,create.dir = TRUE)

# visualize with heatmap
ggplot(adhE.fmt.meta.eckp, aes(SampleID, gen_spec, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3")+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_grid(Condition~., scales = "free_y",space="free_y")+coord_flip()

# heatmap with genus & species, faceted by condition
fmt.ec.kp.hm1<-ggplot(adhE.fmt.meta.eckp, aes(SampleID, interaction(genus,species,sep=" "), fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=40,breaks=c(0,20,40,60,80),
                       labels=c("0","20","40","60","80"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(face="italic"),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") + facet_grid(Condition~., scales = "free_y",space="free_y") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(fmt.ec.kp.hm1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_Ecoli_Kpneumoniae_Condition_heatmap.png",width = 5500,height = 7000,units = "px",dpi = 200,create.dir = TRUE)

fmt.ec.kp.hm2<-ggplot(adhE.fmt.meta.eckp, aes(SampleID, interaction(genus,species,sep=" "), fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=40,breaks=c(0,20,40,60,80),
                       labels=c("0","20","40","60","80"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=30),axis.text = element_text(size=26),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01,face="italic"),
                        legend.title = element_text(hjust=0.5,size=30),
                        legend.text = element_text(size=30),plot.title = element_text(size=30), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=26, face = "bold"),strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") + facet_grid(Patient~FMT_Label, scales = "free_y",space="free_y") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(fmt.ec.kp.hm2,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_Ecoli_Kpneumoniae_Condition_heatmap2.png",width = 5500,height = 7000,units = "px",dpi = 200,create.dir = TRUE)

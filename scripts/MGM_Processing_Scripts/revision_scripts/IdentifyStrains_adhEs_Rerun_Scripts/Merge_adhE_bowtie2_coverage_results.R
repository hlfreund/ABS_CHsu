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
setwd("/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/data/Metagenomes/Revisions_3.3.2025/adhE_Gene_Search")
list.files() # sanity check

# #### Import the adhE Read Depth Data ####
# # first the adhE depth in ABS mgms
# adhE.depth<-fread(file = 'ABS_adhE_bowtie2_Depth_Results.tsv', sep='\t',header = TRUE)
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
# adhE.HC.depth<-fread(file = 'ABS_HCs_adhE_bowtie2_Depth_Results.tsv', sep='\t',header = TRUE)
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
adhE.cov<-read.table("ABS_adhE_bowtie2_Coverage_Results.tsv",header=TRUE,sep="\t")
names(adhE.cov) # look at column headers to ensure this imported correctly

# remove some sequencing info from sample ID names
adhE.cov$SampleID<-gsub("_L00[0-9]","",adhE.cov$SampleID)

# drop Seq18_adhE (from pathogenic Ecoli)
adhE.cov<-adhE.cov[adhE.cov$gene_id!="Seq18_adhE",]

# import metadata about adhE genes and the genomes they came from
adhE.meta<-read.table("adhE_gene_genome_metadata.txt",sep="\t")
names(adhE.meta) # look at column headers
colnames(adhE.meta)[which(names(adhE.meta) == "V1")] <- "gene_id"
colnames(adhE.meta)[which(names(adhE.meta) == "V2")] <- "genus"
colnames(adhE.meta)[which(names(adhE.meta) == "V3")] <- "species"

head(adhE.cov)

# then import the adhE coverage in propensity-matched HC mgms
adhE.HC.cov<-read.table("ABS_HCs_adhE_bowtie2_Coverage_Results.tsv",header=TRUE,sep="\t")
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
write.table(adhE.subset,file="ABS_and_Matched_HCs_adhE_bowtie2_coverage_with_metadata_updated.tsv",sep="\t",col.names=TRUE,row.names=FALSE)
## NOTE: saving subsetted table but also want to visualize with all adhE results; be mindful of the dfs you're referring to below

# more info on breath vs depth of coverage here:
## https://www.metagenomics.wiki/pdf/qc/coverage-read-depth
## https://www.metagenomics.wiki/tools/samtools/breadth-of-coverage

#### Visualize adhE by Taxa & Sample ####

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

# visualize with stacked barplot
adhe1<-ggplot(adhE.melt2, aes(x=SampleID, y=meandepth, fill=genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),axis.text = element_text(size=22),
        axis.text.x = element_text(hjust=1), legend.title = element_text(size=22,hjust=0),legend.text = element_text(size=22),
        plot.title = element_text(size=25),plot.subtitle = element_text(size=18))+
  guides(fill=guide_legend(ncol=1)) + 
  labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus")+coord_flip()

ggsave(adhe1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_barplot.png", width = 6000,height = 7000,units = "px",dpi = 200,create.dir = TRUE)

adhe1a<-ggplot(adhE.melt2, aes(x=SampleID, y=meandepth, fill=interaction(genus,species,sep="_")))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),axis.text = element_text(size=22),
        axis.text.x = element_text(hjust=1), legend.title = element_text(size=22,hjust=0),legend.text = element_text(size=22),
        plot.title = element_text(size=25),plot.subtitle = element_text(size=18))+
  guides(fill=guide_legend(ncol=1)) + 
  labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus")+coord_flip()

ggsave(adhe1a,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_species_barplot.png", width = 6000,height = 7000,units = "px",dpi = 200,create.dir = TRUE)

# visualize with heatmap
adhe1b<-ggplot(adhE.melt2, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=100,breaks=c(0,50,100,150,200),
                       labels=c("0","50","100","150","200"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=25),axis.text = element_text(size=22),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=25),
                        legend.text = element_text(size=25),plot.title = element_text(size=25), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=22, face = "bold"), plot.subtitle = element_text(size=22), strip.background=element_rect(fill="white")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(adhe1b,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_heatmap.png", width = 6000,height = 7500,units = "px",dpi = 200,create.dir = TRUE)

adhe1c<-ggplot(adhE.melt2, aes(SampleID, interaction(genus,species,sep="_"), fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=100,breaks=c(0,50,100,150,200),
                       labels=c("0","50","100","150","200"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=25),axis.text = element_text(size=22),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=25),
                        legend.text = element_text(size=25),plot.title = element_text(size=25), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=22, face = "bold"), plot.subtitle = element_text(size=22), strip.background=element_rect(fill="white")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(adhe1c,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_coverage_bowtie2_all_genus_species_heatmap.png", width = 6000,height = 7500,units = "px",dpi = 200,create.dir = TRUE)

#### Visualize Higher Coverage adhE Results ####

# subset out samples with coverage over 50%, meandepth > 5
adhE.subset<-adhE.all[ which( adhE.all$coverage > 50 & adhE.all$meandepth > 5) , ]
# note: seem to lose HC samples when we do this filtering becuase they had such low depth of coverage for adhE genes...
min(adhE.subset$meandepth) # sanity check that our filtering work...
min(adhE.subset$coverage) # sanity check that our filtering work...

# create long (dcast) version of adhE.all
adhE.subset.table<-reshape2::dcast(adhE.subset, SampleID ~ gene_id, value.var = "meandepth",sum)
rownames(adhE.subset.table)<-adhE.subset.table$SampleID

# melt down long version of adhE data for visualizing
adhE.subset.melt<-reshape2::melt(adhE.subset.table,by="SampleID")
head(adhE.subset.melt)
# rename column names for merging
colnames(adhE.subset.melt)[which(names(adhE.subset.melt) == "variable")] <- "gene_id"
colnames(adhE.subset.melt)[which(names(adhE.subset.melt) == "value")] <- "meandepth"
head(adhE.subset.melt)

# drop the zeros that came from creating subset table
adhE.subset.melt<-adhE.subset.melt[adhE.subset.melt$meandepth!=0,]
head(adhE.subset.melt)

# merge subsetted adhE cov data with adhE metadata
adhE.subset.melt2<-merge(adhE.meta,adhE.subset.melt,by="gene_id")
head(adhE.subset.melt2)

# stacked barplot
adhe2<-ggplot(adhE.subset.melt2, aes(x=SampleID, y=meandepth, fill=genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),axis.text = element_text(size=22),
        axis.text.x = element_text(hjust=1), legend.title = element_text(size=22,hjust=0),legend.text = element_text(size=22),
        plot.title = element_text(size=25),plot.subtitle = element_text(size=18))+
  guides(fill=guide_legend(ncol=1)) + 
  labs(x="Sample ID", y="Mean Depth of Coverage", title="adhE mean depth of coverage",fill="Genus",subtitle="Contains adhE genes with > 50% breadth of coverage and > 5 mean depth of coverage")+coord_flip()

ggsave(adhe2,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_higher_coverage_bowtie2_subset_barplot.png", width = 3000,height = 4000,units = "px",dpi = 200,create.dir = TRUE)

# heatmap
adhe2b<-ggplot(adhE.subset.melt2, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=100,breaks=c(0,50,100,150,200),
                       labels=c("0","50","100","150","200"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=25),axis.text = element_text(size=22),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=25),
                        legend.text = element_text(size=25),plot.title = element_text(size=25), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=22, face = "bold"), plot.subtitle = element_text(size=22), strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains adhE genes with > 50% breadth of coverage and > 5 mean depth of coverage") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+coord_flip()

ggsave(adhe2b,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_MGMs_adhE_higher_coverage_bowtie2_subset_heatmap.png", width = 4000,height = 5000,units = "px",dpi = 200,create.dir = TRUE)

#### Import FMT Metadata ####

fmt.metadata<-as.data.frame(read_xlsx("/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/data/ABS_FMT_Metadata_CLH.xlsx", sheet="Sheet1"))
fmt.metadata$SampleID<-gsub("_L00[0-9]_R1_001.fastq.gz","",fmt.metadata$filename)
unique(fmt.metadata$SampleID)

sample.order<-unique(fmt.metadata$SampleID)

# which FMT patient/HHP/Donor samples are found in the adhE data that had higher coverage?
fmt.metadata$SampleID[which(unique(fmt.metadata$SampleID) %in% unique(adhE.subset$SampleID))]
# ^ this doesn't include all samples, so let's visualize with all samples first

# view the df containing all adhE cov data used for visualization above...
head(adhE.melt2)

# subset out only FMT/FMT-HHP/Donor adhE cov data
adhE.fmt<-adhE.melt2[which(adhE.melt2$SampleID %in% fmt.metadata$SampleID),]
unique(adhE.fmt$SampleID) %in% unique(fmt.metadata$SampleID) # sanity check, should all be TRUE

# drop rows with meandepth=0
adhE.fmt<-adhE.fmt[adhE.fmt$meandepth>0,]
min(adhE.fmt$meandepth) # sanity check

# merge FMT/FMT-HHP/Donor adhE cov data with FMT metadata
adhE.fmt.meta<-merge(adhE.fmt,fmt.metadata,by="SampleID")

# turn adhE.fmt.meta$FMT_Label & Condition columns into factors for ordering
unique(adhE.fmt.meta$FMT_Label)
adhE.fmt.meta$FMT_Label<-factor(adhE.fmt.meta$FMT_Label,levels=c("Pre-FMT","Pre-FMT1 Abx",
                                                                 "Post-FMT1","Pre-FMT2 Abx",
                                                                 "Post-FMT2","Donor"))
unique(adhE.fmt.meta$Condition)
adhE.fmt.meta$Condition<-factor(adhE.fmt.meta$Condition,levels=c("Flare","Remission","HHP","Donor"))

adhE.fmt.meta$SampleID<-factor(adhE.fmt.meta$SampleID,levels=c(sample.order))

adhE.fmt.meta<-adhE.fmt.meta[order(adhE.fmt.meta$SampleID),]

# stacked barplot
adhe.fmt1<-ggplot(adhE.fmt.meta, aes(x=forcats::fct_rev(SampleID), y=meandepth, fill=genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  theme(axis.title.x = element_text(size=24),axis.title.y = element_text(size=24),axis.text = element_text(size=22),
        axis.text.x = element_text(hjust=1), legend.title = element_text(size=22,hjust=0),legend.text = element_text(size=22),
        plot.title = element_text(size=25),plot.subtitle = element_text(size=18),strip.text.y = element_text(size=18, face = "bold"),
        strip.background=element_rect(fill="ghostwhite"))+
  guides(fill=guide_legend(ncol=1)) + 
  labs(x="Sample ID", y="Mean Depth of Coverage", title="Contains only FMT Patient + FMT HHP + Donor Data",fill="Genus")+
  facet_grid(Condition~., scales = "free_y",space="free_y")+coord_flip()

ggsave(adhe.fmt1,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_barplot.png", width = 3000,height = 4000,units = "px",dpi = 200,create.dir = TRUE)

# heatmap
adhe.fmt2<-ggplot(adhE.fmt.meta, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
                       labels=c("0","50","100","150"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=25),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=28),
                        legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=25, face = "bold"), plot.subtitle = element_text(size=25), strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_grid(Condition~., scales = "free_y",space="free_y")+coord_flip()

ggsave(adhe.fmt2,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_heatmap.png", width = 5000,height = 6000,units = "px",dpi = 200,create.dir = TRUE)

# heatmap 2
adhe.fmt2b<-ggplot(adhE.fmt.meta, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
                       labels=c("0","50","100","150"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=25),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),
                        legend.title = element_text(hjust=0.5,size=28),
                        legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=25, face = "bold"), plot.subtitle = element_text(size=25), strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_grid(.~FMT_Label, scales = "free_y",space="free_y")+coord_flip()

ggsave(adhe.fmt2b,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_heatmap2.png", width = 5000,height = 6000,units = "px",dpi = 200,create.dir = TRUE)

# heatmap 3
adhe.fmt2c<-ggplot(adhE.fmt.meta, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
                       labels=c("0","50","100","150"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=25),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),axis.text.y = element_text(face="italic"),
                        legend.title = element_text(hjust=0.5,size=28),
                        legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=25, face = "bold"), plot.subtitle = element_text(size=25), strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_wrap(Condition~FMT_Label, scales = "free")

ggsave(adhe.fmt2c,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_heatmap3.png", width = 8000,height = 6000,units = "px",dpi = 200,create.dir = TRUE)

# heatmap 4
adhe.fmt2d<-ggplot(adhE.fmt.meta, aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
                       labels=c("0","50","100","150"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=25),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),axis.text.y = element_text(face="italic"),
                        legend.title = element_text(hjust=0.5,size=28),
                        legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=25, face = "bold"), plot.subtitle = element_text(size=25), strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_grid(Condition~FMT_Label, scales = "free",space="free")

ggsave(adhe.fmt2d,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_heatmap4.png", width = 9000,height = 8000,units = "px",dpi = 200,create.dir = TRUE)

# heatmap without HHP
# heatmap 5
adhe.fmt2e<-ggplot(adhE.fmt.meta[adhE.fmt.meta$Condition!="HHP",], aes(SampleID, genus, fill= meandepth)) +geom_tile(color="white")+
  scale_fill_gradient2(low="lightgray",mid="red1",high="purple3",midpoint=75,breaks=c(0,50,100,148),
                       labels=c("0","50","100","150"))+
  theme_classic()+theme(axis.title.x = element_blank(),axis.title.y = element_text(size=28),axis.text = element_text(size=25),
                        axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),axis.text.y = element_text(face="italic"),
                        legend.title = element_text(hjust=0.5,size=28),
                        legend.text = element_text(size=28),plot.title = element_text(size=28), strip.text.x = element_text(size = 20, face = "bold"),
                        strip.text.y = element_text(size=25, face = "bold"), plot.subtitle = element_text(size=25), strip.background=element_rect(fill="ghostwhite")) +
  labs(x="Sample ID", y="Genus", title="adhE mean depth of coverage",fill="Mean Depth of Coverage",subtitle="Contains only FMT Patient + FMT HHP + Donor Data") +
  guides(fill = guide_colourbar(barwidth = 0.7,barheight = 10))+scale_x_discrete(expand = c(0,0))+scale_y_discrete(expand=c(0,0))+
  facet_grid(Condition~FMT_Label, scales = "free",space="free")

ggsave(adhe.fmt2e,filename = "/Volumes/HsuLab/LF_HsuLabProjects/ABS_new/figures/adhE_MGMs/ABS_FMT_MGMs_adhE_coverage_bowtie2_NoHHP_heatmap5.png", width = 6000,height = 7000,units = "px",dpi = 200,create.dir = TRUE)

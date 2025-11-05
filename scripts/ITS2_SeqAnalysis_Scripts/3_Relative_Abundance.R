#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
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
  library(ggbiplot)
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
  library(NatParksPalettes)
  library(shades)
  library(ALDEx2)
  library(rstatix)
  library(decontam)
  library(ggvegan)
  library(microbiome)
  library(pairwiseAdonis)
  library(corrplot)
  library(fst)
  library(plotly)
  library(htmlwidgets)
  library(MoMAColors)
  library(microshades)
  library(lmtest)
  library(forcats)
  library(Maaslin2)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/ITS2/ITS2_DADA2/ABS_ITS2_DataReady.Rdata") # save global env to Rdata file

its2.ASV_table[1:4,1:4]
its2.ASV_table[(nrow(its2.ASV_table)-4):(nrow(its2.ASV_table)),(ncol(its2.ASV_table)-4):(ncol(its2.ASV_table))] # last 4 rows & cols
head(metadata)

# drop Donor sample
metadata.v2<-metadata[!grepl("DonorSamp", metadata$SampleID),]
head(metadata.v2)

#### List of Important R Objects for Reference ####

its2.ASV_table[1:5,1:5] 
# ^ this is the updated ASV table where singletons have been removed
## sample names were also updated

head(metadata) # updated metadata that includes Flare group colors
head(metadata.v2) # updated metadata with Flare group colors & FMT donor sample dropped
its2.ASV_meta[1:5,1:10] # combined metadata + ASV counts

head(its2.ASV_tax.clean) # updated taxa table with singleton ASVs removed
head(its2.AllTax.melt) # melted ASVs + taxonomy table -- for relative abundance calcs

head(its2.ALL.dat) # ALL data combined
# ^ metadata, ASVs + counts, & taxonomy data

#### Phyla Relative Abundance ####

# use dcast to count up ASVs within each Phylum across all of the samples
f.phyla_counts <- as.data.frame(dcast(its2.AllTax.melt, SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(f.phyla_counts) # counts by phyla per sample

# drop unknown & unclassified phyla
f.phyla_counts<-subset(f.phyla_counts,select=-c(`Fungi phy Incertae sedis`,Unknown))
head(f.phyla_counts) # counts by phyla per sample

dim(f.phyla_counts)

rownames(f.phyla_counts)<-f.phyla_counts$SampleID
dim(f.phyla_counts)
f.phyla_counts<-f.phyla_counts[,colSums(f.phyla_counts[,-1])>0] # drop phyla that are not represented
dim(f.phyla_counts) # sanity check that we dropped taxa with no hits

f.phyla_RelAb<-data.frame(decostand(f.phyla_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.phyla_RelAb) # sanity check to make sure the transformation worked!

f.phyla_RelAb$SampleID<-rownames(f.phyla_RelAb)
head(f.phyla_RelAb)
#write.csv(f.phyla_RelAb,"16S_Phyla_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# create a df with relative abundance data for phyla where each phyla is a column, plus the metadata
f.phyla_RA_meta<-merge(f.phyla_RelAb,metadata.v2, by="SampleID")

## now to create df with melted relative abundance data + metadata...
# melt down relativized data to merge with metadata.v2
f.phyla_m<-melt(f.phyla_RelAb)

head(f.phyla_m)
colnames(f.phyla_m)[which(names(f.phyla_m) == "variable")] <- "Phylum"
colnames(f.phyla_m)[which(names(f.phyla_m) == "value")] <- "Count"
head(f.phyla_m) ## relative abundance based on sum of counts by phyla!

f.phyla_melt_meta<-merge(f.phyla_m,metadata.v2, by="SampleID")
head(f.phyla_melt_meta) ## relative abundance based on sum of counts by phyla!
#f.phyla_melt_meta$SampleID = factor(f.phyla_melt_meta$SampleID, levels=unique(f.phyla_melt_meta$SampleID[order(f.phyla_melt_meta$Flare,f.phyla_melt_meta$Seas_Coll_Year)]), ordered=TRUE)

# Order phyla from most to least abundant for plots
phyla.list<-names(sort(colSums(f.phyla_RelAb[,-ncol(f.phyla_RelAb)]),decreasing=TRUE)) # see most abundant phyla
f.phyla_melt_meta$Phylum<-factor(f.phyla_melt_meta$Phylum,levels=phyla.list,ordered=TRUE)
f.phyla_melt_meta.v2<-f.phyla_melt_meta[order(f.phyla_melt_meta$Phylum,f.phyla_melt_meta$Count),]

# Barplot by SampleID
p.b1<-ggplot(f.phyla_melt_meta.v2, aes(x=fct_inorder(SampleID), y=Count, fill=Phylum))+geom_bar(stat="identity")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance", x="SampleID", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1))+scale_fill_manual(values=natparks.pals("Triglav", 11))+
  facet_wrap(vars(Flare), scales = "free")

ggsave(p.b1,filename = "figures/RelativeAbundance/Phylum/ABS_ITS2_Phyla.RA_barplot_mine.png", width=14, height=8, dpi=600,create.dir = TRUE)

p.b2<-ggplot(f.phyla_melt_meta.v2, aes(x = fct_inorder(SampleID), y = Count, fill = Phylum)) +
  geom_col(position = "fill")+
  facet_grid(~Flare, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance", labels = scales::percent)+
  scale_fill_manual(values=natparks.pals("Triglav", 11))+
  theme(text = element_text(colour = "black", size=11),
        axis.text=element_text(colour="black", size=11),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=7,hjust=1,angle=45))

ggsave(p.b2,filename = "figures/RelativeAbundance/Phylum/ABS_ITS2_Phyla.RA_barplot_withIDs.png", width = 2000,height = 730,units = "px",dpi = 200,create.dir = TRUE)

p.b2b<-ggplot(f.phyla_melt_meta.v2, aes(x = fct_inorder(SampleID), y = Count, fill = Phylum)) +
  geom_col(position = "fill")+
  facet_grid(~Flare, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance", labels = scales::percent)+
  scale_fill_manual(values=natparks.pals("Triglav", 11))+
  theme(text = element_text(colour = "black", size=11),
        axis.text=element_text(colour="black", size=11),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_blank())

ggsave(p.b2b,filename = "figures/RelativeAbundance/Phylum/ABS_ITS2_Phyla.RA_barplot_NoIDs.png", width = 1800,height = 720,units = "px",dpi = 200,create.dir = TRUE)

# Heatmap by SampleID
p.h1<-ggplot(f.phyla_melt_meta.v2, aes(SampleID, Phylum, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.45)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Fungal Phyla", title="Fungal Phyla & Flare Group",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ Flare, scales = "free")

ggsave(p.h1,filename = "figures/RelativeAbundance/Phylum/ABS_ITS2_Phyla.RA_heatmap.png", width=16, height=10, dpi=600,create.dir = TRUE)

#### Create Color Palette for Fungi ####

phyla.list # list of fungal phyla from most to least abundant in data
melt(as.character(natparks.pals("Triglav", 11))) # gives us the hex codes of the colors used in this palette!
fun.cols<-data.frame(melt(as.character(natparks.pals("Triglav", 11)))) # hex codes for color palette we are using for fungal phyla

fun.cols.list<-cbind(phyla.list,fun.cols)
fun.cols.list
colnames(fun.cols.list)[which(names(fun.cols.list) == "phyla.list")] <- "Phylum"
colnames(fun.cols.list)[which(names(fun.cols.list) == "value")] <- "FungalPhy_Color"
fun.cols.list

# ^^ will use this later when we drop more abudant and absent phyla but want to ensure the phyla colors do not change in the plot! 

#### Phyla Relative Abundance with Donor ####

# use dcast to count up ASVs within each Phylum across all of the samples
f.phyla_counts <- as.data.frame(dcast(its2.AllTax.melt, SampleID~Phylum, value.var="Count", fun.aggregate=sum)) ###
head(f.phyla_counts) # counts by phyla per sample

# drop unknown & unclassified phyla
f.phyla_counts<-subset(f.phyla_counts,select=-c(`Fungi phy Incertae sedis`,Unknown))
head(f.phyla_counts) # counts by phyla per sample

dim(f.phyla_counts)

rownames(f.phyla_counts)<-f.phyla_counts$SampleID
dim(f.phyla_counts)
f.phyla_counts<-f.phyla_counts[,colSums(f.phyla_counts[,-1])>0] # drop phyla that are not represented
dim(f.phyla_counts) # sanity check that we dropped taxa with no hits

f.phyla_RelAb<-data.frame(decostand(f.phyla_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.phyla_RelAb) # sanity check to make sure the transformation worked!

f.phyla_RelAb$SampleID<-rownames(f.phyla_RelAb)
head(f.phyla_RelAb)
#write.csv(f.phyla_RelAb,"16S_Phyla_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# create a df with relative abundance data for phyla where each phyla is a column, plus the metadata
f.phyla_RA_meta<-merge(f.phyla_RelAb,metadata, by="SampleID")

## now to create df with melted relative abundance data + metadata...
# melt down relativized data to merge with metadata
f.phyla_m<-melt(f.phyla_RelAb)

head(f.phyla_m)
colnames(f.phyla_m)[which(names(f.phyla_m) == "variable")] <- "Phylum"
colnames(f.phyla_m)[which(names(f.phyla_m) == "value")] <- "Count"
head(f.phyla_m) ## relative abundance based on sum of counts by phyla!

f.phyla_melt_meta<-merge(f.phyla_m,metadata, by="SampleID")
head(f.phyla_melt_meta) ## relative abundance based on sum of counts by phyla!
#f.phyla_melt_meta$SampleID = factor(f.phyla_melt_meta$SampleID, levels=unique(f.phyla_melt_meta$SampleID[order(f.phyla_melt_meta$Flare,f.phyla_melt_meta$Seas_Coll_Year)]), ordered=TRUE)

# subset specific patients and order samples by date
subset.patients<-c("A020", "A021", "Donor")
f.phyla_patient.subset<-f.phyla_melt_meta[grepl(paste(subset.patients,collapse="|"),f.phyla_melt_meta$Patient),]

# merge subsetted ASV rel ab data with the fungal phyla color palette we made earlier
fun.cols.list
f.phyla_patient.sbst.cols<-merge(f.phyla_patient.subset,fun.cols.list,by.x="Phylum",by.y="Phylum")

# Order dates chronologically
f.phyla_patient.sbst.cols$Date<-factor(f.phyla_patient.sbst.cols$Date,levels=c("10/14/20","10/15/20","11/10/22","5/9/23","5/12/23",
                                                                         "5/22/23","8/6/23","8/18/23","11/15/23","11/16/23",
                                                                         ""),ordered=TRUE)
# make Dx into a factor so it has specific order
f.phyla_patient.sbst.cols$Dx<-factor(f.phyla_patient.sbst.cols$Dx,levels=c("HHP","ABS","Donor"))

# order df by Date then Dx so the samples appear in chronological order and subsetted in the right order
f.phyla_patient.sbst.cols<-f.phyla_patient.sbst.cols[order(f.phyla_patient.sbst.cols$Date,f.phyla_patient.sbst.cols$Dx),]

# drop counts that are > 0, < 1, and associated with phylum Ascomycota
## this is so that we can see the rarer fungal taxa more clearly across time
## if we do not specifically drop Ascomycota hits, the phyla will appear in the legend even though it's not shown on the actual stacked bar plot
#f.phyla_patient.sbst.cols.rare<-subset(f.phyla_patient.sbst.cols, Count > 0.000000 & Count < 1 & Phylum!="Ascomycota")
f.phyla_patient.sbst.cols.rare<-subset(f.phyla_patient.sbst.cols, Count > 0.000000)

# Barplot by Date & Dx
p.b1<-ggplot(f.phyla_patient.sbst.cols.rare, aes(x=fct_inorder(Date), y=Count, fill=Phylum))+geom_bar(stat="identity")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Phylum Relative Abundance", x="", y="Relative Abundance", fill="Phylum")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=90),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(name = "Relative Abundance", labels = scales::percent,limits = c(0,0.15))+
  scale_fill_manual(name ="Phylum",values=unique(f.phyla_patient.sbst.cols.rare$FungalPhy_Color[order(f.phyla_patient.sbst.cols.rare$Phylum)]))+
  facet_grid(~Dx, scales = "free_x", space = "free_x") + rremove("x.ticks")

ggsave(p.b1,filename = "figures/RelativeAbundance_Donor/Phylum/ABS_ITS2_Phyla.RA_Subset_with_Donor_Dx_barplot_mine.png", width=14, height=8, dpi=600,create.dir = TRUE)

p.b2<-ggplot(f.phyla_patient.sbst.cols.rare, aes(x = fct_inorder(Date), y = Count, fill = Phylum)) +
  geom_col(position = "stack")+ # << changed position to "stack" from "fill" because otherwise it uses the scale_y_continuous limits as a way to rule out other taxa
  facet_grid(~Dx, scales = "free_x", space = "free_x")+
  theme_classic()+
  rremove("x.ticks")+
  scale_y_continuous(name = "Relative Abundance",breaks=c(0,0.05,0.10,0.15,0.90,0.95,1),labels = scales::percent)+
  scale_fill_manual(name ="Phylum",values=unique(f.phyla_patient.sbst.cols.rare$FungalPhy_Color[order(f.phyla_patient.sbst.cols.rare$Phylum)]))+
  theme(text = element_text(colour = "black", size=11),
        axis.text=element_text(colour="black", size=11),
        axis.title.x = element_blank(),
        legend.title = element_blank(),
        axis.text.x = element_text(size=7,hjust=1,angle=90)) + coord_cartesian(ylim = c(0,0.15))

ggsave(p.b2,filename = "figures/RelativeAbundance_Donor/Phylum/ABS_ITS2_Phyla.RA_Subset_with_Donor_Dx_barplot_withDates_Ascomycota.png", width = 2000,height = 730,units = "px",dpi = 200,create.dir = TRUE)

# p.b2b<-ggplot(f.phyla_patient.subset, aes(x = fct_inorder(Date), y = Count, fill = Phylum)) +
#   geom_col(position = "fill")+
#   facet_grid(~Flare, scales = "free_x", space = "free_x")+
#   theme_classic()+
#   rremove("x.ticks")+
#   scale_y_continuous(name = "Relative Abundance", labels = scales::percent)+
#   scale_fill_manual(values=natparks.pals("Triglav", 11))+
#   theme(text = element_text(colour = "black", size=11),
#         axis.text=element_text(colour="black", size=11),
#         axis.title.x = element_blank(),
#         legend.title = element_blank(),
#         axis.text.x = element_blank())
# 
# ggsave(p.b2b,filename = "figures/RelativeAbundance_Donor/Phylum/ABS_ITS2_Phyla.RA_Donor_barplot_NoIDs.png", width = 1800,height = 720,units = "px",dpi = 200,create.dir = TRUE)

# Heatmap by SampleID
p.h1<-ggplot(f.phyla_melt_meta.v2, aes(SampleID, Phylum, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.45)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Fungal Phyla", title="Fungal Phyla & Flare Group",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ Flare, scales = "free")

ggsave(p.h1,filename = "figures/RelativeAbundance_Donor/Phylum/ABS_ITS2_Phyla.RA_heatmap.png", width=16, height=10, dpi=600,create.dir = TRUE)

#### Class Relative Abundance ####

# use dcast to count up ASVs within each Class across all of the samples
f.class_counts <- as.data.frame(dcast(its2.AllTax.melt, SampleID~Class, value.var="Count", fun.aggregate=sum)) ###
head(f.class_counts) # counts by class per sample
dim(f.class_counts)
rownames(f.class_counts)<-f.class_counts$SampleID
f.class_counts<-subset(f.class_counts, select=-c(Unknown))
dim(f.class_counts)
f.class_counts<-f.class_counts[,colSums(f.class_counts[,-1])>0] # drop classes that are not represented
dim(f.class_counts) # sanity check that we dropped taxa with no hits

f.class_RelAb<-data.frame(decostand(f.class_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.class_RelAb) # sanity check to make sure the transformation worked!

f.class_RelAb$SampleID<-rownames(f.class_RelAb)
head(f.class_RelAb)
#write.csv(f.class_RelAb,"16S_class_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
f.class_m<-melt(f.class_RelAb)

head(f.class_m)
colnames(f.class_m)[which(names(f.class_m) == "variable")] <- "Class"
colnames(f.class_m)[which(names(f.class_m) == "value")] <- "Count"
head(f.class_m) ## relative abundance based on sum of counts by class!

f.class_RA_meta<-merge(f.class_m,metadata, by="SampleID")
head(f.class_RA_meta) ## relative abundance based on sum of counts by class!
f.class_RA_meta$SampleID = factor(f.class_RA_meta$SampleID, levels=unique(f.class_RA_meta$SampleID[order(f.class_RA_meta$Flare,f.class_RA_meta$Seas_Coll_Year)]), ordered=TRUE)

# Barplot by SampleID

c.b1<-ggplot(f.class_RA_meta, aes(x=SampleID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Classes", x="SampleID", y="Relative Abundance", fill="Class")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=4)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(c.b1,filename = "figures/RelativeAbundance/Class/ABS_ITS2_Class.RA_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

c.b2<-ggplot(f.class_RA_meta[f.class_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Taxa with Relative Abundance > 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(c.b2,filename = "figures/RelativeAbundance/Class/ABS_ITS2_Class.RA_1perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

c.b3<-ggplot(f.class_RA_meta[f.class_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Class))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Taxa with Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(c.b3,filename = "figures/RelativeAbundance/Class/ABS_ITS2_Class.RA_5perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

head(f.class_RA_meta)

# Heatmap by SampleID
c.h1<-ggplot(f.class_RA_meta, aes(SampleID, Class, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Fungal Class", title="Fungal Class & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ SampDate,scales="free")

ggsave(c.h1,filename = "figures/RelativeAbundance/Class/ABS_ITS2_class.RA_heatmap.png", width=16, height=15, dpi=600,create.dir = TRUE)

its2.AllTax.melt[1:4,1:4]


#### Order Relative Abundance ####

# use dcast to count up ASVs within each Order across all of the samples
f.ord_counts <- as.data.frame(dcast(its2.AllTax.melt, SampleID~Order, value.var="Count", fun.aggregate=sum)) ###
head(f.ord_counts) # counts by class per sample
dim(f.ord_counts)
rownames(f.ord_counts)<-f.ord_counts$SampleID
f.ord_counts<-subset(f.ord_counts, select=-c(Unknown))
dim(f.ord_counts)
f.ord_counts<-f.ord_counts[,colSums(f.ord_counts[,-1])>0] # drop orders that are not represented
dim(f.ord_counts) # sanity check that we dropped taxa with no hits

f.ord_RelAb<-data.frame(decostand(f.ord_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.ord_RelAb) # sanity check to make sure the transformation worked!

f.ord_RelAb$SampleID<-rownames(f.ord_RelAb)
head(f.ord_RelAb)
#write.csv(f.ord_RelAb,"16S_class_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
f.ord_m<-melt(f.ord_RelAb)

head(f.ord_m)
colnames(f.ord_m)[which(names(f.ord_m) == "variable")] <- "Order"
colnames(f.ord_m)[which(names(f.ord_m) == "value")] <- "Count"
head(f.ord_m) ## relative abundance based on sum of counts by order!

f.ord_RA_meta<-merge(f.ord_m,metadata, by="SampleID")
head(f.ord_RA_meta) ## relative abundance based on sum of counts by order!
f.ord_RA_meta$SampleID = factor(f.ord_RA_meta$SampleID, levels=unique(f.ord_RA_meta$SampleID[order(f.ord_RA_meta$Flare,f.ord_RA_meta$Seas_Coll_Year)]), ordered=TRUE)

# Barplot by SampleID

o.b1<-ggplot(f.ord_RA_meta, aes(x=SampleID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Orders", x="SampleID", y="Relative Abundance", fill="Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=4)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(o.b1,filename = "figures/RelativeAbundance/Order/ABS_ITS2_Order.RA_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

o.b2<-ggplot(f.ord_RA_meta[f.ord_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Orders", x="SampleID", y="Relative Abundance", fill="Order",subtitle="Only Taxa with Relative Abundance > 1%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(o.b2,filename = "figures/RelativeAbundance/Order/ABS_ITS2_Order.RA_1perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

o.b3<-ggplot(f.ord_RA_meta[f.ord_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Orders", x="SampleID", y="Relative Abundance", fill="Order",subtitle="Only Taxa with Relative Abundance > 5%")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(o.b3,filename = "figures/RelativeAbundance/Order/ABS_ITS2_Order.RA_5perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

burkholderiales.only<-ggplot(f.ord_RA_meta[f.ord_RA_meta$Order=="Burkholderiales",], aes(x=SampleID, y=Count, fill=Order))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Burkholderiales", x="SampleID", y="Relative Abundance", fill="Order",subtitle="Only Burkholderiales Order")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(burkholderiales.only,filename = "figures/RelativeAbundance/Order/ABS_ITS2_Order.RA_Burkholderiales_Only_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

head(f.ord_RA_meta)

# Heatmap by SampleID
o.h1<-ggplot(f.ord_RA_meta, aes(SampleID, Order, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Fungal Order", title="Fungal Order & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ SampDate,scales="free")

ggsave(o.h1,filename = "figures/RelativeAbundance/Order/ABS_ITS2_class.RA_heatmap.png", width=16, height=15, dpi=600,create.dir = TRUE)

its2.AllTax.melt[1:4,1:4]


#### Family Relative Abundance ####

# use dcast to count up ASVs within each Family across all of the samples
f.fam_counts <- as.data.frame(dcast(its2.AllTax.melt, SampleID~Family, value.var="Count", fun.aggregate=sum)) ###
head(f.fam_counts) # counts by fam per sample
dim(f.fam_counts)
rownames(f.fam_counts)<-f.fam_counts$SampleID
colnames(f.fam_counts)<-gsub(" ", ".",colnames(f.fam_counts))
f.fam_counts<-subset(f.fam_counts, select=-c(Unknown, Unknown.Family))
dim(f.fam_counts)
f.fam_counts<-f.fam_counts[,colSums(f.fam_counts[,-1])>0] # drop families that are not represented
dim(f.fam_counts) # sanity check that we dropped taxa with no hits

f.fam_RelAb<-data.frame(decostand(f.fam_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.fam_RelAb) # sanity check to make sure the transformation worked!

f.fam_RelAb$SampleID<-rownames(f.fam_RelAb)
head(f.fam_RelAb)
#write.csv(f.fam_RelAb,"16S_fam_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
f.fam_m<-melt(f.fam_RelAb)

head(f.fam_m)
colnames(f.fam_m)[which(names(f.fam_m) == "variable")] <- "Family"
colnames(f.fam_m)[which(names(f.fam_m) == "value")] <- "Count"
head(f.fam_m) ## relative abundance based on sum of counts by fam!
f.fam_m$Family<-gsub("^X.","",f.fam_m$Family) # get rid of leading X. in Family names
f.fam_m$Family<-gsub("\\.\\."," ",f.fam_m$Family) # get rid of .. in family name --> . is regex
f.fam_m$Family<-gsub("\\."," ",f.fam_m$Family) # get rid of . in family name --> . is regex

f.fam_RA_meta<-merge(f.fam_m,metadata, by="SampleID")
head(f.fam_RA_meta) ## relative abundance based on sum of counts by fam!
f.fam_RA_meta$SampleID = factor(f.fam_RA_meta$SampleID, levels=unique(f.fam_RA_meta$SampleID[order(f.fam_RA_meta$Flare,f.fam_RA_meta$Seas_Coll_Year)]), ordered=TRUE)

# Barplot by SampleID

f.b1<-ggplot(f.fam_RA_meta, aes(x=SampleID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Families", x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=3))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(f.b1,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_barplot.png", width=17, height=10, dpi=600,create.dir = TRUE)

f.b1a<-ggplot(f.fam_RA_meta[f.fam_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Families", subtitle="Only Taxa with Relative Abundance > 1%",x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=4))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(f.b1a,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_1perc_barplot.png", width=17, height=10, dpi=600,create.dir = TRUE)

f.b1b<-ggplot(f.fam_RA_meta[f.fam_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Families", subtitle="Only Taxa with Relative Abundance > 5%",x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(f.b1b,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_5perc_barplot.png", width=17, height=10, dpi=600,create.dir = TRUE)

f.b1c<-ggplot(f.fam_RA_meta[f.fam_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Fungal Families", subtitle="Only Taxa with Relative Abundance > 10%",x="SampleID", y="Relative Abundance", fill="Family")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+scale_y_continuous(expand = c(0,0),limits = c(0,1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(f.b1c,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_10perc_barplot.png", width=17, height=10, dpi=600,create.dir = TRUE)

head(f.fam_RA_meta)

# Heatmap by SampleID

f.h1<-ggplot(f.fam_RA_meta, aes(SampleID, Family, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3",mid="white",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Fungal Family", title="Fungal Families & Sample Type",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0)) +
  facet_grid(. ~ SampDate, scales="free")

ggsave(f.h1,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_heatmap.png", width=12, height=10, dpi=600,create.dir = TRUE)

# Taxonomic summary by Sample ID + Collection Period

f.ts1<-ggplot(f.fam_RA_meta[f.fam_RA_meta$Count>0.025,], aes(Family, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.fam_RA_meta$SampDate_Color[order(f.fam_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Families", y="Relative Abundance", title="Salton Sea Dust: Fungal Families & Sample Date",subtitle="Includes taxa with Relative Abundance > 2.5%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(f.ts1,filename = "figures/RelativeAbundance/Family/ABS_ITS2_Family.RA_2.5perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

f.ts2<-ggplot(f.fam_RA_meta[f.fam_RA_meta$Count>0.05,], aes(Family, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.fam_RA_meta$SampDate_Color[order(f.fam_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Families", y="Relative Abundance", title="Salton Sea Dust: Fungal Families & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(f.ts2,filename = "figures/RelativeAbundance/Family/ABS_ITS2_Family.RA_5perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

its2.AllTax.melt[1:4,1:4]


# # by Family + depth
# its2.fam.dep <- as.data.frame(dcast(its2.AllTax.melt,Depth_m~Family, value.var="Count", fun.aggregate=sum)) ###
# head(its2.fam.dep) # counts by Family + sample depe
# rownames(its2.fam.dep)<-its2.fam.dep$Depth_m
# colnames(its2.fam.dep)<-gsub(" ", ".",colnames(its2.fam.dep))
# its2.fam.dep<-subset(its2.fam.dep, select=-c(Unknown, Unknown.Family))
#
# f.RA_fam.dep<-data.frame(decostand(its2.fam.dep[,-1], method="total", MARGIN=1, na.rm=TRUE))
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(f.RA_fam.dep) # sanity check
# f.RA_fam.dep$Depth_m<-rownames(f.RA_fam.dep) # Depth_m is now a character, not a factor!
# head(f.RA_fam.dep)
#
# #melt down relativized data to merge with metadata
# f.fam.dep_m<-melt(f.RA_fam.dep, by="Depth_m")
#
# head(f.fam.dep_m)
# colnames(f.fam.dep_m)[which(names(f.fam.dep_m) == "variable")] <- "Family"
# colnames(f.fam.dep_m)[which(names(f.fam.dep_m) == "value")] <- "Count"
# head(f.fam.dep_m) ## relative abundance based on sum of counts by Family!
#
# #dep_meta<-unique(data.frame("Depth_m"=metadata$Depth_m, "Sample_Color"=metadata$Sample_Color))
# #p_dep_meta<-merge(dep_meta,f.fam.dep_m, by="Depth_m")
#
# # Barplot by Depth
#
# fd1<-ggplot(f.fam.dep_m, aes(x=Depth_m, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Fungal Classes", x="Depth (m)", y="Relative Abundance", fill="Class")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=3))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))
# #+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fd1,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_barplot_depth.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# fd1a<-ggplot(f.fam.dep_m[f.fam.dep_m$Count>0.05,], aes(x=Depth_m, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Fungal Classes", x="SampleID", y="Relative Abundance", fill="Class",subtitle="Only Relative Abundance > 5%")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))
# #+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fd1a,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_barplot_depth_5percent.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# # Taxonomic Summary by Depth
#
# fd2<-ggplot(f.fam.dep_m, aes(Family, Count)) +
#   geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
#   scale_colour_gradient2(high="blue3",low="red",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Family", y="Relative Abundance", title="Fungal Families & Depth",color="Depth (m)")
#
# ggsave(fd2,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_depth_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# fd2a<-ggplot(f.fam.dep_m[f.fam.dep_m$Count>0.05,], aes(Family, Count)) +
#   geom_jitter(aes(color=as.numeric(Depth_m)), size=2, width=0.15, height=0) +
#   scale_colour_gradient2(low="red",high="blue3",midpoint=5.5,guide = guide_colourbar(reverse = TRUE)) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Family", y="Relative Abundance", title="Fungal Families & Depth",color="Depth (m)")
#
# ggsave(fd2a,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_depth_5percent_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# # by Family + Sampling Date
# its2.fam.date <- as.data.frame(dcast(its2.AllTax.melt,SampDate~Family, value.var="Count", fun.aggregate=sum)) ###
# head(its2.fam.date) # counts by Family + sample depe
# rownames(its2.fam.date)<-its2.fam.date$SampDate
# colnames(its2.fam.date)<-gsub(" ", ".",colnames(its2.fam.date))
# its2.fam.date<-subset(its2.fam.date, select=-c(Unknown, Unknown.Family))
#
# f.RA_fam.date<-data.frame(decostand(its2.fam.date[,-1], method="total", MARGIN=1, na.rm=TRUE))
# # relative abundance of taxa data where everything is divided by margin total (default MARGIN = 1 = rows) -- rows = samples
# rowSums(f.RA_fam.date) # sanity check
# f.RA_fam.date$SampDate<-rownames(f.RA_fam.date)
# head(f.RA_fam.date)
#
# #melt down relativized data to merge with metadata
# f.fam.date_m<-melt(f.RA_fam.date, by="SampDate")
#
# head(f.fam.date_m)
# colnames(f.fam.date_m)[which(names(f.fam.date_m) == "variable")] <- "Family"
# colnames(f.fam.date_m)[which(names(f.fam.date_m) == "value")] <- "Count"
# head(f.fam.date_m) ## relative abundance based on sum of counts by Family!
#
# f.fam.date_m$SampDate<-factor(f.fam.date_m$SampDate, levels=c("August.2021","December.2021","April.2022"))
#
# # Barplot by Sample Date
#
# fsd1<-ggplot(f.fam.date_m, aes(x=SampDate, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Fungal Families", x="SampleID", y="Relative Abundance", fill="Family")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=3))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fsd1,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_barplot_sampdate.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# fsd1a<-ggplot(f.fam.date_m[f.fam.date_m$Count>0.01,], aes(x=SampDate, y=Count, fill=Family))+geom_bar(stat="identity",colour="black")+theme_classic()+
#   labs(title = "Relative Abundance of Fungal Families", x="SampleID", y="Relative Abundance", fill="Family",subtitle="Only Relative Abundance > 1%")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=2))+
#   scale_y_continuous(expand = c(0,0),limits = c(0,1))+ scale_x_discrete(labels=c("August.2021"="August 2021","December.2021"="December 2021","April.2022"="April 2022"))
#
# ggsave(fsd1a,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_barplot_sampdate_1percent.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# # Taxonomic Summary by Sample Date
#
# #f.fam.date_m2<-merge(f.fam.date_m, metadata, by="SampDate")
#
# colorset1 # remember which date goes with each color
#
# fsd2<-ggplot(f.fam.date_m, aes(Family, Count)) +
#   geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
#   scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Family", y="Relative Abundance", title="Fungal Families & Sample Date")
#
# ggsave(fsd2,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_date_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# fsd3<-ggplot(f.fam.date_m[f.fam.date_m$Count>0.1,], aes(Family, Count)) +
#   geom_jitter(aes(color=factor(SampDate)), size=2, width=0.15, height=0) +
#   scale_color_manual(name ="Sample Date", values=c("#ef781c","#03045e","#059c3f"), labels=c("August 2021","December 2021","April 2022")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Family", y="Relative Abundance", title="Fungal Families & Sample Date")
#
# ggsave(fsd3,filename = "figures/RelativeAbundance/Family/ABS_ITS2_fam.RA_date_taxasum_5percent.png", width=15, height=10, dpi=600,create.dir = TRUE)


#### Genus Relative Abundance ####
# use dcast to count up ASVs within each Genus across all of the samples
head(its2.AllTax.melt)
its2.AllTax.melt.g<-subset(its2.AllTax.melt, !Genus %in% c("Unknown","Fungi gen Incertae sedis")) # drop unknown genera so they don't skew analyses
"Unknown" %in% its2.AllTax.melt.g$Genus

f.genus_counts <- as.data.frame(dcast(its2.AllTax.melt.g, SampleID~Genus, value.var="Count", fun.aggregate=sum)) ###
head(f.genus_counts) # counts by genus per sample
dim(f.genus_counts)
rownames(f.genus_counts)<-f.genus_counts$SampleID
f.genus_counts[1:4,1:4]
f.genus_counts<-f.genus_counts[,colSums(f.genus_counts[,-1])>0] # drop classes that are not represented
dim(f.genus_counts) # sanity check that we dropped taxa with no hits

f.genus_RelAb<-data.frame(decostand(f.genus_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.genus_RelAb) # sanity check to make sure the transformation worked!

f.genus_RelAb$SampleID<-rownames(f.genus_RelAb)
head(f.genus_RelAb)
#write.csv(f.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case
names(f.genus_RelAb)<-gsub("\\."," ",names(f.genus_RelAb)) # replace "." with " " in genus names
head(f.genus_RelAb)

# melt down relativized data to merge with metadata
f.genus_m<-melt(f.genus_RelAb)

head(f.genus_m)
colnames(f.genus_m)[which(names(f.genus_m) == "variable")] <- "Genus"
colnames(f.genus_m)[which(names(f.genus_m) == "value")] <- "Count"
head(f.genus_m) ## relative abundance based on sum of counts by genus!
# f.genus_m$Genus<-gsub("^X.","",f.genus_m$Genus) # get rid of leading X. in Genus names
# f.genus_m$Genus<-gsub("\\.\\."," ",f.genus_m$Genus) # get rid of .. in species name --> . is regex
# f.genus_m$Genus<-gsub("\\."," ",f.genus_m$Genus) # get rid of . in species name --> . is regex
# f.genus_m$Genus<-gsub("_"," ",f.genus_m$Genus) #
head(f.genus_m) ## relative abundance based on sum of counts by genus!

f.genus_melt_meta<-merge(f.genus_m,metadata, by="SampleID")
head(f.genus_melt_meta) ## relative abundance based on sum of counts by genus!

# Order genera from most to least abundant for plots
genera.list<-names(sort(colSums(f.genus_RelAb[,-ncol(f.genus_RelAb)]),decreasing=TRUE)) # see most abundant phyla
f.genus_melt_meta$Genus<-factor(f.genus_melt_meta$Genus,levels=genera.list,ordered=TRUE)
f.genus_melt_meta.v2<-f.genus_melt_meta[order(f.genus_melt_meta$Genus,f.genus_melt_meta$Count),]

saveRDS(f.genus_melt_meta, file = "data/Amplicon/ABS_ITS2_GenusOnly_RelativeAbundance_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

#### Maaslin2 with Fungal Genera ####

head(f.genus_RelAb)
head(metadata.v2)

# drop Donor sample
f.genus_RelAb.v2<-f.genus_RelAb[!grepl("DonorSamp", rownames(f.genus_RelAb)),]
head(f.genus_RelAb.v2)
dim(f.genus_RelAb.v2) # should be 62 samples now that we dropped the Donor

# reorder metadata to match fungal genera rel.ab. df
## make sure they both have sample IDs as rows
metadata.v2=metadata.v2[rownames(f.genus_RelAb.v2),] ## will drop rows that are not shared by both dataframes!

# arcsin square root transformation on fungal genera relative abundances
## default method for model is linear regression
## if we add fixed & random effects, uses lme4 linear mixed-effect model: https://cran.r-project.org/web/packages/lme4/index.html
## more here: https://github.com/biobakery/biobakery/wiki/maaslin2#31-maaslin-2-input
fit_dataFlare  = Maaslin2(input_data = f.genus_RelAb.v2, # < input would be fungal genera relative abundances
                          input_metadata = metadata.v2,
                          min_abundance  = 0.00001,
                          min_prevalence = 0.2,
                          transform = "AST",      # < the addition of the AST
                          normalization  = "NONE",
                          output         = "data/ABS_ITS2_Maaslin2_Results",
                          fixed_effects  = c("Flare"),
                          random_effects = c("House"),
                          reference      = c('Flare,Flare'),
                          max_significance = 0.10,
                          standardize = FALSE)

# adding Patient as a random effect since patients were sampled twice...
# ^ results of adding Patient as a random effect were example the same
fit_dataFlare2  = Maaslin2(input_data = f.genus_RelAb.v2, # < input would be fungal genera relative abundances
                          input_metadata = metadata.v2,
                          min_abundance  = 0.00001,
                          min_prevalence = 0.2,
                          transform = "AST",      # < the addition of the AST
                          normalization  = "NONE",
                          output         = "data/ABS_ITS2_Maaslin2_Results_Patient",
                          fixed_effects  = c("Flare"),
                          random_effects = c("House","Patient"),
                          reference      = c('Flare,Flare'),
                          max_significance = 0.10,
                          standardize = FALSE)

# Maaslin2 Fixed Effects Notes:
## The following command runs MaAsLin 2 on the HMP2 data, running a multivariable regression model to test for the association between microbial species abundance versus IBD diagnosis and dysbiosis scores (fixed_effects = c("diagnosis", "dysbiosis")). 
## For any categorical variable with more than 2 levels, you will also have to specify which variable should be the reference level by using (reference = c("diagnosis,nonIBD")). 
## NOTE: adding a space between the variable and level might result in the wrong reference level being used. Output are generated in a folder called demo_output under the current working directory (output = "demo_output"). 
## In this case, the example input data has been pre-normalized and pre-filtered, so we turn off the default normalization and prevalence filtering as well.

# Maaslin2 Data output files
# significant_results.tsv
## full list of associations that pass MaAsLin 2's significance threshold, ordered by increasing q-values
# all_results.tsv
## Same format as significant_results.tsv, but include all association results (instead of just the significant ones).
## You can also access this table within R using fit_data$results.
# residuals.rds
## This file contains a data frame with residuals for each feature.
# maaslin2.log
## This file contains all log information for the run.
## It includes all settings, warnings, errors, and steps run.

# Visualization output files
# heatmap.pdf
## This file contains a heatmap of the significant associations.
# [a-z/0-9]+.pdf
## A plot is generated for each significant association.
## Scatter plots are used for continuous metadata.
## Box plots are for categorical data.
## Data points plotted are after normalization, filtering, and transform.

#### Genus + Species Relative Abundance ####
# use dcast to count up ASVs within each Genus across all of the samples
head(its2.AllTax.melt)
its2.AllTax.melt.g<-subset(its2.AllTax.melt, its2.AllTax.melt$Genus!="Unknown") # drop unknown genera so they don't skew analyses
"Unknown" %in% its2.AllTax.melt.g$Genus

f.gen.spec_counts <- as.data.frame(dcast(its2.AllTax.melt.g, SampleID~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
head(f.gen.spec_counts) # counts by gen.spec per sample
dim(f.gen.spec_counts)
rownames(f.gen.spec_counts)<-f.gen.spec_counts$SampleID
f.gen.spec_counts[1:4,1:4]
f.gen.spec_counts<-f.gen.spec_counts[,colSums(f.gen.spec_counts[,-1])>0] # drop classes that are not represented
dim(f.gen.spec_counts) # sanity check that we dropped taxa with no hits

f.gen.spec_RelAb<-data.frame(decostand(f.gen.spec_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.gen.spec_RelAb) # sanity check to make sure the transformation worked!

f.gen.spec_RelAb$SampleID<-rownames(f.gen.spec_RelAb)
head(f.gen.spec_RelAb)
#write.csv(f.gen.spec_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
f.gen.spec_m<-melt(f.gen.spec_RelAb)

head(f.gen.spec_m)
colnames(f.gen.spec_m)[which(names(f.gen.spec_m) == "variable")] <- "Genus_species"
colnames(f.gen.spec_m)[which(names(f.gen.spec_m) == "value")] <- "Count"
head(f.gen.spec_m) ## relative abundance based on sum of counts by gen.spec!
f.gen.spec_m$Genus_species<-gsub("^X.","",f.gen.spec_m$Genus_species) # get rid of leading X. in Genus_species names
f.gen.spec_m$Genus_species<-gsub("\\.\\."," ",f.gen.spec_m$Genus_species) # get rid of .. in species name --> . is regex
f.gen.spec_m$Genus_species<-gsub("\\."," ",f.gen.spec_m$Genus_species) # get rid of . in species name --> . is regex
f.gen.spec_m$Genus_species<-gsub("_"," ",f.gen.spec_m$Genus_species) #
head(f.gen.spec_m) ## relative abundance based on sum of counts by gen.spec!

f.gen.spec_RA_meta<-merge(f.gen.spec_m,metadata, by="SampleID")
head(f.gen.spec_RA_meta) ## relative abundance based on sum of counts by gen.spec!
max(f.gen.spec_RA_meta$Count)
f.gen.spec_RA_meta$SampleID = factor(f.gen.spec_RA_meta$SampleID, levels=unique(f.gen.spec_RA_meta$SampleID[order(f.gen.spec_RA_meta$Flare,f.gen.spec_RA_meta$Seas_Coll_Year)]), ordered=TRUE)
f.gen.spec_RA_meta$Sample_Type<-"Dust"

# separate genera RelAb data by site for downstream figs
WI.gen.spec.RA<-subset(f.gen.spec_RA_meta,f.gen.spec_RA_meta$Flare=="WI")
BDC.gen.spec.RA<-subset(f.gen.spec_RA_meta,f.gen.spec_RA_meta$Flare=="BDC")
PD.gen.spec.RA<-subset(f.gen.spec_RA_meta,f.gen.spec_RA_meta$Flare=="PD")
DP.gen.spec.RA<-subset(f.gen.spec_RA_meta,f.gen.spec_RA_meta$Flare=="DP")

saveRDS(f.gen.spec_RA_meta, file = "data/Amplicon/SSD_GenusSpecies_RelativeAbundance_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# Barplot by SampleID

f.gen_RAall<-ggplot(f.gen.spec_RA_meta, aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=10)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(f.gen_RAall,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

f.gen_RA0.0<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.005,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 0.05%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=10)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(f.gen_RA0.0,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

f.gen_RA0<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=5))

ggsave(f.gen_RA0,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_1perc_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

f.gen_RA0v2<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.gen_RA0v2,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_1perc_barplot_v2.png", width=30, height=25, dpi=600,create.dir = TRUE)

f.gen_RA01<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.02,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 2%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=5))

ggsave(f.gen_RA01,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_2perc_barplot.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.gen_RA01v2<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.02,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 2%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+
  facet_wrap(vars(Flare), scales = "free")

ggsave(f.gen_RA01v2,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_2perc_barplot_v2.png", width=30, height=25, dpi=600,create.dir = TRUE)

f.gen_RA1a<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))

ggsave(f.gen_RA1a,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_5perc_barplot_v1.png", width=20, height=10, dpi=600,create.dir = TRUE)

f.gen_RA1v2<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.gen_RA1v2,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_5perc_barplot_v2.png", width=20, height=25, dpi=600,create.dir = TRUE)

f.gen_RA1v3<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.gen_RA1v3,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_5perc_barplot_v3.png", width=20, height=10, dpi=600,create.dir = TRUE)

f.gen_RA2v1<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(f.gen_RA2v1,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_10perc_barplot_v1.png", width=16, height=10, dpi=600,create.dir = TRUE)

f.gen_RA2v2<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  facet_wrap(vars(Flare), scales = "free")

ggsave(f.gen_RA2v2,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_10perc_barplot_v2.png", width=16, height=25, dpi=600,create.dir = TRUE)

f.gen_RA2v3<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.gen_RA2v3,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_10perc_barplot_v3.png", width=16, height=10, dpi=600,create.dir = TRUE)

f.gen_RA3a<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.15,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 15%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(f.gen_RA3a,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_15perc_barplot_v1.png", width=16, height=10, dpi=600,create.dir = TRUE)

f.gen_RA3b<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.15,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 15%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.gen_RA3b,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_15perc_barplot_v2.png", width=16, height=10, dpi=600,create.dir = TRUE)

f.gen_RA4<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.25,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 25%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(f.gen_RA4,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_25perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

f.gen_RA5<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.35,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 35%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(f.gen_RA5,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_35perc_barplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

# plot just Massilia
massilia.relab<-f.gen.spec_RA_meta[grepl("Massilia",f.gen.spec_RA_meta$Genus_species),]
massilia.only<-ggplot(massilia.relab[!grepl("unknown",massilia.relab$Genus_species),], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Massilia Species", x="SampleID", y="Relative Abundance", fill="Taxa",subtitle="Only species in Massilia - excluding unknown species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  guides(fill=guide_legend(ncol=1)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(massilia.only,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_Massilia_Only_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

# plot just Sphingomonas
sphingomonas.relab<-f.gen.spec_RA_meta[grepl("Sphingomonas",f.gen.spec_RA_meta$Genus_species),]
sphingomonas.only<-ggplot(sphingomonas.relab[!grepl("unknown",sphingomonas.relab$Genus_species),], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Relative Abundance of Sphingomonas Species", x="SampleID", y="Relative Abundance", fill="Taxa",subtitle="Only species in Sphingomonas - excluding unknown species")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5,
        legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  guides(fill=guide_legend(ncol=2)) +
  facet_wrap(vars(Flare), scales = "free")

ggsave(sphingomonas.only,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.Spec.RA_Sphingomonas_Only_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

# prep for heatmap
max(f.gen.spec_RA_meta$Count)
mean(f.gen.spec_RA_meta$Count)
max(f.gen.spec_RA_meta$Count)/2 # what is the mid point of the RA here?

# Heatmap by SampleID

# g.h1<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.01,], aes(SampleID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.35)+
#   theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample ID", y="Fungal Genera", title="Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#
# ggsave(g.h1,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_1perc_heatmap_A.png", width=20, height=15, dpi=600,create.dir = TRUE)
#
g.h2<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.05,], aes(SampleID, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue",mid="pink",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Fungal Genera", title="Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.h2,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_5perc_heatmap_B.png", width=16, height=10, dpi=600,create.dir = TRUE)

its2.AllTax.melt[1:4,1:4]

# Taxonomic Summary by Sample ID + Collection Date

ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.01,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.gen.spec_RA_meta$SampDate_Color[order(f.gen.spec_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

tg0<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.02,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.gen.spec_RA_meta$SampDate_Color[order(f.gen.spec_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 2%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg0,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_2perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)


tg1<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.05,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.gen.spec_RA_meta$SampDate_Color[order(f.gen.spec_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg1,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_5perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

tg2<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.10,], aes(Genus_species, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.gen.spec_RA_meta$SampDate_Color[order(f.gen.spec_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 10%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg2,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_10perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

# Taxonomic summary by Sample ID + Collection Period

# tg1<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.01,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")+coord_flip()
#
# ggsave(tg1,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_1perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)
#
# tg1a<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Genus_species == "Massilia unknown",], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),
#         axis.text.x = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="", y="Relative Abundance", title="Bacterial Genus Massilia Across Samples")
#
# ggsave(tg1a,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Massilia.RA_only_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# tg1b<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.02,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 2%")+coord_flip()
#
# ggsave(tg1b,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_2perc_taxasum.png", width=15, height=23, dpi=600,create.dir = TRUE)
#
# tg1a2<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.05,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2.5),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+coord_flip()
#
# ggsave(tg1a2,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_5perc_taxasum.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# tg1b<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.1,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 10%")+coord_flip()
#
# ggsave(tg1b,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_10perc_taxasum.png", width=18, height=10, dpi=600,create.dir = TRUE)
#
# tg1c<-ggplot(f.gen.spec_RA_meta[f.gen.spec_RA_meta$Count>0.15,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=5, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 15%")+coord_flip()
#
# ggsave(tg1c,filename = "figures/RelativeAbundance/Genus_species/ABS_ITS2_Genera.RA_15perc_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)

#### Genus (by Kingdom) Relative Abundance ####
# use dcast to count up ASVs within each Genus across all of the samples
head(its2.AllTax.melt)
# its2.AllTax.melt.g<-subset(its2.AllTax.melt, its2.AllTax.melt$Genus!="Unknown") # drop unknown genera so they don't skew analyses
# "Unknown" %in% its2.AllTax.melt.g$Genus

f.k.genus_counts <- as.data.frame(dcast(its2.AllTax.melt.g, SampleID~Genus+Species+Kingdom, value.var="Count", fun.aggregate=sum)) ###
head(f.k.genus_counts) # counts by genus per sample
dim(f.k.genus_counts)
rownames(f.k.genus_counts)<-f.k.genus_counts$SampleID
f.k.genus_counts[1:4,1:4]
f.k.genus_counts<-f.k.genus_counts[,colSums(f.k.genus_counts[,-1])>0] # drop classes that are not represented
dim(f.k.genus_counts) # sanity check that we dropped taxa with no hits

f.k.genus_RelAb<-data.frame(decostand(f.k.genus_counts[,-1], method="total", MARGIN=1, na.rm=TRUE))
# relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
rowSums(f.k.genus_RelAb) # sanity check to make sure the transformation worked!

f.k.genus_RelAb$SampleID<-rownames(f.k.genus_RelAb)
head(f.k.genus_RelAb)
#write.csv(f.k.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case

# melt down relativized data to merge with metadata
f.k.genus_m<-melt(f.k.genus_RelAb)

head(f.k.genus_m)
colnames(f.k.genus_m)[which(names(f.k.genus_m) == "variable")] <- "Genus_species_Kingdom"
colnames(f.k.genus_m)[which(names(f.k.genus_m) == "value")] <- "Count"
head(f.k.genus_m) ## relative abundance based on sum of counts by genus!
f.k.genus_m$Genus_species_Kingdom<-gsub("^X.","",f.k.genus_m$Genus_species_Kingdom) # get rid of leading X. in Genus_species_Kingdom names
f.k.genus_m$Genus_species_Kingdom<-gsub("\\.\\."," ",f.k.genus_m$Genus_species_Kingdom) # get rid of .. in species name --> . is regex
f.k.genus_m$Genus_species_Kingdom<-gsub("\\."," ",f.k.genus_m$Genus_species_Kingdom) # get rid of . in species name --> . is regex

#f.k.genus_m$Genus_species_Kingdom<-gsub("_"," ",f.k.genus_m$Genus_species_Kingdom) #
head(f.k.genus_m) ## relative abundance based on sum of counts by genus!

f.k.genus_RA_meta<-merge(f.k.genus_m,metadata, by="SampleID")
head(f.k.genus_RA_meta) ## relative abundance based on sum of counts by genus!
max(f.k.genus_RA_meta$Count)
f.k.genus_RA_meta$SampleID = factor(f.k.genus_RA_meta$SampleID, levels=unique(f.k.genus_RA_meta$SampleID[order(f.k.genus_RA_meta$Flare,f.k.genus_RA_meta$Seas_Coll_Year)]), ordered=TRUE)
f.k.genus_RA_meta$Sample_Type<-"Dust"

saveRDS(f.k.genus_RA_meta, file = "data/Amplicon/SSD_GenusSpecies_RelativeAbundance_Robject.rds", ascii = FALSE, version = NULL,
        compress = TRUE, refhook = NULL)

# separate genera RelAb data by site for downstream figs
WI.k.gen.RA<-subset(f.k.genus_RA_meta,f.k.genus_RA_meta$Flare=="WI")
BDC.k.gen.RA<-subset(f.k.genus_RA_meta,f.k.genus_RA_meta$Flare=="BDC")
PD.k.gen.RA<-subset(f.k.genus_RA_meta,f.k.genus_RA_meta$Flare=="PD")
DP.k.gen.RA<-subset(f.k.genus_RA_meta,f.k.genus_RA_meta$Flare=="DP")

# Barplot by SampleID

f.k.gen_RAall<-ggplot(f.k.genus_RA_meta, aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=10))

ggsave(f.k.gen_RAall,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

f.k.gen_RA0.0<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.005,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 0.05%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=10))

ggsave(f.k.gen_RA0.0,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

f.k.gen_RA0<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=5))

ggsave(f.k.gen_RA0,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_1perc_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

f.k.gen_RA0v2<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.k.gen_RA0v2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_1perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA01<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.02,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 2%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=5))

ggsave(f.k.gen_RA01,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_2perc_barplot.png", width=30, height=10, dpi=600,create.dir = TRUE)

f.k.gen_RA01v2<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.02,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 2%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.k.gen_RA01v2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_2perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA1a<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))

ggsave(f.k.gen_RA1a,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_5perc_barplot_v1.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA1v2<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.k.gen_RA1v2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_5perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA1v3<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=2))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.k.gen_RA1v3,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_5perc_barplot_v3.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA2v1<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(f.k.gen_RA2v1,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_10perc_barplot_v1.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA2v2<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.k.gen_RA2v2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_10perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA2v3<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.10,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 10%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.k.gen_RA2v3,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_10perc_barplot_v3.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA3a<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.15,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 15%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(f.k.gen_RA3a,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_15perc_barplot_v1.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA3b<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.15,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 15%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Flare), scales = "free")

ggsave(f.k.gen_RA3b,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_15perc_barplot_v2.png", width=30, height=20, dpi=600,create.dir = TRUE)

f.k.gen_RA4<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.25,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 25%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(f.k.gen_RA4,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_25perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

f.k.gen_RA5<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.35,], aes(x=SampleID, y=Count, fill=Genus_species_Kingdom))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Includes taxa with Relative Abundance > 35%",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))

ggsave(f.k.gen_RA5,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.Spec.RA_35perc_barplot.png", width=12, height=10, dpi=600,create.dir = TRUE)

# prep for heatmap
max(f.k.genus_RA_meta$Count)
mean(f.k.genus_RA_meta$Count)
max(f.k.genus_RA_meta$Count)/2 # what is the mid point of the RA here?

# Heatmap by SampleID

# g.h1<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.01,], aes(SampleID, Genus_species_Kingdom, fill= Count)) +geom_tile()+scale_fill_gradient2(low="orange",mid="white",high="purple",midpoint=0.35)+
#   theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample ID", y="Fungal Genera", title="Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#
# ggsave(g.h1,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_1perc_heatmap_A.png", width=20, height=15, dpi=600,create.dir = TRUE)
#
g.h2<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.05,], aes(SampleID, Genus_species_Kingdom, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue",mid="pink",high="red",midpoint=0.3)+
  theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  labs(x="Sample ID", y="Fungal Genera", title="Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))

ggsave(g.h2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_5perc_heatmap_B.png", width=30, height=20, dpi=600,create.dir = TRUE)

its2.AllTax.melt[1:4,1:4]

# Taxonomic Summary by Sample ID + Collection Date

ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.01,], aes(Genus_species_Kingdom, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.k.genus_RA_meta$SampDate_Color[order(f.k.genus_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

tg0<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.02,], aes(Genus_species_Kingdom, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.k.genus_RA_meta$SampDate_Color[order(f.k.genus_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 2%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg0,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_2perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)


tg1<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.05,], aes(Genus_species_Kingdom, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.k.genus_RA_meta$SampDate_Color[order(f.k.genus_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg1,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_5perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

tg2<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.10,], aes(Genus_species_Kingdom, Count)) +
  geom_jitter(aes(color=factor(SampDate),shape=Flare), size=4, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(f.k.genus_RA_meta$SampDate_Color[order(f.k.genus_RA_meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
  labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 10%")+
  scale_shape_manual(values = c(7,10, 15,16)) + coord_flip()

ggsave(tg2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_10perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)

# Taxonomic summary by Sample ID + Collection Period

# tg1<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.01,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 1%")+coord_flip()
#
# ggsave(tg1,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_1perc_taxasum.png", width=20, height=23, dpi=600,create.dir = TRUE)
#
# tg1a<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Genus_species_Kingdom == "Massilia unknown",], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),
#         axis.text.x = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="", y="Relative Abundance", title="Bacterial Genus Massilia Across Samples")
#
# ggsave(tg1a,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Massilia.RA_only_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# tg1b<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.02,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 2%")+coord_flip()
#
# ggsave(tg1b,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_2perc_taxasum.png", width=15, height=23, dpi=600,create.dir = TRUE)
#
# tg1a2<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.05,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2.5),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%")+coord_flip()
#
# ggsave(tg1a2,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_5perc_taxasum.png", width=25, height=15, dpi=600,create.dir = TRUE)
#
# tg1b<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.1,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 10%")+coord_flip()
#
# ggsave(tg1b,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_10perc_taxasum.png", width=18, height=10, dpi=600,create.dir = TRUE)
#
# tg1c<-ggplot(f.k.genus_RA_meta[f.k.genus_RA_meta$Count>0.15,], aes(Genus_species_Kingdom, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=5, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Salton Sea Dust: Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 15%")+coord_flip()
#
# ggsave(tg1c,filename = "figures/RelativeAbundance/Genus_by_Kingdom/ABS_ITS2_Genera.RA_15perc_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)


#### Save Just Relative Abundance Data ####
save.image("data/Amplicon/SSDust_RelAb_Scaled_by_DeploymentDays_Data.Rdata")

#### Core Microbiome Analysis ####
# code for this section is from here: https://microbiome.githuf.io/tutorials/Core.html
# detection = detection threshold for absence/presence (strictly greater by default)
## aka detection threshold is based on relative abundance
# prevalence: how many samples have this taxa at the detection (relative abundance) threshold or higher?
f.gs.df <- as.data.frame(dcast(f.gen.spec_m, SampleID~Genus_species, value.var="Count", fun.aggregate=sum)) ###
f.gs.df[1:4,1:10] # using relative abundance of genus_species data
rownames(f.gs.df)<-f.gs.df$SampleID

f.gs.mat<-as.matrix(t(f.gs.df[,-1]))
f.gs.mat[1:4,1:4]

# relative population frequencies aka at 1% relative abundance threshold
head(prevalence(f.gs.mat, detection = 1/100, sort = TRUE),10)
# ^ output of prevalence(): For each OTU, the fraction of samples where a given OTU is detected. The output is readily given as a percentage.

# absolute population frequencies based on sample counts
head(prevalence(f.gs.mat, detection = 1/100, sort = TRUE, count = TRUE),10)

core.taxa.names<-core_members(f.gs.mat,detection=0.005,prevalence=40/100)
core.taxa.names # taxa that appear in half the samples

# With compositional (relative) abundances
det <- c(0.1, 0.5, 1, 2, 5, 10, 20)/100
det
prevalences <- seq(.25, 1, .05)
prevalences
#ggplot(d) + geom_point(aes(x, y)) + scale_x_continuous(trans="log10", limits=c(NA,1))


plot_core(f.gs.mat,
          prevalences = prevalences,
          detections = det,
          plot.type = "lineplot") +
  xlab("Relative Abundance (%)")

## a nicer core mircrobiome heatmap...

# create sequence of prevalences for heatmap, going from 0.05 to 1 in increments of 0.05
prevalences <- seq(.25, 1, .05)
prevalences

# create detection levels
det <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3, 3.5, 4)/100
det
#min-prevalence gets the 100th highest prevalence
core.p <- plot_core(f.gs.mat,
               plot.type = "heatmap",
               colours = brewer.pal(5, "Reds"),
               prevalences = prevalences,
               detections = det,
               min.prevalence = prevalence(f.gs.mat, sort = TRUE)[100]) +
  labs(x = "Detection Threshold\n(Relative Abundance (%))") +
  #Adjusts axis text size and legend bar height
  theme_bw() +
  theme(axis.text.y= element_text(size=14, face="italic"),
        axis.text.x.bottom=element_text(size=14),
        axis.title = element_text(size=16),
        legend.text = element_text(size=10),
        legend.title = element_text(size=14))

print(core.p)

ggsave(core.p,filename = "figures/RelativeAbundance/CoreMicrobiome/ABS_ITS2_CoreMicrobiomeGenera_heatmap.png", width=20, height=15, dpi=600,create.dir = TRUE)

# add prevalence labels to heatmap
core.p2<-core.p + geom_text(aes(label = round(core.p$data$Prevalence,2)), size = 6) +
  labs(title="Core Microbiome Heatmap",subtitle="Includes Overall Prevalence of Taxa") +
  theme(plot.title = element_text(size=15),plot.subtitle = element_text(size=13))
ggsave(core.p2,filename = "figures/RelativeAbundance/CoreMicrobiome/ABS_ITS2_CoreMicrobiomeGenera_heatmap2.png", width=22, height=15, dpi=600,create.dir = TRUE)


# pull out core microbiome taxa names
core.taxa.names<-sapply(core.p[["data"]][["Taxa"]],levels)[,1]
f.gen.spec_RA_meta[1:4,]

core.gen.meta<-f.gen.spec_RA_meta[(f.gen.spec_RA_meta$Genus_species %in% core.taxa.names),]

core.barplot1<-ggplot(core.gen.meta, aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Core Salton Sea Dust Microbiome", x="SampleID", y="Relative Abundance", subtitle="",fill="Taxa")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  facet_wrap(vars(Flare), scales = "free")

ggsave(core.barplot1,filename = "figures/RelativeAbundance/CoreMicrobiome/ABS_ITS2_CoreDustMicrobiome_byFlare_barplot.png", width=20, height=15, dpi=600,create.dir = TRUE)

core.barplot2<-ggplot(core.gen.meta, aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Core Salton Sea Dust Microbiome", x="SampleID", y="Relative Abundance", subtitle="",fill="Taxa")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  facet_wrap(vars(SampDate), scales = "free")

ggsave(core.barplot2,filename = "figures/RelativeAbundance/CoreMicrobiome/ABS_ITS2_CoreDustMicrobiome_bySampDate_barplot.png", width=30, height=35, dpi=600,create.dir = TRUE)

#### Look at Specific Shared Taxa Across Flares ####
head(f.genus_RA_meta)

gencol = melt(c(Massilia="cyan3",Sphingomonas="purple1",Planomicrobium="greenyellow", Hymenobacter="darkgoldenrod1", Planococcus="tomato1",
                Nibribacter="slateblue1", Devosia="red", `Allorhizobium Neorhizobium Pararhizobium Rhizobium`="brown", Pseudomonas="darkslategray2",
                Novosphingobium="springgreen2", Roseomonas="deeppink", Paracoccus="orchid1",Kocuria="yellow2"))
head(gencol)
gencol$Genus<-rownames(gencol)
gencol$Genus<-factor(gencol$Genus, levels=c("Massilia","Sphingomonas","Planomicrobium","Hymenobacter","Planococcus",
                                            "Nibribacter","Devosia","Allorhizobium Neorhizobium Pararhizobium Rhizobium","Pseudomonas",
                                            "Novosphingobium","Roseomonas","Paracoccus","Kocuria"))
gencol
colnames(gencol)[which(names(gencol) == "value")] <- "Genus_Color"

core.genus.meta<-merge(f.genus_RA_meta, gencol, by="Genus")
head(core.genus.meta)
length(unique(core.genus.meta$Genus)) # should be 13, 13 core genera

core.gen.barplot<-ggplot(core.genus.meta, aes(x=SampleID, y=Count, fill=Genus))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Fungal Genus Relative Abundance in Salton Sea Dust", x="SampleID", y="Relative Abundance", subtitle="Core Microbiome Taxa Only",fill="Genus")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+facet_wrap(vars(Flare), scales = "free") +
  scale_fill_manual(name ="Genus",values=unique(core.genus.meta$Genus_Color[order(core.genus.meta$Genus)]))

ggsave(core.gen.barplot,filename = "figures/RelativeAbundance/CoreMicrobiome/ABS_ITS2_CoreMicrobiomeGenera_barplot_v2.png", width=25, height=30, dpi=600,create.dir = TRUE)


#### Find Unique Genera Per Flare ####
head(metadata)
head(f.genus_m)

site_list<-data.frame(Flare=metadata$Flare,SampleID=metadata$SampleID)

f.g.RA.site<-merge(site_list, f.genus_m,by="SampleID")
head(f.g.RA.site)

# finding shared genera...
n_occur <- data.frame(table(f.g.RA.site$Genus_species)) # find frequencies of genera to see which are shared between sample types
n_occur[n_occur$Freq > 1,] # shows us which genera have a greater frequency than 2
g_shared_site<-g_site_meta[g_site_meta$Genus %in% n_occur$Var1[n_occur$Freq > 1],]
#write.csv(g_shared_site,"16S_Genera_SampleType_Shared.csv",row.names=FALSE)

g_not.shared_site<-subset(g_site_meta, !(g_site_meta$Genus %in% g_shared_site$Genus)) # subset based off of what is NOT in one dataframe from another data frame

WI.g1<-subset(f.g.RA.site,Flare=="WI")
PD.g1<-subset(f.g.RA.site,Flare=="PD")
BDC.g1<-subset(f.g.RA.site,Flare=="BDC")
DP.g1<-subset(f.g.RA.site,Flare=="DP")

# pull out taxa only in WI that's not in PD, BDC
WI.g<-subset(WI.g1, !(which(WI.g1$Genus_species %in% PD.g1$Genus_species))) # subset based off of what is NOT in one dataframe from another data frame

## Comparing Genera in Dust vs Seawater
# X Axis Breaks and Labels
lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
brks=seq(-0.05,0.05,0.01)
g_shared_site1<-subset(g_shared_site, Flare!="Soil")
g_shared_site1$Count2 <- ifelse(g_shared_site1$Flare == "Seawater", -1*g_shared_site1$Count, g_shared_site1$Count)
g_shared_site1<-g_shared_site1[order(-g_shared_site1$Count2,g_shared_site1$Genus),]
g_shared_site1$GenSamp<-interaction(g_shared_site1$Genus,g_shared_site1$Flare)
g_shared_site1$GenSamp<-factor(g_shared_site1$GenSamp, levels=g_shared_site1$GenSamp)
class(g_shared_site1$GenSamp)
g_shared_site1$Genus<-factor(g_shared_site1$Genus, levels=unique(g_shared_site1$Genus[sort(g_shared_site1$GenSamp)]))

share1<-ggplot(g_shared_site1, aes(x = Genus, y = -Count2, fill = Flare)) +
  geom_bar(data = subset(g_shared_site1[g_shared_site1$Count2<=-0.0005,], Flare == "Seawater"), stat = "identity") +
  geom_bar(data = subset(g_shared_site1[g_shared_site1$Count2>=0.0005,], Flare == "Dust"), stat = "identity") +
  coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_site1$Sample_Color[rev(order(g_shared_site1$Flare))]))+ylab("Relative Abundance")+
  theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Fungal Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")

ggsave(share1,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.PD_population.pyramid.png", width=12, height=10, dpi=600,create.dir = TRUE)

# #### Shared Genus Relative Abundance ####
# # first let's get RelAb of taxa by site
# f.g.site <- as.data.frame(dcast(its2.AllTax.melt.g, Flare+SampDate~Genus+Species, value.var="Count", fun.aggregate=sum)) ###
# f.g.site[1:7,1:7] # counts by genus per sample
# dim(f.g.site)
# rownames(f.g.site)<-interaction(f.g.site$Flare, f.g.site$SampDate,sep="-")
# f.g.site[1:7,1:7]
# f.g.site<-f.g.site[,colSums(f.g.site[,-c(1:2)])>0] # drop genera that are not present
# dim(f.g.site) # sanity check that we dropped taxa with no hits
#
# f.g.site_RA<-data.frame(decostand(f.g.site[,-c(1:2)], method="total", MARGIN=1, na.rm=TRUE))
# # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
# rowSums(f.g.site_RA) # sanity check to make sure the transformation worked!
#
# f.g.site_RA$Flare.SampDate<-rownames(f.g.site_RA)
# f.g.site_RA[1:4,1:4]
# #write.csv(f.genus_RelAb,"16S_Genera_Relative_Abundance.csv",row.names=TRUE) # good to save just in case
#
# # melt down relativized data to merge with metadata
# f.g.site_RA.m<-melt(f.g.site_RA)
#
# head(f.g.site_RA.m)
# colnames(f.g.site_RA.m)[which(names(f.g.site_RA.m) == "variable")] <- "Genus_species"
# colnames(f.g.site_RA.m)[which(names(f.g.site_RA.m) == "value")] <- "Count"
# head(f.g.site_RA.m) ## relative abundance based on sum of counts by genus!
# f.g.site_RA.m$Genus_species<-gsub("^X.","",f.g.site_RA.m$Genus_species) # get rid of leading X. in Genus_species names
# f.g.site_RA.m$Genus_species<-gsub("\\.\\."," ",f.g.site_RA.m$Genus_species) # get rid of .. in species name --> . is regex
# f.g.site_RA.m$Genus_species<-gsub("\\."," ",f.g.site_RA.m$Genus_species) # get rid of . in species name --> . is regex
# f.g.site_RA.m$Genus_species<-gsub("_"," ",f.g.site_RA.m$Genus_species) #
# head(f.g.site_RA.m) ## relative abundance based on sum of counts by genus!
# f.g.site_RA.m <- separate(data=f.g.site_RA.m, col="Flare.SampDate", sep="-", into=c("Flare","Samp.Date"),remove=FALSE) # use separate to create new columns (Flare, SampDate) from 1 column
# f.g.site_RA.m$Samp.Date<-factor(f.g.site_RA.m$Samp.Date, levels=c("July.2020","August.2020","October.2020","November.2020",
#                                                                       "July.2021","August.2021","September.2021","December.2021"))
# f.g.site_RA.m$Flare<-factor(f.g.site_RA.m$Flare, levels=c("PD","BDC","DP","WI"))
# f.g.site_RA.m$Flare.SampDate<-factor(f.g.site_RA.m$Flare.SampDate, levels=unique(f.g.site_RA.m$Flare.SampDate[order(f.g.site_RA.m$Samp.Date,f.g.site_RA.m$Flare)]))
# f.g.site_RA.m<-f.g.site_RA.m[f.g.site_RA.m$Count>0,]
#
# ggplot(f.g.site_RA.m[f.g.site_RA.m$Count>0.025,], aes(Flare.SampDate, Genus_species, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue",mid="pink",high="red",midpoint=0.3)+
#   theme_classic()+theme(axis.title.x = element_text(size=13,vjust=-0.5),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample ID", y="Fungal Genera", title="Fungal Genera & Sample Date",subtitle="Includes taxa with Relative Abundance > 5%",fill="Relative Abundance")+scale_x_discrete(expand = c(0,0))
#
# ggsave(g.h2,filename = "figures/RelativeAbundance/ABS_ITS2_Genera.RA_5perc_heatmap_B.png", width=16, height=10, dpi=600,create.dir = TRUE)
#
# # merge metadata and RA data
# site_meta<-unique(data.frame("Flare"=metadata$Flare, "Flare_Color"=metadata$Flare_Color))
# site_meta
# g_site_meta<-merge(site_meta,f.g.site_RA.m, by="Flare")
# g_site_meta<-subset(g_site_meta, Genus_species!="Unknown unknown")
# g_site_meta<-subset(g_site_meta, Count!=0)
#
# # finding shared genera...
#
# n_occur <- data.frame(table(g_site_meta$Genus_species)) # find frequencies of genera to see which are shared between site
# n_occur[n_occur$Freq > 1,] # shows us which genera have a greater frequency than 2
# g_shared_site<-g_site_meta[g_site_meta$Genus_species %in% n_occur$Var1[n_occur$Freq > 1],]
# #write.csv(g_shared_site,"16S_Genera_SampleType_Shared.csv",row.names=FALSE)
#
# g_not.shared_site<-subset(g_site_meta, !(g_site_meta$Genus_species %in% g_shared_site$Genus_species)) # subset based off of what is NOT in one dataframe from another data frame
#
# # sh.t1<-ggplot(g_shared_site, aes(Genus_species, Count)) +
# #   geom_jitter(aes(color=factor(Flare)), size=2, width=0.15, height=0) +
# #   scale_color_manual(name ="Flare", values=unique(g_shared_site$Flare_Color[order(g_shared_site$Flare)])) +
# #   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
# #   theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
# #   labs(x="Fungal Species", y="Relative Abundance", title="Fungal Species Shared Across Flare",subtitle="Only Includes Genera Shared Across All Flares")
# #
# # ggsave(sh.t1,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_site_taxasum.png", width=23, height=10, dpi=600,create.dir = TRUE)
#
# ggplot(g_shared_site[g_shared_site$Count>=0.01,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Flare)), size=2, width=0.15, height=0) +
#   scale_color_manual(name ="Flare", values=unique(g_shared_site$Flare_Color[order(g_shared_site$Flare)])) + theme_classic() +
#   theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Species", y="Relative Abundance", title="Fungal Species Shared Across Flares",subtitle="Only Includes Genera Shared Across All Flares")
#
# #ggsave(sh.t1a,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_site_1perc_taxasum.png", width=23, height=10, dpi=600,create.dir = TRUE)
#
# ggplot(g_shared_site[g_shared_site$Count>=0.01,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Flare)), size=3, width=0.15, height=0) +
#   scale_color_manual(name ="Flare", values=unique(g_shared_site$Flare_Color[order(g_shared_site$Flare)])) + theme_classic() +
#   geom_boxplot(fill=NA, outlier.color=NA) +
#   theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Species", y="Relative Abundance", title="Fungal Species Shared Across Flares")
#
# ggplot(g_shared_site[g_shared_site$Count>=0.0005,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Flare)), size=2, width=0.15, height=0) +
#   scale_color_manual(name ="Flare", values=unique(g_shared_site$Flare_Color[order(g_shared_site$Flare)])) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Species", y="Relative Abundance", title="Fungal Species Shared Across Flares")
#
# sh.t2<-ggplot(g_shared_site[g_shared_site$Count>=0.005,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Flare)), size=3, width=0.15, height=0) +
#   scale_color_manual(name ="Flare", values=unique(g_shared_site$Flare_Color[order(g_shared_site$Flare)])) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Species", y="Relative Abundance", title="Fungal Species Shared Across Flares",subtitle="Only Includes Taxa with Relative Abundance > 0.5%")
#
# ggsave(sh.t2,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_site_0.5perc_taxasum.png", width=23, height=10, dpi=600,create.dir = TRUE)
#
# sh.t3<-ggplot(g_shared_site[g_shared_site$Count>=0.01,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Flare)), size=3, width=0.15, height=0) +
#   scale_color_manual(name ="Flare", values=unique(g_shared_site$Flare_Color[order(g_shared_site$Flare)])) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13, vjust=-1),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Species", y="Relative Abundance", title="Fungal Species Shared Across Flares",subtitle="Only Includes Taxa with Relative Abundance > 1%")
#
# ggsave(sh.t3,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_site_0.25perc_taxasum.png", width=23, height=10, dpi=600,create.dir = TRUE)
#
# ## Comparing Shared Genera in WI vs PD
# # X Axis Breaks and Labels
# lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
# brks=seq(-0.05,0.05,0.01)
# g_shared_site1<-subset(g_shared_site, Flare=="WI" | Flare=="PD")
# g_shared_site1$Count2 <- ifelse(g_shared_site1$Flare == "WI", -1*g_shared_site1$Count, g_shared_site1$Count)
# g_shared_site1<-g_shared_site1[order(-g_shared_site1$Count2,g_shared_site1$Genus_species),]
# g_shared_site1$GenSamp<-interaction(g_shared_site1$Genus_species,g_shared_site1$Flare)
# g_shared_site1$GenSamp<-factor(g_shared_site1$GenSamp, levels=g_shared_site1$GenSamp)
# class(g_shared_site1$GenSamp)
# g_shared_site1$Genus_species<-factor(g_shared_site1$Genus_species, levels=unique(g_shared_site1$Genus_species[sort(g_shared_site1$GenSamp)]))
#
# share1<-ggplot(g_shared_site1, aes(x = Genus_species, y = -Count2, fill = Flare)) +
#   geom_bar(data = subset(g_shared_site1[g_shared_site1$Count2<=-0.0005,], Flare == "WI"), stat = "identity") +
#   geom_bar(data = subset(g_shared_site1[g_shared_site1$Count2>=0.0005,], Flare == "PD"), stat = "identity") +
#   coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Flare", values=unique(g_shared_site1$Flare_Color[rev(order(g_shared_site1$Flare))]))+ylab("Relative Abundance")+
#   theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Palm Desert", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")
#
# ggsave(share1,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.PD_population.pyramid.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# pp2<-ggplot(g_shared_site1[g_shared_site1$Count>=0.0005,], aes(x = reorder(Genus_species,Count), fill = Flare,y = ifelse(test = Flare == "PD",yes = Count, no = -Count))) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = abs, limits = max(g_shared_site1$Count) * c(-1,1)) +
#   coord_flip()+scale_fill_manual(name ="Flare", values=unique(g_shared_site1$Flare_Color[order(g_shared_site1$Flare)]))+ylab("Relative Abundance")+theme_classic()+
#   theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 0.05%")+xlab("Genus species")
#
# ggsave(pp2,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.PD_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# pp3<-ggplot(g_shared_site1[g_shared_site1$Count>=0.005,], aes(x = reorder(Genus_species,Count), fill = Flare,y = ifelse(test = Flare == "PD",yes = Count, no = -Count))) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = abs, limits = max(g_shared_site1$Count) * c(-1,1)) +
#   coord_flip()+scale_fill_manual(name ="Flare", values=unique(g_shared_site1$Flare_Color[order(g_shared_site1$Flare)]))+ylab("Relative Abundance")+theme_classic()+
#   theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 0.5%")+xlab("Genus species")
#
# ggsave(pp3,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.PD_0.5perc.population.pyramid.pretty.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# pp4<-ggplot(g_shared_site1[g_shared_site1$Count>=0.010,], aes(x = reorder(Genus_species,Count), fill = Flare,y = ifelse(test = Flare == "PD",yes = Count, no = -Count))) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = abs, limits = max(g_shared_site1$Count) * c(-1,1)) +
#   coord_flip()+scale_fill_manual(name ="Flare", values=unique(g_shared_site1$Flare_Color[order(g_shared_site1$Flare)]))+ylab("Relative Abundance")+theme_classic()+
#   theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 1%")+xlab("Genus species")
#
# ggsave(pp4,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.PD_1perc_population.pyramid.pretty.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# ## Lollipop chart
#
# lg1<-ggplot(g_shared_site1[g_shared_site1$Count>=0.005,], aes(x = reorder(Genus_species,Count),
#                                                               y = ifelse(test = Flare == "PD",yes = Count, no = -Count),color=g_shared_site1$Flare[g_shared_site1$Count>=0.005])) +
#   geom_point(stat='identity',size=3)  +
#   geom_segment(aes(y = 0,
#                    x = Genus_species,
#                    yend = ifelse(test = Flare == "PD",yes = Count, no = -Count),
#                    xend = Genus_species),color = "black") +
#   coord_flip()+scale_color_manual(name ="Flare", values=unique(g_shared_site1$Flare_Color[order(g_shared_site1$Flare)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 0.5%")+xlab("Genus species")
#
# ggsave(lg1,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.PD_0.5perc_lollipop_chart.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# lg2<-ggplot(g_shared_site1[g_shared_site1$Count>=0.01,], aes(x = reorder(Genus_species,Count),
#                                                               y = ifelse(test = Flare == "PD",yes = Count, no = -Count),color=g_shared_site1$Flare[g_shared_site1$Count>=0.01])) +
#   geom_point(stat='identity',size=3)  +
#   geom_segment(aes(y = 0,
#                    x = Genus_species,
#                    yend = ifelse(test = Flare == "PD",yes = Count, no = -Count),
#                    xend = Genus_species),color = "black") +
#   coord_flip()+scale_color_manual(name ="Flare", values=unique(g_shared_site1$Flare_Color[order(g_shared_site1$Flare)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Palm Desert",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 0.5%")+xlab("Genus species")
#
# ggsave(lg2,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.PD_1perc_lollipop_chart.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# ## Comparing Shared Genera in WI vs BDC
# # X Axis Breaks and Labels
# lbls = paste0(as.character(c(seq(0.05, 0, -0.01), seq(0.01, 0.05, 0.01)))) # labels
# brks=seq(-0.05,0.05,0.01)
# g_shared_site2<-subset(g_shared_site, Flare=="WI" | Flare=="BDC")
# g_shared_site2$Count2 <- ifelse(g_shared_site2$Flare == "WI", -1*g_shared_site2$Count, g_shared_site2$Count)
# g_shared_site2<-g_shared_site2[order(-g_shared_site2$Count2,g_shared_site2$Genus_species),]
# g_shared_site2$GenSamp<-interaction(g_shared_site2$Genus_species,g_shared_site2$Flare)
# g_shared_site2$GenSamp<-factor(g_shared_site2$GenSamp, levels=g_shared_site2$GenSamp)
# class(g_shared_site2$GenSamp)
# g_shared_site2$Genus_species<-factor(g_shared_site2$Genus_species, levels=unique(g_shared_site2$Genus_species[sort(g_shared_site2$GenSamp)]))
#
# share1<-ggplot(g_shared_site2, aes(x = Genus_species, y = -Count2, fill = Flare)) +
#   geom_bar(data = subset(g_shared_site2[g_shared_site2$Count2<=-0.0005,], Flare == "WI"), stat = "identity") +
#   geom_bar(data = subset(g_shared_site2[g_shared_site2$Count2>=0.0005,], Flare == "BDC"), stat = "identity") +
#   coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Flare", values=unique(g_shared_site2$Flare_Color[rev(order(g_shared_site2$Flare))]))+ylab("Relative Abundance")+
#   theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Fungal Genera by Flare", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")
#
# ggsave(share1,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.BDC_population.pyramid.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# # pp2<-ggplot(g_shared_site2[g_shared_site2$Count>=0.0005,], aes(x = reorder(Genus_species,Count), fill = Flare,y = ifelse(test = Flare == "BDC",yes = Count, no = -Count))) +
# #   geom_bar(stat = "identity") +
# #   scale_y_continuous(labels = abs, limits = max(g_shared_site2$Count) * c(-1,1)) +
# #   coord_flip()+scale_fill_manual(name ="Flare", values=unique(g_shared_site2$Flare_Color[order(g_shared_site2$Flare)]))+ylab("Relative Abundance")+theme_classic()+
# #   theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
# #   labs(title="Fungal Genera by Flare",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 0.05%")+xlab("Genus species")
# #
# # ggsave(pp2,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.BDC_population.pyramid.pretty.05percent.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# pp3a<-ggplot(g_shared_site2[g_shared_site2$Count>=0.005,], aes(x = reorder(Genus_species,Count), fill = Flare,y = ifelse(test = Flare == "BDC",yes = Count, no = -Count))) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = abs, limits = max(g_shared_site2$Count) * c(-1,1)) +
#   coord_flip()+scale_fill_manual(name ="Flare", values=unique(g_shared_site2$Flare_Color[order(g_shared_site2$Flare)]))+ylab("Relative Abundance")+theme_classic()+
#   theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Boyd Deep Canyon",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 0.5%")+xlab("Genus species")
#
# ggsave(pp3a,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.BDC_0.5perc.population.pyramid.pretty.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# pp4a<-ggplot(g_shared_site2[g_shared_site2$Count>=0.010,], aes(x = reorder(Genus_species,Count), fill = Flare,y = ifelse(test = Flare == "BDC",yes = Count, no = -Count))) +
#   geom_bar(stat = "identity") +
#   scale_y_continuous(labels = abs, limits = max(g_shared_site2$Count) * c(-1,1)) +
#   coord_flip()+scale_fill_manual(name ="Flare", values=unique(g_shared_site2$Flare_Color[order(g_shared_site2$Flare)]))+ylab("Relative Abundance")+theme_classic()+
#   theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Boyd Deep Canyon",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 1%")+xlab("Genus species")
#
# ggsave(pp4a,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.BDC_1perc_population.pyramid.pretty.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# ## Lollipop chart
#
# lg1a<-ggplot(g_shared_site2[g_shared_site2$Count>=0.005,], aes(x = reorder(Genus_species,Count),
#                                                               y = ifelse(test = Flare == "BDC",yes = Count, no = -Count),color=g_shared_site2$Flare[g_shared_site2$Count>=0.005])) +
#   geom_point(stat='identity',size=3)  +
#   geom_segment(aes(y = 0,
#                    x = Genus_species,
#                    yend = ifelse(test = Flare == "BDC",yes = Count, no = -Count),
#                    xend = Genus_species),color = "black") +
#   coord_flip()+scale_color_manual(name ="Flare", values=unique(g_shared_site2$Flare_Color[order(g_shared_site2$Flare)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Boyd Deep Canyon",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 0.5%")+xlab("Genus species")
#
# ggsave(lg1a,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.BDC_0.5perc_lollipop_chart.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# lg2a<-ggplot(g_shared_site2[g_shared_site2$Count>=0.01,], aes(x = reorder(Genus_species,Count),
#                                                              y = ifelse(test = Flare == "BDC",yes = Count, no = -Count),color=g_shared_site2$Flare[g_shared_site2$Count>=0.01])) +
#   geom_point(stat='identity',size=3)  +
#   geom_segment(aes(y = 0,
#                    x = Genus_species,
#                    yend = ifelse(test = Flare == "BDC",yes = Count, no = -Count),
#                    xend = Genus_species),color = "black") +
#   coord_flip()+scale_color_manual(name ="Flare", values=unique(g_shared_site2$Flare_Color[order(g_shared_site2$Flare)]))+ylab("Relative Abundance")+theme_classic()+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+
#   labs(title="Shared Fungal Genera - Wister vs. Boyd Deep Canyon",subtitle="Includes Shared Taxa with Relative Abundance by Flare of at least 0.5%")+xlab("Genus species")
#
# ggsave(lg2a,filename = "figures/RelativeAbundance/ABS_ITS2_shared_Genera_WI.v.BDC_1perc_lollipop_chart.png", width=13, height=10, dpi=600,create.dir = TRUE)
#
# g.typ.05p<-na.omit(subset(g_site_meta, Count>=0.005))
#
# g.t.05<-ggplot(g.typ.05p, aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Flare)), size=2, width=0.15, height=0) + geom_boxplot(fill=NA, outlier.color=NA) +scale_color_manual(name ="Flare", values=c(unique(g.typ.05p$Flare_Color[order(g.typ.05p$Flare)])),c("Seawater","Soil", "Dust")) +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Bacteria/Archaea & Flare")
#
# ggsave(g.t.05,filename = "figures/RelativeAbundance/ABS_ITS2_Gen.0.5percRA.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# g.t.05a<-ggplot(g.typ.05p, aes(Genus_species, Count)) +
#   geom_jitter(aes(color=ifelse(Count>0.01,factor(Flare),"grey")), size=2, width=0.15, height=0) + geom_boxplot(fill=NA, outlier.color=NA) +scale_color_manual(name ="Flare", values=c(unique(g.typ.05p$Flare_Color[order(g.typ.05p$Flare)]),"grey"),c("Seawater","Soil", "Dust","<1% RA")) +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(angle=40, vjust=.93, hjust=1.01),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Bacteria/Archaea & Flare")
#
# ggsave(g.t.05a,filename = "figures/RelativeAbundance/ABS_ITS2_Gen.0.5percRA.v2.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# g.typ.1p<-subset(g_site_meta, Count>=0.01)
# g.typ.1p<-subset(g.typ.1p, Genus_species!="Unknown")
#
# g2<-ggplot(g.typ.1p, aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Flare)), size=2, width=0.15, height=0) +
#   scale_color_manual(name ="Flare", values=unique(g.typ.1p$Flare_Color[order(g.typ.1p$Flare)])) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=0.5),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Fungal Genera by Flare", subtitle="Includes Taxa with a Relative Abundance of at least 1%")+coord_flip()
#
# ggsave(g2,filename = "figures/RelativeAbundance/ABS_ITS2_Genera.RA_1percent_v1.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# f.gen_RA.st<-ggplot(g.typ.1p, aes(x=Flare, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Fungal Genera by Flare", x="Flare", y="Relative Abundance", fill="Genus_species", subtitle="Includes Taxa with a Relative Abundance of at least 1%")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=0.5),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+guides(fill=guide_legend(ncol=1))
#
# ggsave(f.gen_RA.st,filename = "figures/RelativeAbundance/bacterial_genera_1percent_RA_by_SampleType.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
#
# head(g_site_meta)
# g_site_meta.no0<-subset(g_site_meta, Count!=0)
#
# tg.h1<-ggplot(g_site_meta.no0, aes(Flare, Genus, fill= Count)) +geom_tile()+scale_fill_gradient2(low="blue3", mid="white",high="red",midpoint=.025)+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.y=element_text(margin = margin(0,0)),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Sample Type", y="Fungal Genera", title="Fungal Genera & Sample Type",fill="Relative Abundance")+theme_classic()+scale_x_discrete(expand = c(0,0))
#
# ggsave(tp.h1,filename = "figures/RelativeAbundance/ABS_ITS2_Phyla.RA_heatmap.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# #### Look at Specific Shared Taxa Across Flares ####
#
# head(g_site_meta)
# length(g_site_meta$Flare)
#
# d.gen<-subset(g_site_meta, Flare=="Dust")
# sw.gen<-subset(g_site_meta, Flare=="Seawater")
#
# colnames(d.gen)[which(names(d.gen) == "Count")] <- "D.RA"
# d.gen$D.RA<-as.numeric(d.gen$D.RA)
# #d.gen<-subset(d.gen, D.RA>0.0001)
#
# colnames(sw.gen)[which(names(sw.gen) == "Count")] <- "SW.RA"
# sw.gen$SW.RA<-as.numeric(sw.gen$SW.RA)
# #sw.gen<-subset(sw.gen, SW.RA>0.0001)
#
# head(d.gen)
# head(sw.gen)
#
# tgen.comp<-merge(d.gen, sw.gen, by="Genus")
# colnames(tgen.comp)[which(names(tgen.comp) == "Flare.x")] <- "D_Samp"
# colnames(tgen.comp)[which(names(tgen.comp) == "Flare.y")] <- "SW_Samp"
# colnames(tgen.comp)[which(names(tgen.comp) == "Sample_Color.x")] <- "D_color"
# colnames(tgen.comp)[which(names(tgen.comp) == "Sample_Color.y")] <- "SW_color"
# write.csv(tgen.comp,"16S_Genera_SampleType_Shared_RA.Separated.csv",row.names=FALSE)
#
# tgen.comp$Genus[tgen.comp$D.RA==max(tgen.comp$D.RA)] # most relatively abundant genus in dust
# tgen.comp$Genus[tgen.comp$SW.RA==max(tgen.comp$SW.RA)] # most relatively abundant genus in seawater
# tgen.comp$Gen2 <- ifelse(tgen.comp$SW.RA>0.01 | tgen.comp$D.RA>0.01, tgen.comp$Genus, "Other")
# unique(tgen.comp$Gen2) # all genera names
# length(unique(tgen.comp$Gen2)) # how many colors do we need
#
# gencol = melt(c(Bacillus="darkgreen",Blastococcus="orange",Geodermatophilus="green3", Halomonas="darkslategray4", Marinobacter="blue3", Other="black", Paracoccus="firebrick1", Roseovarius="deeppink3", Salinicoccus="cornflowerblue", Sphingomonas="purple", Spiroplasma="deepskyblue1", Truepera="darkgoldenrod2"))
# head(gencol)
# gencol$Gen2<-rownames(gencol)
# gencol
#
# tgen.comp<-merge(tgen.comp, gencol, by="Gen2")
# head(tgen.comp)
# tgen.comp$Gen_Col<-as.character(tgen.comp$Gen_Col)
#
# tgen.comp$SumRA<-tgen.comp$SW.RA+tgen.comp$D.RA
# tgen.comp<-tgen.comp[order(-tgen.comp$SumRA),]
# tgen.comp$SW.RA <- -1*tgen.comp$SW.RA
# tgen.comp$Genus<-factor(tgen.comp$Genus, levels=tgen.comp$Genus)
# class(tgen.comp$Genus)
#
# test1<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Gen2))+geom_point(size=2.5) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+scale_colour_manual(values=unique(tgen.comp$Gen_Col[order(tgen.comp$Gen2)]))+
#   labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Fungal Genera & Sample Type",subtitle="Points labeled 'Other' include Genera with a Relative Abundance of < 1%")+guides(shape = guide_legend(override.aes = list(size = 3)),col=guide_legend(title="Genus"))
#
# ggsave(test1,filename = "figures/RelativeAbundance/ABS_ITS2_Gen.RA_sample.type_scatterplot.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# test2<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Genus))  +geom_point(aes(color=ifelse(SW.RA>0.01 | D.RA>0.01,Genus,"black")),size=2) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Fungal Genera & Sample Type")
#
# ggsave(test2,filename = "figures/RelativeAbundance/ABS_ITS2_Gen.RA_site_linear1_colortest.png", width=12, height=10, dpi=600,create.dir = TRUE)
#
# ggplot(tgen.comp, aes(SW.RA, D.RA, color=ifelse(SW.RA>0.01 | D.RA>0.01,Genus,"black")))  +geom_point() + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Fungal Genera & Sample Type")
#
# tgc1<-ggplot(tgen.comp, aes(SW.RA, D.RA, color=Genus)) +geom_point(aes(color=D.RA>0.01 | SW.RA>0.01))+ theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="Seawater Relative Abundance", y="Dust Relative Abundance", title="Fungal Genera & Sample Type", subtitle="Only shows shared taxa w/ RA > 0.0001")
#
# ggsave(tgc1,filename = "figures/RelativeAbundance/ABS_ITS2_Gen.RA_site_linear1.png", width=12, height=10, dpi=600,create.dir = TRUE)
# ## Date of note: 11/1/21 vvv
# # 16S_Gen.RA_site_linear1 - cutoff is 0.0001
# # 16S_Gen.RA_site_linear2 - cutoff is 0.00001
# # 16S_Gen.RA_site_linear 3 - cutoff is 0.000001
# # ** only 2 genera shared at 0.001 cutoff: Halomonas, Truepera
#
# ggplot(g_shared_site, aes(x = Genus, y = Count2, fill = Flare)) +
#   geom_bar(data = subset(g_shared_site, Flare == "Dust"), stat = "identity") +
#   geom_bar(data = subset(g_shared_site, Flare == "Seawater"), stat = "identity") +
#   coord_flip()+scale_y_continuous(labels = lbls,breaks=brks)+theme_classic()+scale_fill_manual(name ="Sample Type", values=unique(g_shared_site$Sample_Color[order(g_shared_site$Flare)]))+ylab("Relative Abundance")+theme(axis.title.x = element_text(size=14, vjust=-1),axis.title.y = element_text(size=14, vjust=1),axis.text = element_text(size=11),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=13),plot.title = element_text(size=15))+labs(title="Fungal Genera by Sample Type", subtitle="Includes Taxa with a Relative Abundance of at least 0.05%")
#
#
#
#
# #### Look at Unique Taxa In Flares ####
# head(g_not.shared_site) # see section at 1271 to see where this df came from & how it was calculated
# # ^ contains unique taxa by site
#
# WI.uniq.b1<-subset(g_not.shared_site,Flare=="WI")
# WI.g.RA.meta<-subset(f.genus_RA_meta, Flare=="WI")
# WI.uniq.f.meta<-WI.g.RA.meta[WI.g.RA.meta$Genus_species %in% WI.uniq.b1$Genus_species,]
#
# # Barplot by SampleID
#
# WI.uniq.gen1<-ggplot(WI.uniq.f.meta[WI.uniq.f.meta$Count>0.0025,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Fungal Genera Unique to Wister", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 0.25%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
#
# ggsave(WI.uniq.gen1,filename = "figures/RelativeAbundance/Wister/ABS_ITS2_WI_Genera.Spec.RA_0.25perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# WI.uniq.gen2<-ggplot(WI.uniq.f.meta[WI.uniq.f.meta$Count>0.005,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Fungal Genera Unique to Wister", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 0.5%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
#
# ggsave(WI.uniq.gen2,filename = "figures/RelativeAbundance/Wister/ABS_ITS2_WI_Genera.Spec.RA_0.5perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# WI.uniq.gen3<-ggplot(WI.uniq.f.meta[WI.uniq.f.meta$Count>0.01,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Fungal Genera Unique to Wister", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 1%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(WI.uniq.gen3,filename = "figures/RelativeAbundance/Wister/ABS_ITS2_WI_Genera.Spec.RA_1perc_barplot.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# WI.uniq.gen4<-ggplot(WI.uniq.f.meta[WI.uniq.f.meta$Count>0.05,], aes(x=SampleID, y=Count, fill=Genus_species))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
#   labs(title = "Fungal Genera Unique to Wister", x="SampleID", y="Relative Abundance", subtitle="Includes Taxa with Relative Abundance > 5%",fill="Genus")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   guides(fill=guide_legend(ncol=1))
#
# ggsave(WI.uniq.gen4,filename = "figures/RelativeAbundance/Wister/ABS_ITS2_WI_Genera.Spec.RA_5perc_barplot.png", width=16, height=10, dpi=600,create.dir = TRUE)
#
# # Taxonomic summary by Sample ID + Collection Period
#
# wi.uniq.tg1<-ggplot(WI.uniq.f.meta[WI.uniq.f.meta$Count>0,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=3, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Fungal Genera Unique to Wister by Collection",subtitle="Includes taxa with Relative Abundance > 1%")
#
# ggsave(wi.uniq.tg1,filename = "figures/RelativeAbundance/Wister/ABS_ITS2_WI_Genera.RA_taxasum_1perc.png", width=23, height=10, dpi=600,create.dir = TRUE)
#
# wi.uniq.tg1a<-ggplot(WI.uniq.f.meta[WI.uniq.f.meta$Genus_species == "Massilia unknown",], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=4, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),
#         axis.text.x = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
#   labs(x="", y="Relative Abundance", title="Bacterial Genus Massilia Across Samples")
#
# ggsave(wi.uniq.tg1a,filename = "figures/RelativeAbundance/Wister/ABS_ITS2_WI_Massilia.RA_only_taxasum.png", width=15, height=10, dpi=600,create.dir = TRUE)
#
# wi.uniq.tg1a2<-ggplot(WI.uniq.f.meta[WI.uniq.f.meta$Count>0.05,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=3, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 2.5),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Fungal Genera Unique to Wister by Collection",subtitle="Includes taxa with Relative Abundance > 5%")
#
# ggsave(wi.uniq.tg1a2,filename = "figures/RelativeAbundance/Wister/ABS_ITS2_WI_Genera.RA_taxasum_5perc.png", width=25, height=10, dpi=600,create.dir = TRUE)
#
# wi.uniq.tg1b<-ggplot(WI.uniq.f.meta[WI.uniq.f.meta$Count>0.1,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=3, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Fungal Genera Unique to Wister by Collection",subtitle="Includes taxa with Relative Abundance > 10%")
#
# ggsave(wi.uniq.tg1b,filename = "figures/RelativeAbundance/Wister/ABS_ITS2_WI_Genera.RA_taxasum_10perc.png", width=18, height=10, dpi=600,create.dir = TRUE)
#
# wi.uniq.tg1c<-ggplot(WI.uniq.f.meta[WI.uniq.f.meta$Count>0.15,], aes(Genus_species, Count)) +
#   geom_jitter(aes(color=factor(Seas_Coll_Year),shape=Flare), size=3, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Period", values=c("#14c9cb","#2962ff","#9500ff","#ff0059","#ff8c00","#0B6623","#ffd500"), labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   geom_boxplot(fill=NA, outlier.color=NA) + theme_classic() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15),plot.margin = unit(c(1, 1, 1, 3),"cm")) +
#   labs(x="Fungal Genera", y="Relative Abundance", title="Fungal Genera Unique to Wister by Collection",subtitle="Includes taxa with Relative Abundance > 15%")
#
# ggsave(wi.uniq.tg1c,filename = "figures/RelativeAbundance/Wister/ABS_ITS2_WI_Genera.RA_taxasum_15perc.png", width=15, height=10, dpi=600,create.dir = TRUE)

#### Save Everything ####
save.image("data/Amplicon/ABS_ITS2_RelAb_DataAll.Rdata")

# Version Information
sessionInfo()

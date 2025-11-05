#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaDust")
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
  #library(ggbiplot)
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
  #library(wesanderson)
  #library(nationalparkcolors)
  library(shades)
  #library(ALDEx2)
  library(rstatix)
  #library(decontam)
  #library(ggvegan)
  library(microbiome)
  library(pairwiseAdonis)
  library(corrplot)
  library(fst)
  library(plotly)
  library(htmlwidgets)
  library(rbiom)
})

## NOTE FOR THIS SCRIPT:
# In order to generate Unifrac distances, you need to have created a Multiple Sequence Alignment & built a Phylogenetic Tree
## if you do not have this yet, go back and run 5a_Prep_Data_for_MSA_PhylogeneticTree.R, then 5b_MultipleSequenceAlignment_and_PhylogeneticTree.sh

#### Load Global Env to Import Count/ASV Tables ####
load("data/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file
#load("data/SSD_16S_CLR_WeightedUnifracDist_Ready.Rdata")

# here is your guide: https://f1000research.com/articles/5-1492

bac.ASV_round.table[1:4,1:4]
bac.ASV_round.table[(nrow(bac.ASV_round.table)-4):(nrow(bac.ASV_round.table)),(ncol(bac.ASV_round.table)-4):(ncol(bac.ASV_round.table))] # last 4 rows & cols
head(bac.ASV_tax)

rownames(bac.ASV_tax)<-bac.ASV_tax$ASV_ID # set rownames as ASV IDs

#### CLR Transform All Comp Data ####
rownames(bac.ASV_round.table)
bac.ASV_round.table[1:4,1:4]

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
b.clr<-decostand(bac.ASV_round.table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]

b.clr.mat.t<-as.matrix(t(b.clr))
b.clr.mat.t[1:4,1:4]

#### Import Phylogenetic Tree File ####

phylo.tree<-read.tree("data/SSD_16SV3V4_ASVs_Updated_ssualign.bacteria.FastTree.tre")
# ^ this object is class "phylo"
phylo.tree$tip.label[length(phylo.tree$tip.label)] # shows us the last tip in the tree, could use as an outgroup
# ASV 2423 will serve as outgroup for rooting tree to get Unifrac distances
phylo.tree$tip.label[length(phylo.tree$tip.label)-1] # shows us the second to last tip in the tree, could use as an outgroup

phylo.tree$tip.label[1] # first tip in tree

# use the ASV outgroup to root your phylogenetic tree -- need this for calculating Unifrac distance
phylo.tree.rooted<-root(phylo.tree,outgroup="ASV 2423", resolve.root = TRUE)
is.rooted(phylo.tree.rooted) # should be TRUE

phylo.tree.rooted$tip.label<-gsub("ASV ","ASV_",phylo.tree.rooted$tip.label)
phylo.tree.rooted$tip.label

#cophenetic(phylo.tree)
#dist.nodes(phylo.tree)

#### Drop ASVs from Table that are NOT in Phylo Tree ####

# first from transformed ASV tables
b.clr.mat.t[1:4,1:4]
phylo.tree.rooted$tip.label

b.clr.mat.t<-b.clr.mat.t[rownames(b.clr.mat.t) %in% phylo.tree.rooted$tip.label,]

# check dimensions and if they're the same
dim(as.data.frame(phylo.tree.rooted$tip.label))
dim(b.clr.mat.t)

# then from raw ASV table which we will use for Unweighted Unifrac Distance
bac.ASV_round.table[1:4,1:4]
bac.ASV_mat.t<-as.matrix(t(bac.ASV_round.table[,-1])) # transpose ASV table and turn into a matrix
bac.ASV_mat.t<-bac.ASV_mat.t[rownames(bac.ASV_mat.t) %in% phylo.tree.rooted$tip.label,] # drop ASVs not in tree

#### Calculate Weighted Unifrac Distance ####

b.wunifrac.dist<-beta.div(b.clr.mat.t, 'unifrac', weighted = TRUE, tree = phylo.tree.rooted)
# weighted Unifrac considers the ASV abundances, non-weighted unifrac looks at presence/absence data

class(b.wunifrac.dist) #check if this is now class dist, and it is!

#### Beta Diversity - Weighted Unifrac ####
meta.all.scaled=meta.all.scaled[colnames(b.clr.mat.t),] ## reorder metadata to match order of CLR data used to generate Unifrac distances

# creating our hierarcical clustering dendrogram
b.wunifrac_clust <- hclust(b.wunifrac.dist, method="ward.D2")

# let's make it a little nicer...
b.wunifrac_dend <- as.dendrogram(b.wunifrac_clust, hang=0.02)
b.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(b.wunifrac_dend)])
labels_colors(b.wunifrac_dend) <- b.dend_cols

colorset8 # color dendrogram by collection date
#png(filename="figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_WeightedUnifracDist_Dendrogram1.png",width = 7, height = 7, units = "in",res = 800)
plot(b.wunifrac_dend, ylab="CLR Weighted Unifrac Distance", cex = 0.1,horiz=TRUE) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topleft",legend = colorset8$SampDate,cex=.8,col = colorset8$SampDate_Color,pch = 15, bty = "n")
#dev.off()

# b.wunifrac_dend1 <- as.dendrogram(b.wunifrac_clust, hang=0.06)
# b.dend_cols1 <- as.character(meta.all.scaled$SampMonth_Color[order.dendrogram(b.wunifrac_dend1)])
# labels_colors(b.wunifrac_dend1) <- b.dend_cols1
#
# colorset2 # color dendrogram by month of collection
# png(filename="figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_WeightedUnifracDist_Dendrogram2.png",width = 7, height = 7, units = "in",res = 800)
# plot(b.wunifrac_dend1, ylab="CLR Euclidean Distance", cex = 0.6) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
# legend("topright",legend = colorset2$SampleMonth,cex=.8,col = colorset2$SampMonth_Color,pch = 15, bty = "n")
# dev.off()

# PCOA w/ Weighted Unifrac distance matrix (of CLR data)
b.wunifrac.pcoa <- pcoa(b.wunifrac.dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("SSD_16S_CLR_WeightedUnifracDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.wunifrac.pcoa$values

# extract principal coordinates
b.wunifrac.pcoa.vectors<-data.frame(b.wunifrac.pcoa$vectors)
b.wunifrac.pcoa.vectors$SampleID<-rownames(b.wunifrac.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.wunifrac.pcoa.meta<-merge(b.wunifrac.pcoa.vectors, meta.all.scaled, by.x="SampleID", by.y="SampleID")
b.wunifrac.pcoa.meta$SampleMonth
b.wunifrac.pcoa.meta$SampDate

head(b.wunifrac.pcoa.meta)
rownames(b.wunifrac.pcoa.meta)<-b.wunifrac.pcoa.meta$SampleID

head(b.wunifrac.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 41.35%, PC2 = 15.96%

save.image("data/SSD_16S_CLR_WeightedUnifracDist_Ready.Rdata")

#### Visualize PCoAs - Weighted Unifrac ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints

# create PCoA ggplot fig
pcoa1<-ggplot(b.wunifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate,shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data & Weighted Unifrac Distance",color="Sample Date",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_size_manual(values = c("2020" = 7, "2021"=4),labels=c("2020","2021")) +
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.wunifrac.pcoa.meta$SampDate_Color[order(b.wunifrac.pcoa.meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_SampDate_Site_Year_PCOA1.png", width=14, height=10, dpi=600)

pcoa2<-ggplot(b.wunifrac.pcoa.meta, aes(x=Axis.2, y=Axis.3)) +geom_point(aes(color=SampDate,shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data & Weighted Unifrac Distance",color="Sample Date",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_size_manual(values = c("2020" = 7, "2021"=4),labels=c("2020","2021")) +
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.wunifrac.pcoa.meta$SampDate_Color[order(b.wunifrac.pcoa.meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  xlab("PC2 [15.96%]") + ylab("PC3 [11.63%]")

ggsave(pcoa2,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_SampDate_Site_Year_PCOA2.png", width=14, height=10, dpi=600)

# pcoa1a<-ggplot(b.wunifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Date",values=unique(b.wunifrac.pcoa.meta$Site_Color[order(b.wunifrac.pcoa.meta$Site)]),labels=c(unique(b.wunifrac.pcoa.meta$Site[order(b.wunifrac.pcoa.meta$Site)]))) +
#   xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")
#
# ggsave(pcoa1a,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_Site_Year_PCOA1.png", width=12, height=10, dpi=600)
#
# # specific season
# pcoa3<-ggplot(b.wunifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Season",values=unique(b.wunifrac.pcoa.meta$SeasonSpec_Color[order(b.wunifrac.pcoa.meta$Season_Specific)]),
#                      labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
#   xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")
#
# ggsave(pcoa3,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600)
#
# # by collection year & site
# pcoa4a<-ggplot(b.wunifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(CollectionYear),shape=Site), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Year",values=unique(b.wunifrac.pcoa.meta$Year_Color[order(b.wunifrac.pcoa.meta$CollectionYear)])) +
#   xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")
#
# ggsave(pcoa4a,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_CollectionYear_Site_PCOA1.png", width=12, height=10, dpi=600)
#
# # by collection period & site
# pcoa5<-ggplot(b.wunifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Period",values=unique(b.wunifrac.pcoa.meta$SCY_Color[order(b.wunifrac.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")
#
# ggsave(pcoa5,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_SeasCollYr_Site_PCOA1.png", width=12, height=10, dpi=600)

# 3D PCoA

pltly.all.a<-plot_ly(b.wunifrac.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~SampDate, colors = c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     symbol=~Site,symbols = c("square-open", "circle-open","circle","diamond")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 41.35%'),
                      yaxis = list(title = 'PC2 15.96%'),
                      zaxis = list(title = 'PC3 11.63%')))

saveWidget(widget = pltly.all.a, #the plotly object,
           file = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_WeightedUnifrac_SampDate_CollYr_Site_3D_PCOA1.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

plot_ly(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$Site!="WI",], x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~SampDate, colors = c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
        symbol=~Site,symbols = c("square-open", "circle-open","circle","diamond")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 41.35%'),
                      yaxis = list(title = 'PC2 15.96%'),
                      zaxis = list(title = 'PC3 11.63%')))

## visualizing only specific sites

# WI - by collection period
wi.pc1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$Site=="WI",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate), size=5)+theme_bw()+
  labs(title="PCoA: Wister Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.wunifrac.pcoa.meta$SampDate_Color[order(b.wunifrac.pcoa.meta$SampDate)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(wi.pc1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_SeasCollYr_Wister_PCOA1.png", width=12, height=10, dpi=600)

dp.pc1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$Site=="DP",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate), size=5)+theme_bw()+
  labs(title="PCoA: Dos Palmas Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.wunifrac.pcoa.meta$SampDate_Color[order(b.wunifrac.pcoa.meta$SampDate)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(dp.pc1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_SeasCollYr_DosPalmas_PCOA1.png", width=12, height=10, dpi=600)

pd.pc1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$Site=="PD",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate), size=5)+theme_bw()+
  labs(title="PCoA: Palm Desert Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.wunifrac.pcoa.meta$SampDate_Color[order(b.wunifrac.pcoa.meta$SampDate)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(pd.pc1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_SeasCollYr_PD_PCOA1.png", width=12, height=10, dpi=600)

bdc.pc1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$Site=="BDC",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate), size=5)+theme_bw()+
  labs(title="PCoA: Boyd Deep Canyon Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.wunifrac.pcoa.meta$SampDate_Color[order(b.wunifrac.pcoa.meta$SampDate)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(bdc.pc1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_SeasCollYr_BDC_PCOA1.png", width=12, height=10, dpi=600)

# visualizing pcoa by month per year
unique(b.wunifrac.pcoa.meta$SampDate)

j2020.p1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$SampDate=="July.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.wunifrac.pcoa.meta$Site_Color[order(b.wunifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(j2020.p1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_AllSites_July2020_PCOA1.png", width=12, height=10, dpi=600)

a2020.p1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$SampDate=="August.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: August 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.wunifrac.pcoa.meta$Site_Color[order(b.wunifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(a2020.p1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_AllSites_Aug2020_PCOA1.png", width=12, height=10, dpi=600)

o2020.p1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$SampDate=="October.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: October 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.wunifrac.pcoa.meta$Site_Color[order(b.wunifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(o2020.p1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_AllSites_Oct2020_PCOA1.png", width=12, height=10, dpi=600)

n2020.p1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$SampDate=="November.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: November 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.wunifrac.pcoa.meta$Site_Color[order(b.wunifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(n2020.p1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_AllSites_Nov2020_PCOA1.png", width=12, height=10, dpi=600)

ja2021.p1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$SampDate=="July.2021" | b.wunifrac.pcoa.meta$SampDate=="August.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July & August 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.wunifrac.pcoa.meta$Site_Color[order(b.wunifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(ja2021.p1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_AllSites_July_Aug_2021_PCOA1.png", width=12, height=10, dpi=600)

s2021.p1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$SampDate=="September.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: September 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.wunifrac.pcoa.meta$Site_Color[order(b.wunifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(s2021.p1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_AllSites_Sept2021_PCOA1.png", width=12, height=10, dpi=600)

d2021.p1<-ggplot(b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$SampDate=="December.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: December 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.wunifrac.pcoa.meta$Site_Color[order(b.wunifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(d2021.p1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_AllSites_Dec2021_PCOA1.png", width=12, height=10, dpi=600)

## visualizing by year
b.wunifrac.pcoa1.20<-b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$CollectionYear=="2020",]
b.wunifrac.pcoa1.21<-b.wunifrac.pcoa.meta[b.wunifrac.pcoa.meta$CollectionYear=="2021",]

# 2020
twntytwnty.pc1<-ggplot(b.wunifrac.pcoa1.20, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Microbial Composition in 2020",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.wunifrac.pcoa1.20$SCY_Color[order(b.wunifrac.pcoa1.20$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(twntytwnty.pc1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_SeasCollYr_2020_PCOA1.png", width=12, height=10, dpi=600)

# 2021
twntytwnty1.pc1<-ggplot(b.wunifrac.pcoa1.21, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Microbial Composition in 2021",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.wunifrac.pcoa1.21$SCY_Color[order(b.wunifrac.pcoa1.21$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(twntytwnty1.pc1,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_SeasCollYr_2021_PCOA1.png", width=12, height=10, dpi=600)

#### Visualize Location & Date with PC Axes ####
# this idea was suggested by Dr. Will Porter as a way to see if we can identify compositional similiarites or dissimilarites by site and sample collection date
# using PC axes as input data, where x axis is site and y axis is the sample date in the heat map

dust.time.site<-subset(meta.all.scaled, select=c(SampleID, Site, SampDate))
b.wunifrac.pcoa.dts<-merge(b.wunifrac.pcoa.vectors, dust.time.site, by.x="SampleID", by.y="SampleID")

head(b.wunifrac.pcoa.dts)

ggplot(b.wunifrac.pcoa.dts, aes(Site, SampDate)) +
  geom_tile(aes(fill = Axis.1)) +
  geom_text(aes(label = round(Axis.1, 1))) +
  scale_fill_gradient(low = "white", high = "red")

# For heatmap color gradient, PC1
max(b.wunifrac.pcoa.dts$Axis.1, na.rm=TRUE)
max(b.wunifrac.pcoa.dts$Axis.1, na.rm=TRUE)/2
min(b.wunifrac.pcoa.dts$Axis.1, na.rm=TRUE)

pcoa.axis1.hm<-ggplot(b.wunifrac.pcoa.dts, aes(Site, SampDate, Axis.1, fill=Axis.1)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="blue", high="red",labels=c("25","-36.5","-95"),breaks=c(25,-36,-95)) + labs(title="Salton Sea Dust Microbial PCoA Axis 1 by Sample Date & Site",subtitle="Using CLR-Transformed 16S Data",fill="PC1 Values") +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0)) + geom_text(aes(label = round(Axis.1, 2)),size=9)

ggsave(pcoa.axis1.hm,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_PCoA_PC1_Site_by_SampDate.png", width=18, height=13, dpi=600)

# For heatmap color gradient, PC2
max(b.wunifrac.pcoa.dts$Axis.2, na.rm=TRUE)
max(b.wunifrac.pcoa.dts$Axis.2, na.rm=TRUE)/2
min(b.wunifrac.pcoa.dts$Axis.2, na.rm=TRUE)

pcoa.axis2.hm<-ggplot(b.wunifrac.pcoa.dts, aes(Site, SampDate, Axis.2, fill=Axis.2)) +
  geom_tile(colour="white",size=0.25) +
  scale_fill_gradient(low="blue", high="red",labels=c("45","-36","-115"),breaks=c(45,-36,-115)) + labs(title="Salton Sea Dust Microbial PCoA Axis 2 by Sample Date & Site",subtitle="Using CLR-Transformed 16S Data",fill="PC2 Values") +
  geom_text(aes(label = round(Axis.2, 2)),size=9) +
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
        axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
        axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
  xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

ggsave(pcoa.axis2.hm,filename = "figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_PCoA_PC2_Site_by_SampDate.png", width=18, height=13, dpi=600)

## Loop to Generate Heat Map for Each PC Axis

pc.plot.list<-list() # create empty list for each plot to be stored in
pc.axes<-names(b.wunifrac.pcoa.dts)[grepl("Axis",names(b.wunifrac.pcoa.dts))] # pull out names of columns in df that contain "Axis" in name

# heatmap function that you will use to generate each heatmap
hm.fxn<-function(df, x_var, y_var, f_var) {
  a <- ggplot(df, aes(x = x_var, y = y_var, fill = f_var)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "blue", high = "red") +
    geom_text(aes(label = round(f_var,2)), color = "white", size = 4) +
    coord_fixed() +
    theme_classic() +
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
          axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
          axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),text=element_text(size=14)) +
    labs(x = "",
         y = "",
         fill = "R",
         title = as.character(f_var))


  return(a)
}

# loop through variable containing col names in df of interest
# create heatmap, then adjust title and legend title, add plot to a list of plots, then save plot to file
for (i in pc.axes) {
  pc.heatmap=hm.fxn(b.wunifrac.pcoa.dts, b.wunifrac.pcoa.dts$Site, b.wunifrac.pcoa.dts$SampDate, b.wunifrac.pcoa.dts[,i])
  hm_titled = pc.heatmap + ggtitle(as.character(i)) + guides(fill=guide_legend(title="PC Values"))
  pc.plot.list[[i]]=hm_titled
  ggsave(hm_titled,filename = paste("figures/BetaDiversity/WeightedUnifrac/PCoA_Axes_Heatmaps/SSD_16S_PCoA_Site_SampDate_heatmap_",i,".png",sep=""), width=18, height=13, dpi=600) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

}


# call each plot in list by index to view plots
pc.plot.list[[1]] # PC 1
pc.plot.list[[2]] # PC 2
pc.plot.list[[9]] # PC 9
pc.plot.list[[26]] # PC 27

#
# loop.hm<-function(df, x_var, y_var, f_var) {
#   a <- ggplot(df, aes(x = x_var, y = y_var, fill = f_var)) +
#     geom_tile(color = "black") +
#     scale_fill_gradient(low = "blue", high = "red") +
#     geom_text(aes(label = round(f_var,2)), color = "white", size = 4) +
#     coord_fixed() +
#     theme_minimal() +
#     labs(x = "",
#          y = "",
#          fill = "R", # Want the legend title to be each of the column names that are looped
#          title = as.character(f_var))
#
#   #ggsave(a, file = paste0("figures/BetaDiversity/WeightedUnifrac/SSD_16S.PCoA_heatmap_", f_var,".png"), device = png, width = 15, height = 15, units = "cm")
#
#   return(a)
# }
# # loop with heatmap function to create heatmap

#### Subset Count Data & Calculate Weighted Unifrac Dist by Site ####
head(meta.all.scaled)

site_list<-unique(meta.all.scaled$Site) #define an array of string values
# go through metadata & create a list of data frames
## when metadata$Variable == element in site_list (aka x in this case), subset metadata by said element into elements of a list

# here the function(x) is using site_list aka x to subset metadata, when $Variable column == site_list
# Run the function so it's stored in Global Env
site_subsets<-lapply(site_list, function(x) {subset(meta.all.scaled, Site==x)})

site_subsets # sanity check1 (should see all elements in list)
site_subsets[[1]] # sanity check2 (see 1st element in list)
#rename the list elements

# name each element in list
names(site_subsets)<-site_list # * only do this if the order of names in site_list match order of the elements in site_subsets!
site_subsets$BDC # sanity check3 - should be able to pull dataframes by names rather than index now

# example of subsetting
site_subsets[[2]][1:3]
site_subsets$BDC[1:3] # should produce same ouptut as line above

site_subsets[[2]][1:2,1:2] # another example

# ^ subsetting to [[second dataframe]], [[row #, column #]]
site_subsets[[2]][[1,2]] # [[second dataframe]], [[row 1, column 2]]

# set up the function and run this to store it in our Global environment
df_specific.subset<-function(var_vec,var_subsets){
  # var_vec = vector of variable elements from specific categorical variable;
  ## e.g. vector of names from Site categorical variable (metadata sites)
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$Site names
  for(i in seq_along(var_vec)){
    # print(var_vec[i]) -- var_vec[i] = each element in var_vec
    # print(var_subsets[[i]]) -- var_subsets[[i]] = each sub
    df<-paste(var_vec[i])
    #print(df)
    assign(df, var_subsets[[i]], envir = .GlobalEnv)
    print(paste("Dataframe", var_vec[i] ,"done"))

  }

}

# run the function
df_specific.subset(site_list, site_subsets) # used scaled metadata quantitative values

head(WI) # sanity check
WI[1:5,] # double check that our new Variable (here Site) data frames still have scaled chemical data
rownames(WI)

head(b.clr.mat.t)
# matching data with user defined function -- here is the function, must run to store function in Global env
# note: metadata must have rownames for this to work!
match_transposed_dat<-function(t.compdata, subset_metadata){
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data = pullrow<-(is.element(colnames(t.compdata), row.names(subset_metadata)))
  subset_comp_data=t.compdata[,pullrow]
  return(subset_comp_data)
}

# loop through list containing each site's metadata and use match_dat to pair with CLR data
## also calculate weighted unifrac distance per site
# NEED phylogenetic tree for this!
for (i in seq_along(site_subsets)){
  print(site_subsets[[i]]) # shows what is in each element within list
  #print(names(site_subsets[i]))
  new.t.clr.df<-match_transposed_dat(b.clr.mat.t,site_subsets[[i]]) # first match CLR data to metadata by site

  wunifrac.dist<-beta.div(new.t.clr.df, 'unifrac', weighted = TRUE, tree = phylo.tree.rooted) # then generate Weighted Unifrac Distance by site!
  assign(paste0("b.clr.t_",names(site_subsets[i])), new.t.clr.df,envir = .GlobalEnv)
  assign(paste0("wunifrac.dist_",names(site_subsets[i])), wunifrac.dist,envir = .GlobalEnv)
}
# names(site_subsets[i]) --> gives us name of each element in list

# did the function work the way we wanted it to?

b.clr.t_WI[1:4,1:4]
rownames(WI) %in% colnames(b.clr.t_WI) # output should be TRUE
wunifrac.dist_WI

save.image("data/SSD_16S_CLR_WeightedUnifracDist_bySite_Ready.Rdata")

#### Calculate Unweighted Unifrac Distance ####
bac.ASV_mat.t[1:4,1:4]

b.unifrac.dist<-beta.div(bac.ASV_mat.t, 'unifrac', weighted = FALSE, tree = phylo.tree.rooted)
# weighted Unifrac considers the ASV abundances, non-weighted unifrac looks at presence/absence data

class(b.unifrac.dist) #check if this is now class dist, and it is!

#### Beta Diversity - Unweighted Unifrac ####
meta.all.scaled=meta.all.scaled[colnames(bac.ASV_mat.t),] ## reorder metadata to match order of CLR data used to generate Unifrac distances

# creating our hierarcical clustering dendrogram
b.unifrac_clust <- hclust(b.unifrac.dist, method="ward.D2")

# let's make it a little nicer...
b.unifrac_dend <- as.dendrogram(b.unifrac_clust, hang=0.02)
b.dend_cols <- as.character(meta.all.scaled$SampDate_Color[order.dendrogram(b.unifrac_dend)])
labels_colors(b.unifrac_dend) <- b.dend_cols

colorset8 # color dendrogram by collection date
#png(filename="figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_WeightedUnifracDist_Dendrogram1.png",width = 7, height = 7, units = "in",res = 800)
plot(b.unifrac_dend, ylab="CLR Unweighted Unifrac Distance", cex = 0.1,horiz=TRUE) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topleft",legend = colorset8$SampDate,cex=.8,col = colorset8$SampDate_Color,pch = 15, bty = "n")
#dev.off()

# b.unifrac_dend1 <- as.dendrogram(b.unifrac_clust, hang=0.06)
# b.dend_cols1 <- as.character(meta.all.scaled$SampMonth_Color[order.dendrogram(b.unifrac_dend1)])
# labels_colors(b.unifrac_dend1) <- b.dend_cols1
#
# colorset2 # color dendrogram by month of collection
# png(filename="figures/BetaDiversity/WeightedUnifrac/SSD_16S_CLR_WeightedUnifracDist_Dendrogram2.png",width = 7, height = 7, units = "in",res = 800)
# plot(b.unifrac_dend1, ylab="CLR Euclidean Distance", cex = 0.6) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
# legend("topright",legend = colorset2$SampleMonth,cex=.8,col = colorset2$SampMonth_Color,pch = 15, bty = "n")
# dev.off()

# PCOA w/ Weighted Unifrac distance matrix (of CLR data)
b.unifrac.pcoa <- pcoa(b.unifrac.dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("SSD_16S_CLR_WeightedUnifracDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.unifrac.pcoa$values

# extract principal coordinates
b.unifrac.pcoa.vectors<-data.frame(b.unifrac.pcoa$vectors)
b.unifrac.pcoa.vectors$SampleID<-rownames(b.unifrac.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.unifrac.pcoa.meta<-merge(b.unifrac.pcoa.vectors, meta.all.scaled, by.x="SampleID", by.y="SampleID")
b.unifrac.pcoa.meta$SampleMonth
b.unifrac.pcoa.meta$SampDate

head(b.unifrac.pcoa.meta)
rownames(b.unifrac.pcoa.meta)<-b.unifrac.pcoa.meta$SampleID

head(b.unifrac.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 30.37%, PC2 = 6.76%

save.image("data/SSD_16S_CLR_UnweightedUnifracDist_Ready.Rdata")

#### Visualize PCoAs - Unweighted Unifrac ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints

# create PCoA ggplot fig
pcoa1<-ggplot(b.unifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate,shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data & Unweighted Unifrac Distance",color="Sample Date",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_size_manual(values = c("2020" = 7, "2021"=4),labels=c("2020","2021")) +
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.unifrac.pcoa.meta$SampDate_Color[order(b.unifrac.pcoa.meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  xlab("PC1 [30.37%]") + ylab("PC2 [6.76%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_SampDate_Site_Year_PCOA1.png", width=14, height=10, dpi=600)

pcoa2<-ggplot(b.unifrac.pcoa.meta, aes(x=Axis.2, y=Axis.3)) +geom_point(aes(color=SampDate,shape=Site,size=CollectionYear))+theme_bw()+
  labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data & Unweighted Unifrac Distance",color="Sample Date",size="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  scale_size_manual(values = c("2020" = 7, "2021"=4),labels=c("2020","2021")) +
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Date",values=unique(b.unifrac.pcoa.meta$SampDate_Color[order(b.unifrac.pcoa.meta$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  xlab("PC1 [30.37%]") + ylab("PC2 [6.76%]")

ggsave(pcoa2,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_SampDate_Site_Year_PCOA2.png", width=14, height=10, dpi=600)

# pcoa1a<-ggplot(b.unifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Date",values=unique(b.unifrac.pcoa.meta$Site_Color[order(b.unifrac.pcoa.meta$Site)]),labels=c(unique(b.unifrac.pcoa.meta$Site[order(b.unifrac.pcoa.meta$Site)]))) +
#   xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")
#
# ggsave(pcoa1a,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_Site_Year_PCOA1.png", width=12, height=10, dpi=600)
#
# # specific season
# pcoa3<-ggplot(b.unifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Season_Specific),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Season")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Season",values=unique(b.unifrac.pcoa.meta$SeasonSpec_Color[order(b.unifrac.pcoa.meta$Season_Specific)]),
#                      labels=c("Early.Summer"="Early Summer","Late.Summer"="Late Summer","Early.Fall"="Early Fall","Late.Fall"="Late Fall","Fall.Winter"="Fall-Winter")) +
#   xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")
#
# ggsave(pcoa3,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_season_specific_PCOA1.png", width=12, height=10, dpi=600)
#
# # by collection year & site
# pcoa4a<-ggplot(b.unifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(CollectionYear),shape=Site), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Year")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Year",values=unique(b.unifrac.pcoa.meta$Year_Color[order(b.unifrac.pcoa.meta$CollectionYear)])) +
#   xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")
#
# ggsave(pcoa4a,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_CollectionYear_Site_PCOA1.png", width=12, height=10, dpi=600)
#
# # by collection period & site
# pcoa5<-ggplot(b.unifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
#   labs(title="PCoA: Bacteria/Archaea in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Period",values=unique(b.unifrac.pcoa.meta$SCY_Color[order(b.unifrac.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")
#
# ggsave(pcoa5,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_SeasCollYr_Site_PCOA1.png", width=12, height=10, dpi=600)

# 3D PCoA

pltly.all.a<-plot_ly(b.unifrac.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~SampDate, colors = c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     symbol=~Site,symbols = c("square-open", "circle-open","circle","diamond")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 30.37%'),
                      yaxis = list(title = 'PC2 6.76%'),
                      zaxis = list(title = 'PC3 4.79%')))

saveWidget(widget = pltly.all.a, #the plotly object,
           file = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_UnweightedUnifrac_SampDate_CollYr_Site_3D_PCOA1.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

## visualizing only specific sites

# WI - by collection period
wi.pc1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$Site=="WI",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate), size=5)+theme_bw()+
  labs(title="PCoA: Wister Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.unifrac.pcoa.meta$SampDate_Color[order(b.unifrac.pcoa.meta$SampDate)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(wi.pc1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_SeasCollYr_Wister_PCOA1.png", width=12, height=10, dpi=600)

dp.pc1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$Site=="DP",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate), size=5)+theme_bw()+
  labs(title="PCoA: Dos Palmas Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.unifrac.pcoa.meta$SampDate_Color[order(b.unifrac.pcoa.meta$SampDate)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(dp.pc1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_SeasCollYr_DosPalmas_PCOA1.png", width=12, height=10, dpi=600)

pd.pc1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$Site=="PD",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate), size=5)+theme_bw()+
  labs(title="PCoA: Palm Desert Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.unifrac.pcoa.meta$SampDate_Color[order(b.unifrac.pcoa.meta$SampDate)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(pd.pc1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_SeasCollYr_PD_PCOA1.png", width=12, height=10, dpi=600)

bdc.pc1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$Site=="BDC",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=SampDate), size=5)+theme_bw()+
  labs(title="PCoA: Boyd Deep Canyon Microbial Composition",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.unifrac.pcoa.meta$SampDate_Color[order(b.unifrac.pcoa.meta$SampDate)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(bdc.pc1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_SeasCollYr_BDC_PCOA1.png", width=12, height=10, dpi=600)

# visualizing pcoa by month per year
unique(b.unifrac.pcoa.meta$SampDate)

j2020.p1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$SampDate=="July.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.unifrac.pcoa.meta$Site_Color[order(b.unifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(j2020.p1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_AllSites_July2020_PCOA1.png", width=12, height=10, dpi=600)

a2020.p1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$SampDate=="August.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: August 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.unifrac.pcoa.meta$Site_Color[order(b.unifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(a2020.p1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_AllSites_Aug2020_PCOA1.png", width=12, height=10, dpi=600)

o2020.p1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$SampDate=="October.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: October 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.unifrac.pcoa.meta$Site_Color[order(b.unifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(o2020.p1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_AllSites_Oct2020_PCOA1.png", width=12, height=10, dpi=600)

n2020.p1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$SampDate=="November.2020",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: November 2020 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.unifrac.pcoa.meta$Site_Color[order(b.unifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(n2020.p1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_AllSites_Nov2020_PCOA1.png", width=12, height=10, dpi=600)

ja2021.p1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$SampDate=="July.2021" | b.unifrac.pcoa.meta$SampDate=="August.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: July & August 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.unifrac.pcoa.meta$Site_Color[order(b.unifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(ja2021.p1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_AllSites_July_Aug_2021_PCOA1.png", width=12, height=10, dpi=600)

s2021.p1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$SampDate=="September.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: September 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.unifrac.pcoa.meta$Site_Color[order(b.unifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(s2021.p1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_AllSites_Sept2021_PCOA1.png", width=12, height=10, dpi=600)

d2021.p1<-ggplot(b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$SampDate=="December.2021",], aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Site)), size=5)+theme_bw()+
  labs(title="PCoA: December 2021 Microbial Composition by Site",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Site",values=unique(b.unifrac.pcoa.meta$Site_Color[order(b.unifrac.pcoa.meta$Site)])) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(d2021.p1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_AllSites_Dec2021_PCOA1.png", width=12, height=10, dpi=600)

## visualizing by year
b.unifrac.pcoa1.20<-b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$CollectionYear=="2020",]
b.unifrac.pcoa1.21<-b.unifrac.pcoa.meta[b.unifrac.pcoa.meta$CollectionYear=="2021",]

# 2020
twntytwnty.pc1<-ggplot(b.unifrac.pcoa1.20, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Microbial Composition in 2020",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.unifrac.pcoa1.20$SCY_Color[order(b.unifrac.pcoa1.20$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(twntytwnty.pc1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_SeasCollYr_2020_PCOA1.png", width=12, height=10, dpi=600)

# 2021
twntytwnty1.pc1<-ggplot(b.unifrac.pcoa1.21, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Site), size=5)+theme_bw()+
  labs(title="PCoA: Microbial Composition in 2021",subtitle="Using Centered-Log Ratio Data",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Collection Period",values=unique(b.unifrac.pcoa1.21$SCY_Color[order(b.unifrac.pcoa1.21$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
  xlab("PC1 [41.35%]") + ylab("PC2 [15.96%]")

ggsave(twntytwnty1.pc1,filename = "figures/BetaDiversity/UnweightedUnifrac/SSD_16S_CLR_SeasCollYr_2021_PCOA1.png", width=12, height=10, dpi=600)

#### Visualize Location & Date with PC Axes ####
# this idea was suggested by Dr. Will Porter as a way to see if we can identify compositional similiarites or dissimilarites by site and sample collection date
# using PC axes as input data, where x axis is site and y axis is the sample date in the heat map

dust.time.site<-subset(meta.all.scaled, select=c(SampleID, Site, SampDate))

## next for Unweighted Unifrac
pc.uw.plot.list<-list() # create empty list for each plot to be stored in
pc.axes.uwu<-names(b.unifrac.pcoa.meta)[grepl("Axis",names(b.unifrac.pcoa.meta))] # pull out names of columns in df that contain "Axis" in name

# heatmap function that you will use to generate each heatmap
hm.fxn<-function(df, x_var, y_var, f_var) {
  a <- ggplot(df, aes(x = x_var, y = y_var, fill = f_var)) +
    geom_tile(color = "black") +
    scale_fill_gradient(low = "blue", high = "red") +
    geom_text(aes(label = round(f_var,2)), color = "white", size = 4) +
    coord_fixed() +
    theme_classic() +
    theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
          axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
          axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),text=element_text(size=14)) +
    labs(x = "",
         y = "",
         fill = "R",
         title = as.character(f_var))


  return(a)
}

# loop through variable containing col names in df of interest with new function called hm.fxn()
# create heatmap, then adjust title and legend title, add plot to a list of plots, then save plot to file
for (i in pc.axes.uwu) {
  pc.heatmap=hm.fxn(b.unifrac.pcoa.meta, b.unifrac.pcoa.meta$Site, b.unifrac.pcoa.meta$SampDate, b.unifrac.pcoa.meta[,i])
  hm_titled = pc.heatmap + ggtitle(as.character(i)) + guides(fill=guide_legend(title="PC Values"))
  pc.uw.plot.list[[i]]=hm_titled
  ggsave(hm_titled,filename = paste("figures/BetaDiversity/UnweightedUnifrac/PCoA_Axes_Heatmaps/SSD_16S_PCoA_Site_SampDate_heatmap_",i,".png",sep=""), width=18, height=13, dpi=600) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).

}


# call each plot in list by index to view plots
pc.uw.plot.list[[1]] # PC 1
pc.uw.plot.list[[2]] # PC 2
pc.uw.plot.list[[9]] # PC 9
pc.uw.plot.list[[26]] # PC 27

#
# loop.hm<-function(df, x_var, y_var, f_var) {
#   a <- ggplot(df, aes(x = x_var, y = y_var, fill = f_var)) +
#     geom_tile(color = "black") +
#     scale_fill_gradient(low = "blue", high = "red") +
#     geom_text(aes(label = round(f_var,2)), color = "white", size = 4) +
#     coord_fixed() +
#     theme_minimal() +
#     labs(x = "",
#          y = "",
#          fill = "R", # Want the legend title to be each of the column names that are looped
#          title = as.character(f_var))
#
#   #ggsave(a, file = paste0("figures/BetaDiversity/UnweightedUnifrac/SSD_16S.PCoA_heatmap_", f_var,".png"), device = png, width = 15, height = 15, units = "cm")
#
#   return(a)
# }
# # loop with heatmap function to create heatmap

#### Subset Count Data & Calculate Unweighted Unifrac Dist by Site ####
head(meta.all.scaled)

site_list<-unique(meta.all.scaled$Site) #define an array of string values
# go through metadata & create a list of data frames
## when metadata$Variable == element in site_list (aka x in this case), subset metadata by said element into elements of a list

# here the function(x) is using site_list aka x to subset metadata, when $Variable column == site_list
# Run the function so it's stored in Global Env
site_subsets<-lapply(site_list, function(x) {subset(meta.all.scaled, Site==x)})

site_subsets # sanity check1 (should see all elements in list)
site_subsets[[1]] # sanity check2 (see 1st element in list)
#rename the list elements

# name each element in list
names(site_subsets)<-site_list # * only do this if the order of names in site_list match order of the elements in site_subsets!
site_subsets$BDC # sanity check3 - should be able to pull dataframes by names rather than index now

# example of subsetting
site_subsets[[2]][1:3]
site_subsets$BDC[1:3] # should produce same ouptut as line above

site_subsets[[2]][1:2,1:2] # another example

# ^ subsetting to [[second dataframe]], [[row #, column #]]
site_subsets[[2]][[1,2]] # [[second dataframe]], [[row 1, column 2]]

# set up the function and run this to store it in our Global environment
df_specific.subset<-function(var_vec,var_subsets){
  # var_vec = vector of variable elements from specific categorical variable;
  ## e.g. vector of names from Site categorical variable (metadata sites)
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$Site names
  for(i in seq_along(var_vec)){
    # print(var_vec[i]) -- var_vec[i] = each element in var_vec
    # print(var_subsets[[i]]) -- var_subsets[[i]] = each sub
    df<-paste(var_vec[i])
    #print(df)
    assign(df, var_subsets[[i]], envir = .GlobalEnv)
    print(paste("Dataframe", var_vec[i] ,"done"))

  }

}

# run the function
df_specific.subset(site_list, site_subsets) # used scaled metadata quantitative values

head(WI) # sanity check
WI[1:5,] # double check that our new Variable (here Site) data frames still have scaled chemical data
rownames(WI)

head(bac.ASV_mat.t)
# matching data with user defined function -- here is the function, must run to store function in Global env
# note: metadata must have rownames for this to work!
match_transposed_dat<-function(t.compdata, subset_metadata){
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data = pullrow<-(is.element(colnames(t.compdata), row.names(subset_metadata)))
  subset_comp_data=t.compdata[,pullrow]
  return(subset_comp_data)
}

# loop through list containing each site's metadata and use match_dat to pair with CLR data
## also calculate weighted unifrac distance per site
# NEED phylogenetic tree for this!
for (i in seq_along(site_subsets)){
  print(site_subsets[[i]]) # shows what is in each element within list
  #print(names(site_subsets[i]))
  new.t.raw.df<-match_transposed_dat(bac.ASV_mat.t,site_subsets[[i]]) # first match CLR data to metadata by site

  unifrac.dist<-beta.div(new.t.raw.df, 'unifrac', weighted = FALSE, tree = phylo.tree.rooted) # then generate uniweighted Unifrac Distance by site!
  assign(paste0("b.raw.t_",names(site_subsets[i])), new.t.raw.df,envir = .GlobalEnv)
  assign(paste0("unifrac.dist_",names(site_subsets[i])), unifrac.dist,envir = .GlobalEnv)
}
# names(site_subsets[i]) --> gives us name of each element in list

# did the function work the way we wanted it to?

b.raw.t_WI[1:4,1:4]
rownames(WI) %in% colnames(b.raw.t_WI) # output should be TRUE
unifrac.dist_WI

save.image("data/SSD_16S_CLR_UnweightedUnifracDist_bySite_Ready.Rdata")

#### Using Shapiro-Wilk test for Checking Normality of PCoA Axes ####
shapiro.test(b.wunifrac.pcoa.vectors$Axis.1) # what is the p-value?
# W = 0.81864, p-value = 0.0002317
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(b.wunifrac.pcoa.vectors$Axis.1, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(b.wunifrac.pcoa.vectors$Axis.1, pch = 1, frame = FALSE)
qqline(b.wunifrac.pcoa.vectors$Axis.1, col = "red", lwd = 2)


#### Homogeneity of Variance of Weighted Unifrac & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(meta.all.scaled)
head(b.clr)
rownames(meta.all.scaled) %in% rownames(as.matrix(b.wunifrac.dist)) #b.clr was used to make the distance matrix b.wunifrac.dist

# first by compare dispersions by sampling date
b.disper1<-betadisper(b.wunifrac.dist, meta.all.scaled$Site)
b.disper1
b.disper1$distances

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#           PD      BDC       DP    WI
# PD           0.807000 0.964000 0.151
# BDC 0.817890          0.837000 0.043
# DP  0.958211 0.835448          0.091
# WI  0.125936 0.040790 0.088662

anova(b.disper1) # p = 0.2037 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sites
# ANOVA adjusted p-value
aov.beta.p1<-anova(b.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(b.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr       upr     p adj
# BDC-PD  6.265834 -66.77865  79.31032 0.9952078
# DP-PD   1.716804 -71.32768  74.76129 0.9998993
# WI-PD  50.509875 -22.53461 123.55436 0.2515066
# DP-BDC -4.549030 -77.59352  68.49546 0.9981457
# WI-BDC 44.244041 -28.80045 117.28853 0.3602692
# WI-DP  48.793071 -24.25142 121.83756 0.2788122

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(b.wunifrac.dist ~ Site,data=meta.all.scaled,by="terms",permutations=1000)
pnova1 # p-value = 0.2677
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 0.8031968

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod1<-pairwise.adonis(b.wunifrac.dist,meta.all.scaled$Site, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1
#       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 BDC vs DP  1  19966.65 1.3797069 0.10311937   0.206      1.000
# 2 BDC vs PD  1  16647.78 1.1640216 0.08842447   0.312      1.000
# 3 BDC vs WI  1  35906.47 1.6348471 0.11990212   0.116      0.696
# 4  DP vs PD  1  12053.79 0.8391818 0.06536100   0.441      1.000
# 5  DP vs WI  1  28050.83 1.2735918 0.09594929   0.224      1.000
# 6  PD vs WI  1  13351.16 0.6108896 0.04844144   0.736      1.000

# Visualize dispersions
png('figures/BetaDiversity/WeightedUnifrac/SSD_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(b.disper1,main = "Centroids and Dispersion based on Weighted Unifrac Distance", labels=FALSE,col=colorset6$Site_Color)
dev.off()

png('figures/BetaDiversity/WeightedUnifrac/SSD_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(b.disper1,xlab="By Site", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

# below we make the same plot as above but with scatter plot instead of boxplot

b.disper1.distances<-data.frame(Disper.Dist=b.disper1$distances)
b.disper1.distances$SampleID<-rownames(b.disper1.distances)
b.dipser1.dist.meta<-merge(meta.all.scaled,b.disper1.distances,by="SampleID")

ggplot(b.dipser1.dist.meta,aes(x=Site, y=Disper.Dist, col=SampDate))+
  geom_jitter(aes(color=factor(SampDate)), size=5, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date", values=c("green2","orange","red","purple","darkgreen","darkorange3","red4","purple4"),
                     labels=c("July 2020","August 2020","October 2020","November 2020",
                              "July 2021","August 2021","September 2021","December 2021")) + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
        axis.text = element_text(size=11),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15)) +
  stat_summary(fun.data=median_hilow,color="black",size=1,fun.args=list(conf.int=0.5),show.legend=FALSE) + ylab("Distance to the Centroid")

### now compare dispersions by site + year
b.disper2<-betadisper(b.wunifrac.dist, interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear,sep="."))
b.disper2

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(b.disper2) # p = 0.4366 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p3<-anova(b.disper2)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p3,method="bonferroni",n=length(aov.beta.p3))
# 0.8732185

TukeyHSD(b.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova3<-adonis2(b.wunifrac.dist ~ Site*CollectionYear,data=meta.all.scaled,by="terms",permutations=1000)
pnova3
#                       Df SumOfSqs      R2      F  Pr(>F)


p.adjust(pnova3$`Pr(>F)`,method="bonferroni",n=length(pnova3$`Pr(>F)`)) # adjusted pval
# 1 1 1        NA        NA

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod3<-pairwise.adonis(b.wunifrac.dist,interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear,sep="."), p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod3
#                   pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig

# Visualize dispersions
png('figures/BetaDiversity/WeightedUnifrac/SSD_pcoa_betadispersion_site_by_year.png',width = 700, height = 600, res=100)
plot(b.disper2,main = "Centroids and Dispersion based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

png('figures/BetaDiversity/WeightedUnifrac/SSD_boxplot_centroid_distance_site_by_year.png',width = 900, height = 600, res=100)
boxplot(b.disper2,xlab="By Site x Collection Year", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset6$Site_Color)
dev.off()

### now compare dispersions by sample month
b.disper3<-betadisper(b.wunifrac.dist, meta.all.scaled$SampleMonth)
b.disper3

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper3, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:

anova(b.disper3) # p = 0.6944 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sample dates
# ANOVA adjusted p-value
aov.beta.p3<-anova(b.disper3)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p3,method="bonferroni",n=length(aov.beta.p3))
# 1

TukeyHSD(b.disper3) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova4<-adonis2(b.wunifrac.dist ~ SampleMonth,data=meta.all.scaled,by="terms",permutations=1000)
pnova4
#                       Df SumOfSqs      R2      F  Pr(>F)
# SampleMonth  5   106141 0.21142 1.1796 0.2577
# Residual    22   395903 0.78858
# Total       27   502045 1.00000

p.adjust(pnova4$`Pr(>F)`,method="bonferroni",n=length(pnova4$`Pr(>F)`)) # adjusted pval
# 0.7732268  NA        NA

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod4<-pairwise.adonis(b.wunifrac.dist,meta.all.scaled$SampleMonth, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod4
#                   pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1    October vs November  1 54227.924 2.8582855 0.32266802   0.070      1.000
# 2    October vs December  1 19753.181 1.2397941 0.17124714   0.358      1.000
# 3        October vs July  1 16308.984 0.8094964 0.08252171   0.443      1.000
# 4      October vs August  1 10298.519 0.5438269 0.07208899   0.848      1.000
# 5   October vs September  1 18091.578 1.1568596 0.16164347   0.334      1.000
# 6   November vs December  1 22135.771 1.3441578 0.18302409   0.124      1.000
# 7       November vs July  1 47922.099 2.3371986 0.20615310   0.044      0.660
# 8     November vs August  1 42521.720 2.1922783 0.23849129   0.055      0.825
# 9  November vs September  1 23390.459 1.4461740 0.19421706   0.191      1.000
# 10      December vs July  1 16894.623 0.9143247 0.09222259   0.452      1.000
# 11    December vs August  1 12515.731 0.7453922 0.09623686   0.555      1.000
# 12 December vs September  1  8322.350 0.6336265 0.09551736   0.540      1.000
# 13        July vs August  1  9761.077 0.4802180 0.04582137   0.831      1.000
# 14     July vs September  1  9823.604 0.5373482 0.05634147   0.757      1.000
# 15   August vs September  1 10998.877 0.6650385 0.08676257   0.690      1.000

#### Rank Distance Comparison with ANOSIM ####
# NOTES:
# The differences observed between ANOSIM and PERMANOVA demonstrate the importance of understanding the nuances of different statistical methods and their underlying assumptions.
# While PERMANOVA assesses the multivariate dispersion between groups as a whole, ANOSIM focuses on rank distances between pairs of samples.
# The pairwise approach of ANOSIM allows for a more detailed examination of specific comparisons.

# Null hypothesis: There is no difference between the means of two or more groups of (ranked) dissimilarities.
# The R-statistic in ANOSIM is a ratio between within-group and between-group dissimilarities
# An R value close to "1.0" suggests dissimilarity between groups while an R value close to "0" suggests an even distribution of high and low ranks within and between groups.
# the higher the R value, the more dissimilar your groups are in terms of microbial community composition
# The P-value is the proportion of permutations that resulted in a value of R as large or larger than that calculated using the actual grouping factor

# check if rownames of metadata are in same order as the distance matrix
rownames(meta.all.scaled) %in% rownames(as.matrix(b.wunifrac.dist))

anosim1<-anosim(x = b.wunifrac.dist, grouping = meta.all.scaled$Site, permutations = 9999)
anosim1
# ANOSIM statistic R: 0.05224; Significance: 0.1656
# ^ sites are slightly dissimilar, result is not significant...can accept null hypothesis
plot(anosim1)

anosim2<-anosim(x = b.wunifrac.dist, grouping = interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear), permutations = 9999)
# Null hypothesis: There is no difference between the means of two or more groups of (ranked) dissimilarities.
# An R value close to "1.0" suggests dissimilarity between groups while an R value close to "0" suggests an even distribution of high and low ranks within and between groups.
# The P-value is the proportion of permutations that resulted in a value of R as large or larger than that calculated using the actual grouping factor

plot(anosim2)

#### Species Contributions to Dissimilarity with SIMPER ####
# NOTES: (some from here https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/simper/)
# When there are multiple sample units, there are also multiple dissimilarities (i.e., pairwise combinations) to consider
# each species can contribute to the dissimilarity between each pair of sample units.
## If our purpose is to quantify the contribution of each species to the differences between two groups, we have to consider all of these dissimilarities.
# The average dissimilarity for each pairwise combination can be calculated directly via the vegan::meandist() function

meandist(dist = b.wunifrac.dist, grouping = meta.all.scaled$Site)
# The distances shown here are a symmetric square matrix where each value is the average distance between two sample units.
# Values on the diagonal are the average distances between observations in the same group;
# all other values are the average distances between an observation in one group and an observation in another group.

# now for SIMPER analysis...
# Similarity percentage (SIMPER) partitions the dissimilarity matrix for every pair of sample units, and then calculates the average contribution of each species to the difference between the sample units.
# These contributions are relativized so that the average contributions of all species sum to 1.
# Statistical significance of these contributions is assessed by permuting the group identities.
# Note that the title is misleading: a high value for SIMPER means that a species has a high contribution to the difference between the two groups, NOT that it has high similarity between them!

simper1<-simper(b.clr,
                meta.all.scaled$Site,
                permutations = 999
)

simper1.summary<-summary(simper1, ordered = TRUE,
                         digits = max(3,getOption("digits") - 3))

# create vector of comparison names from SIMPER
simper1.names<-c(names(simper1.summary))
for (i in seq_along(simper1.names)){
  print(simper1.summary[i])
}

for (i in seq_along(simper1.names)){
  print(as.data.frame(simper1.summary[i]))
}

# funcion to subset SIMPER results and merge it with ASV taxonomy assignments
## then reorder merged df by p value, and output that
subset.simper<-function(simper_object){
  # comp_names = list of the comparison names from SIMPER output (vegan::simper())
  ## e.g. $BDC-WI, $PD-WI, etc
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$Site names
  simper.sum<-summary(simper_object, ordered = TRUE,
                      digits = max(3,getOption("digits") - 3)) # create table of SIMPER outputs
  comp_names<-c(names(simper.sum)) # create vector of comparison names in simper object
  for(i in seq_along(comp_names)){
    # print(comp_names[i]) -- comp_names[i] = each element in comp_names
    # print(simper_object[i]) -- simper_object[i] = each comparison in simper object
    df<-as.data.frame(simper.sum[i])
    df$ASV_ID<-rownames(df)
    df2<-merge(df,bac.ASV_tax,by="ASV_ID")
    df2.order<-df2[order(df2[,8]),]
    #print(df)
    assign(paste0(comp_names[i],"_SIMPER.results"), df, envir = .GlobalEnv)
    assign(paste0(comp_names[i],"_SIMPER_taxaIDs"), df2.order, envir = .GlobalEnv)

    print(paste("Dataframe", comp_names[i] ,"done"))

  }

}

subset.simper(simper1)

# description of output of SIMPER...
# species are in descending order from highest to lowest contribution of differences between sample units
# columns are in this order: average, sd, ratio, ava/avb, cumsum, p
# ave = average contribution of species to average dissimilarity between observations from two groups
# sd = standard deviation of contribution of this species (based on its contribution to all dissimilarities betwen obsevations from two groups)
# ava/avb = average abundance of this species in each of the two groups; only included if there was a group variable used
# cusum = cumulative contribution fo this and all previous species in list; based on average but expressed as a proportion of the average dissimilarity
# p = permutation based p value aka probability of getting a larger or equal ave contribution for each species if the grouping factor was randomly permuted


simper2<-simper(b.clr,
                interaction(meta.all.scaled$Site,meta.all.scaled$CollectionYear),
                permutations = 999
)


#### Correlations with Env Vars & PC Axes ####

# first let's visualize the correlations in a couple ways...
head(b.wunifrac.pcoa.meta)
heatmap(abs(cor(b.wunifrac.pcoa.meta[,c(2:11,32:33,36:37,61:69)])),
        # Compute pearson correlation (note they are absolute values)
        col = rev(heat.colors(6)),
        Colv = NA, Rowv = NA)
legend("topleft",
       title = "Absolute Pearson R",
       legend =  round(seq(0,1, length.out = 6),1),
       y.intersp = 0.7, bty = "n",
       fill = rev(heat.colors(6)))
dev.off()

# Visualize with a corrplot
cor_mat.env1 <- cor(b.wunifrac.pcoa.meta[,c(2:11,32:33,36:37,61:69)], method='pearson')
cor_mat.env1

symnum(cor_mat.env1)

png('figures/BetaDiversity/WeightedUnifrac/SSD_WUnifrac_EnVar_Corrplot.png')
crrplt1<-corrplot.mixed(cor_mat.env1, tl.pos='lt', tl.cex=0.7, sig.level = 0.05, number.cex=0.8,
                        diag='l',cl.ratio = 0.2, tl.srt = 45)
# env variables with a correlation of <|0.7| is a good threshold for determining if predictors correlate
dev.off()

# now let's run the correlations
head(meta.all.scaled)
SurfTypFreq[,3:12]

# create dfs of only surface type freq data, only climate data, & only the pcoa axes of interest respectively
STF_Clim_Only<-meta.all.scaled[,c(5:6,9:10,35:44)]
head(STF_Clim_Only)
pcoa.axes<-b.wunifrac.pcoa.vectors[,-(ncol(b.wunifrac.pcoa.vectors))]
head(pcoa.axes)

dim(STF_Clim_Only) # confirming that both data frames have the same # of rows
dim(pcoa.axes)

rownames(STF_Clim_Only) # check rownames to see if they are in the same order in both data frames
rownames(pcoa.axes)

# reorder data frames so they are in the same order by row (SampleID)
STF_Clim_Only=STF_Clim_Only[rownames(pcoa.axes),] ## reorder metadata to match order of CLR data

rownames(STF_Clim_Only) # check rownames to see if they are in the same order in both data frames after reordering
rownames(pcoa.axes)

multi.univar.corr.fxn<-function(dep.var.df,indep.var.df){
  # create empty lists to store stuff & model number (cornum) to keep track of models each iteration of loop in fxn
  corr_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create empty list where the corr output is stored
  results_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create an empty list where the corr summaries are stored
  sig.results<-vector('list', ncol(dep.var.df) * ncol(indep.var.df))
  cornum <- 1 # counting our model numbers for indexes purposes in the loop

  # run the nested loop that generates corrs from each data frame
  ## dep.var.df[i] is dependent variable (y), indep.var.df[j] is independent variable (x) in corr
  for (i in 1:ncol(dep.var.df)){ # for each column in dep.var.df
    for (j in 1:ncol(indep.var.df)){ # for each column in indep.var.df
      # cor.test(x, y)
      corr_[[cornum]] <- cor.test(indep.var.df[,j],dep.var.df[,i], method="pearson", alternative="two.sided")
      #corr_[[cornum]]$p.value<-p.adjust(corr_[[cornum]]$p.value,method="bonferroni",n=corr_[[cornum]]$parameter)
      results_[[cornum]] <-corr_[[cornum]] # save results of corr into list called results
      names(results_)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j]) # rename list element to contain the name of the columns used in the model

      # save only significant corrs to another list called sig.results
      ## if p-value < 0.05, save to sig.results list
      ifelse((results_[[cornum]]$p.value < 0.05), sig.results[[cornum]]<-results_[[cornum]], "Not Sig")
      names(sig.results)[cornum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j])
      cornum <- cornum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)

    }
  }

  # drop all NULL elements from sig.results list so it only includes significant corrs
  sig.results[sapply(sig.results, is.null)] <- NULL

  # assign lists to global env so they are saved there are function ends
  assign("env.wunifrac.results.corrs", results_,envir = .GlobalEnv)
  assign("env.wunifrac.sig.results.corrs", sig.results,envir = .GlobalEnv)

}

multi.univar.corr.fxn(pcoa.axes,STF_Clim_Only) # test the function!

corrplot(pcoa.axes,STF_Clim_Only)

#### Linear Regression Functions ####
# glm_<- vector('list', ncol(pcoa.axes) * ncol(STF_Clim_Only)) # create empty list where the GLM output is stored
# results_<- vector('list', ncol(pcoa.axes) * ncol(STF_Clim_Only)) # create an empty list where the GLM summaries are stored
# sig.results<-vector('list', ncol(pcoa.axes) * ncol(STF_Clim_Only))
# mdlnum <- 1 # counting our model numbers for indexes purposes in the loop
#
# # use a loop to run a bunch of GLMs
# ## pcoa.axes[i] is dependent variable (y), STF_Clim_Only[j] is independent variable (x) in GLM
# for (i in 1:ncol(pcoa.axes)){ # for each column in pcoa.axes
#   for (j in 1:ncol(STF_Clim_Only)){ # for each column in STF_Clim_Only
#     glm_[[mdlnum]] <-glm(pcoa.axes[,i]~STF_Clim_Only[,j], family=gaussian) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
#     results_[[mdlnum]] <-summary(glm_[[mdlnum]]) # save results of glm into list called results
#     names(results_)[mdlnum]<-paste(names(pcoa.axes)[i],"~",names(STF_Clim_Only)[j]) # rename list element to contain the name of the columns used in the model
#     mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)
#
#   }
# }

## pcoa.axes[i] is dependent variable (y), STF_Clim_Only[j] is independent variable (x) in GLM
# for (i in 1:ncol(pcoa.axes)){ # for each column in pcoa.axes
#   for (j in 1:ncol(STF_Clim_Only)){ # for each column in STF_Clim_Only
#     glm_[[mdlnum]] <-glm(pcoa.axes[,i]~STF_Clim_Only[,j], family=gaussian) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
#     results_[[mdlnum]] <-summary(glm_[[mdlnum]]) # save results of glm into list called results
#     names(results_)[mdlnum]<-paste(names(pcoa.axes)[i],"~",names(STF_Clim_Only)[j]) # rename list element to contain the name of the columns used in the model
#
#     ifelse(coef(results_[[mdlnum]])[,4] < 0.05, sig.results[[mdlnum]]<-results_[[mdlnum]], "Not Sig")
#     names(sig.results)[mdlnum]<-paste(names(pcoa.axes)[i],"~",names(STF_Clim_Only)[j])
#     mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)
#
#   }
# }
# sig.results[sapply(sig.results, is.null)] <- NULL

# multi.univar.glm.fxn<-function(dep.var.df,indep.var.df,distro){
#   # create empty lists to store stuff & model number (mdlnum) to keep track of models each iteration of loop in fxn
#   glm_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create empty list where the GLM output is stored
#   results_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create an empty list where the GLM summaries are stored
#   sig.results<-vector('list', ncol(dep.var.df) * ncol(indep.var.df))
#   near.sig.results<-vector('list', ncol(dep.var.df) * ncol(indep.var.df))
#   mdlnum <- 1 # counting our model numbers for indexes purposes in the loop
#
#   # run the nested loop that generates GLMs from each data frame
#   ## dep.var.df[i] is dependent variable (y), indep.var.df[j] is independent variable (x) in GLM
#   for (i in 1:ncol(dep.var.df)){ # for each column in dep.var.df
#     for (j in 1:ncol(indep.var.df)){ # for each column in indep.var.df
#       glm_[[mdlnum]] <-glm(dep.var.df[,i]~indep.var.df[,j], family=distro) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
#       results_[[mdlnum]] <-summary(glm_[[mdlnum]]) # save results of glm into list called results
#       names(results_)[mdlnum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j]) # rename list element to contain the name of the columns used in the model
#
#       # save only significant GLMs to another list called sig.results
#       ## if p-value < 0.05, save to sig.results list
#       ifelse(coef(results_[[mdlnum]])[,4] < 0.05, sig.results[[mdlnum]]<-results_[[mdlnum]], "Not Sig")
#       names(sig.results)[mdlnum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j])
#
#       # save only near significant GLMs to another list called near.sig.results
#       ## if p-value < 0.05, save to sig.results list
#       ifelse((coef(results_[[mdlnum]])[,4] > 0.05 & coef(results_[[mdlnum]])[,4] < 0.08), near.sig.results[[mdlnum]]<-results_[[mdlnum]], "Not Sig")
#       names(near.sig.results)[mdlnum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j])
#
#       mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)
#
#     }
#   }
#
#   # drop all NULL elements from sig.results list so it only includes significant GLMs
#   sig.results[sapply(sig.results, is.null)] <- NULL
#   near.sig.results[sapply(near.sig.results, is.null)] <- NULL
#
#   # assign lists to global env so they are saved there are function ends
#   assign("results.glms", results_,envir = .GlobalEnv)
#   assign("sig.results.glms", sig.results,envir = .GlobalEnv)
#   assign("near.sig.results.glms", near.sig.results,envir = .GlobalEnv)
#
#
# }

multi.univar.glm.fxn(pcoa.axes,STF_Clim_Only,gaussian) # test the function!

#### Save Everything ####
save.image("data/SSeaDust_BetaDiv_Data.Rdata")

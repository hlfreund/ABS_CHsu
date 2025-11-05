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
  library(factoextra)
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

# #### PCoA on Raw Data - Sanity Check ####
#
# # check rownames of CLR transformed ASV data & metadata
# rownames(its2.ASV_table) %in% rownames(metadata)
# metadata=metadata[rownames(its2.ASV_table[,-1]),] ## reorder metadata to match order of CLR data
#
# # calculate our Euclidean distance matrix using CLR data
# # df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
# raw.f.euc_dist <- dist(its2.ASV_table[,-1], method = "euclidean")
#
# # creating our hierarcical clustering dendrogram
# raw.f.euc_clust <- hclust(raw.f.euc_dist, method="ward.D2")
#
# # let's make it a little nicer...
# raw.f.euc_dend <- as.dendrogram(raw.f.euc_clust, hang=0.2)
# raw.f.dend_cols <- as.character(metadata$Flare_Color[order.dendrogram(raw.f.euc_dend)])
# labels_colors(raw.f.euc_dend) <- raw.f.dend_cols
#
# ## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
# #(August.2021="#ef781c",December.2021="#03045e",April.2022="#059c3f")
#
# plot(raw.f.euc_dend, ylab="CLR Euclidean Distance", cex = 0.1,horiz=TRUE) + title(main = "Fungi Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
# #legend("topright",legend = c("August 2021","December 2021","April 2022"),cex=.8,col = c("#ef781c","#03045e","#059c3f"),pch = 15, bty = "n")
# # Control is dark blue ("#218380"), #Alternaria is light blue ("#73d2de")
# dev.off()
#
# # PCOA w/ Euclidean distance matrix (of CLR data)
# raw.f.pcoa <- pcoa(raw.f.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
# #save.image("data/Amplicon/ABS_ITS2_CLR_EucDist_Ready.Rdata")
#
# # The proportion of variances explained is in its element values$Relative_eig
# raw.f.pcoa$values
#
# # extract principal coordinates
# raw.f.pcoa.vectors<-data.frame(raw.f.pcoa$vectors)
# raw.f.pcoa.vectors$SampleID<-rownames(raw.f.pcoa$vectors)
#
# # merge pcoa coordinates w/ metadata
# raw.f.pcoa.meta<-merge(raw.f.pcoa.vectors, metadata, by.x="SampleID", by.y="SampleID")
# raw.f.pcoa.meta$SampleMonth
# raw.f.pcoa.meta$Flare
#
# head(raw.f.pcoa.meta)
#
# head(raw.f.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# # PC1 = 50.78%, PC2 = 18.29%
#
# # drop outliars
# #outliars<-c("BDC.D.7.27.21","PD.D.7.27.21","WI.D.9.18.21","WI.D.7.10.20")
# #raw.f.pcoa.meta2<-raw.f.pcoa.meta[!(raw.f.pcoa.meta$SampleID %in% outliars),]
#
# #### Visualize RAW PCoAs####
#
# # create PCoA ggplot fig
# ggplot(raw.f.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(SampleMonth),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Fungi in ABS Patients",subtitle="Using Aitchison Distance",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Date",values=unique(raw.f.pcoa.meta$SampMonth_Color[order(raw.f.pcoa.meta$SampleMonth)]),labels=c(unique(raw.f.pcoa.meta$SampleMonth[order(raw.f.pcoa.meta$SampleMonth)]))) +
#   xlab("PC1 [48.24%]") + ylab("PC2 [17.3%]")
#
# ggplot(raw.f.pcoa.meta2, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Flare),shape=CollectionYear), size=5)+theme_bw()+
#   labs(title="PCoA: Fungi in ABS Patients",subtitle="Using Aitchison Distance",color="Sample Date")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Date",values=unique(raw.f.pcoa.meta$Flare_Color[order(raw.f.pcoa.meta$Flare)]),labels=c(unique(raw.f.pcoa.meta$Flare[order(raw.f.pcoa.meta$Flare)]))) +
#   xlab("PC1 [48.24%]") + ylab("PC2 [5.09%]")
#
# # by collection period & site
# ggplot(raw.f.pcoa.meta2, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=factor(Seas_Coll_Year),shape=Flare), size=5)+theme_bw()+
#   labs(title="PCoA: Fungi in ABS Patients",subtitle="Using Aitchison Distance",color="Collection Period")+theme_classic()+ theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title.align=0.5, legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
#   guides(shape = guide_legend(override.aes = list(size = 5)))+
#   scale_color_manual(name ="Collection Period",values=unique(raw.f.pcoa.meta$SCY_Color[order(raw.f.pcoa.meta$Seas_Coll_Year)]),labels=c("S.1.2020"="Summer #1 2020","S.2.2020"="Summer #2 2020","S.3.2020"="Summer #3 2020","F.1.2020"="Fall #1 2020","S.1.2021"="Summer #1 2021","S.2.2021"="Summer #2 2021","F.1.2021"="Fall #1 2021")) +
#   xlab("PC1 [48.24%]") + ylab("PC2 [5.09%]")
#

#### CLR Transform All Comp Data ####
rownames(its2.ASV_table.v2)
its2.ASV_table.v2[1:4,1:4]

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below
f.clr<-decostand(its2.ASV_table.v2[,-1],method = "clr", pseudocount = 1) #CLR transformation
f.clr[1:4,1:4]

#### Beta Diversity - All Data ####

# check rownames of CLR transformed ASV data & metadata.v2
rownames(f.clr) %in% rownames(metadata.v2)
metadata.v2=metadata.v2[rownames(f.clr),] ## reorder metadata.v2 to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
f.euc_dist <- dist(f.clr, method = "euclidean")
# f.euc_dist2 <- vegdist(f.clr, method = "euclidean") # did this to compare the two matrices to see ifthey're identical
fviz_dist(f.euc_dist, gradient = list(low = "blue", mid = "white", high = "red"))

# creating our hierarcical clustering dendrogram
f.euc_clust <- hclust(f.euc_dist, method="ward.D2")

# let's make it a little nicer...
f.euc_dend <- as.dendrogram(f.euc_clust, hang=0.02)
f.dend_cols <- as.character(metadata.v2$Flare_Color[order.dendrogram(f.euc_dend)])
labels_colors(f.euc_dend) <- f.dend_cols

colorset1 # color dendrogram by collection date
png(filename="figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_EucDist_Dendrogram1.png",width = 7, height = 7, units = "in",res = 800)
plot(f.euc_dend, ylab="CLR Euclidean Distance", cex = 0.1,horiz=TRUE) + title(main = "Fungi Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
legend("topleft",legend = unique(metadata.v2$Flare),cex=.4,col = unique(metadata.v2$Flare_Color),pch = 1, bty = "n")
dev.off()

# let's try another hierarchical cluster dendrogram with eclust()
f.clr.hc<-eclust(f.clr, "hclust",nboot = 501) # use eclust() for "enhanced" hierarchical clustering; more here: https://www.datanovia.com/en/blog/cluster-analysis-in-r-simplified-and-enhanced/
#fviz_dend(f.clr.hc, rect = TRUE) # plot the dendrogam

# PCOA w/ Euclidean distance matrix (of CLR data) aka Aitchison distance
f.pcoa <- pcoa(f.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/Amplicon/ABS_ITS2_CLR_EucDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
f.pcoa$values

# extract principal coordinates
f.pcoa.vectors<-data.frame(f.pcoa$vectors)
f.pcoa.vectors$SampleID<-rownames(f.pcoa$vectors)

# merge pcoa coordinates w/ metadata.v2
f.pcoa.meta<-merge(f.pcoa.vectors, metadata.v2, by.x="SampleID", by.y="SampleID")

head(f.pcoa.meta)
rownames(f.pcoa.meta)<-f.pcoa.meta$SampleID

head(f.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 8.04%, PC2 = 5.09%

# Use princomp() to run PCA of Euclidean distances of CLR data (i.e., still Aitchison distances) so we can get to species contributions!
## create PCA of CLR-transformed counts (this will use Euclidean distance and gives the same results as the PCoA; you can compare based on the PC axes %)
f.pca<-prcomp(f.clr)
summary(f.pca) # proportion of variance gives us variance explained by each axis; should match your pcoa results if you're using Aitchison distance
f.pca$rotation # gives the contributions of each ASV to the variance/similarity of samples across the Euclidean space

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

# K-means Clustering - will show us the groups that also appear in our PCoA
## more on K-means clustering: http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization
## more on K-means clustering: https://uc-r.githuf.io/kmeans_clustering#:~:text=The%20gap%20statistic%20compares%20the,simulations%20of%20the%20sampling%20process.

## first we are going to play around with different values for k (aka # of clusters) and plot them to see what seems best
k2 <- kmeans(f.clr, centers = 2, nstart = 25)
k3 <- kmeans(f.clr, centers = 3, nstart = 25)
k4 <- kmeans(f.clr, centers = 4, nstart = 25)
k5 <- kmeans(f.clr, centers = 5, nstart = 25)

# plots to compare different values for k (aka # clusters)
p1 <- fviz_cluster(k2, geom = "point",  data = f.clr) + ggtitle("k = 2")
p2 <- fviz_cluster(k3, geom = "point",  data = f.clr) + ggtitle("k = 3")
p3 <- fviz_cluster(k4, geom = "point",  data = f.clr) + ggtitle("k = 4")
p4 <- fviz_cluster(k5, geom = "point",  data = f.clr) + ggtitle("k = 5")

# next we can calculate the K-means clustering and get back the ideal # of clusters
f.clr.km <- eclust(f.clr, "kmeans", hc_metric="euclid",nboot = 1000)
f.clr.km

# Visualize kmeans clustering with output form eclust()
# use repel = TRUE to avoid overplotting
fviz_cluster(f.clr.km, f.clr, ellipse = TRUE,
             ellipse.level = 0.95, ellipse.alpha = 0.2,ellipse.type = "convex",outlier.color = "black")

# Gap statistic plot to determine the ideal # of clusters k
# The gap statistic compares the total intracluster variation for different values of k with their expected values under null reference distribution of the data (i.e. a distribution with no obvious clustering)
fviz_gap_stat(f.clr.km$gap_stat) # shows 1 cluster
fviz_nbclust(f.clr, kmeans, method = "gap_stat") # also shows 1 clusters
fviz_nbclust(f.clr, kmeans, method = "silhouette") # shows 2 clusters

# Optimal number of clusters using gap statistics = 1
f.clr.km$nbclust

save.image("data/ITS2_DADA2/ABS_ITS2_CLR_EucDist_Ready.Rdata")

#### Visualize PCoAs - All Data ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints

# create PCoA ggplot fig
pval<-"p-value = 0.5423" # pvalue from PERMANOVA below that I already ran comparing Flare groups

pcoa1<-ggplot(f.pcoa.meta, aes(x=Axis.1, y=Axis.2,color=Flare))+ annotate("text", x = 40, y = -26, label = pval,size=3.5) +geom_point(size=3)+theme_bw()+
  labs(color="Flare Group")+theme_classic()+stat_ellipse(geom="polygon",aes(x=Axis.1, y=Axis.2,fill=Flare),linetype=1,alpha = 0.25)+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare Group",values=unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])) +
  scale_fill_manual(name ="Flare Group",values=unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])) +
  xlab("PC1 [8.04%]") + ylab("PC2 [5.09%]")
  
ggsave(pcoa1,filename = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_FlareGroup_Ellipses_Pval_PCOA1.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa2a<-ggplot(f.pcoa.meta, aes(x=Axis.1, y=Axis.2,color=Flare))+ annotate("text", x = 40, y = -26, label = pval,size=3.5) +geom_point(size=3)+theme_bw()+
  labs(color="Flare Group")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare Group",values=unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])) +
  scale_fill_manual(name ="Flare Group",values=unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])) +
  xlab("PC1 [8.04%]") + ylab("PC2 [5.09%]")

ggsave(pcoa2a,filename = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_FlareGroup_NoEllipses_Pval_PCOA1.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa2b<-ggplot(f.pcoa.meta, aes(x=Axis.1, y=Axis.2,color=Flare))+geom_point(size=3)+theme_bw()+
  labs(color="Flare Group")+theme_classic()+stat_ellipse(geom="polygon",aes(x=Axis.1, y=Axis.2,fill=Flare),linetype=1,alpha = 0.25)+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare Group",values=unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])) +
  scale_fill_manual(name ="Flare Group",values=unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])) +
  xlab("PC1 [8.04%]") + ylab("PC2 [5.09%]")

ggsave(pcoa2b,filename = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_FlareGroup_Ellipses_NoPval_PCOA1.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa2c<-ggplot(f.pcoa.meta, aes(x=Axis.1, y=Axis.2,color=Flare)) +geom_point(size=3)+theme_bw()+
  labs(color="Flare Group")+theme_classic()+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare Group",values=unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])) +
  scale_fill_manual(name ="Flare Group",values=unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])) +
  xlab("PC1 [8.04%]") + ylab("PC2 [5.09%]")

ggsave(pcoa2c,filename = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_FlareGroup_NoEllipses_NoPval_PCOA1.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa2d<-ggplot(f.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=Flare,shape=Dx),size=3)+theme_bw()+
  labs(title="PCoA: Fungal Beta Diversity in ABS Patients",subtitle="Using Aitchison Distance",color="Flare Group")+theme_classic()+ 
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare Group",values=unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])) +
  xlab("PC1 [8.04%]") + ylab("PC2 [5.09%]")

ggsave(pcoa2d,filename = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_FlareGroup_Dx_PCOA2.png", 
       width = 1800,height = 1300,units = "px",
       dpi = 200,create.dir = TRUE)


# 3D PCoA

pltly.all.a<-plot_ly(f.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~Flare, colors = c(unique(f.pcoa.meta$Flare_Color[order(f.pcoa.meta$Flare)])),
        symbol=~Flare,symbols = c("square-open", "circle-open","circle")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 8.04%'),
                      yaxis = list(title = 'PC2 5.09%'),
                      zaxis = list(title = 'PC3 4.63%')))

# before you can run save_image(), run the following lines; follow instructions: https://search.r-project.org/CRAN/refmans/plotly/html/save_image.html
#install.packages('reticulate')
# reticulate::install_miniconda()
# reticulate::conda_install('r-reticulate', 'python-kaleido')
# reticulate::conda_install('r-reticulate', 'plotly', channel = 'plotly')
reticulate::use_miniconda('r-reticulate')

save_image(pltly.all.a, "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_Flare_3D_Aitchison_PCOA1.png",width=1200,height=1000)
save_image(pltly.all.a, "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_Flare_3D_Aitchison_PCOA2.png",width=1400,height=1100)


# save 3D plot as an HTml
saveWidget(widget = pltly.all.a, #the plotly object,
           file = "figures/BetaDiversity/Aitchison/ABS_ITS2_CLR_Flare_3D_Aitchison_PCOA1.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

# ## Loop to Generate Heat Map for Each PC Axis
# 
# pc.plot.list<-list() # create empty list for each plot to be stored in
# pc.axes<-names(f.pcoa.dts)[grepl("Axis",names(f.pcoa.dts))] # pull out names of columns in df that contain "Axis" in name
# 
# # heatmap function that you will use to generate each heatmap
# hm.fxn<-function(df, x_var, y_var, f_var) {
#   a <- ggplot(df, aes(x = x_var, y = y_var, fill = f_var)) +
#     geom_tile(color = "black") +
#     scale_fill_gradient(low = "blue", high = "red") +
#     geom_text(aes(label = round(f_var,2)), color = "white", size = 4) +
#     coord_fixed() +
#     theme_classic() +
#     theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#           axis.text = element_text(size=15),axis.text.x = element_text(),legend.text = element_text(size=15),plot.title = element_text(size=22),
#           axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14),text=element_text(size=14)) +
#     labs(x = "",
#          y = "",
#          fill = "R",
#          title = as.character(f_var))
# 
# 
#   return(a)
# }
# 
# # loop through variable containing col names in df of interest
# # create heatmap, then adjust title and legend title, add plot to a list of plots, then save plot to file
# for (i in pc.axes) {
#   pc.heatmap=hm.fxn(f.pcoa.dts, f.pcoa.dts$Flare, f.pcoa.dts$Flare, f.pcoa.dts[,i])
#   hm_titled = pc.heatmap + ggtitle(as.character(i)) + guides(fill=guide_legend(title="PC Values"))
#   pc.plot.list[[i]]=hm_titled
#   ggsave(hm_titled,filename = paste("figures/BetaDiversity/PCoA_Axes_Heatmaps/ABS_ITS2_PCoA_Flare_SampDate_heatmap_",i,".png",sep=""), width=18, height=13, dpi=600,create.dir = TRUE) # i will be the unique identifier. and change the extension from .png to whatever you like (eps, pdf etc).
# 
# }
# 
# 
# # call each plot in list by index to view plots
# pc.plot.list[[1]] # PC 1
# pc.plot.list[[2]] # PC 2
# pc.plot.list[[9]] # PC 9
# pc.plot.list[[27]] # PC 27

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
#   #ggsave(a, file = paste0("figures/BetaDiversity/Aitchison/ABS_ITS2.PCoA_heatmap_", f_var,".png"), device = png, width = 15, height = 15, units = "cm")
#
#   return(a)
# }
# # loop with heatmap function to create heatmap
# 
# 
# s.t.pcoa<-melt(f.pcoa.dts[,-1], by=c("Flare","SampDate"))
# colnames(s.t.pcoa)[which(names(s.t.pcoa) == "variable")] <- "PCoA_Axis"
# colnames(s.t.pcoa)[which(names(s.t.pcoa) == "value")] <- "PC_Axis_Value"
# 
# ggplot(s.t.pcoa, aes(Flare, SampDate, PCoA_Axis, fill=PC_Axis_Value)) +
#   geom_tile(colour="white",size=0.25) +
#   scale_fill_gradient(low="blue", high="red",labels=c(),breaks=c()) + labs(title="PCoA by Sample Date & Flare",subtitle="Using CLR-Transformed 16S Data",fill="PC Axis Value") +
#   theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title.align=0.5, legend.title = element_text(size=18),
#         axis.text = element_text(size=15),axis.text.x = element_text(angle=45, hjust=1),legend.text = element_text(size=15),plot.title = element_text(size=22),
#         axis.ticks=element_line(size=0.4),panel.grid = element_blank(),plot.subtitle=element_text(size=14)) +
#   xlab("") + ylab("") + scale_y_discrete(expand=c(0, 0))+scale_x_discrete(expand=c(0, 0))

#ggsave(sulf.hm1a,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600,create.dir = TRUE)

#### Homogeneity of Variance & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(metadata.v2)
head(f.clr)
rownames(metadata.v2) %in% rownames(f.clr) #f.clr was used to make the distance matrix f.euc_dist

# first by compare dispersions by sampling date
f.disper1<-betadisper((vegdist(f.clr,method="euclidean")), metadata.v2$Flare)
f.disper1

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(f.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#               HHP Remission Flare
# HHP                 0.52800 0.545
# Remission 0.53009           0.218
# Flare     0.55469   0.21674      

anova(f.disper1) # p = 0.4729 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across flare groups
# ANOVA adjusted p-value
aov.beta.p1<-anova(f.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(f.disper1) # tells us which category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj
# Remission-HHP    2.283555  -6.44518 11.012290 0.8048257
# Flare-HHP       -2.369814 -11.65916  6.919528 0.8134133
# Flare-Remission -4.653369 -13.75439  4.447656 0.4408499

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(f.clr ~ Flare,data=metadata.v2,method = "euclidean",by="terms",permutations= 10000)
pnova1 # p-value = 0.5423
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 1

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

f.clr.dist = (vegdist(f.clr, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(f.clr.dist,metadata.v2$Flare, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1
#                 pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1       Flare vs HHP  1  1365.937 1.0812917 0.02839430  0.2084     0.6252    
# 2 Flare vs Remission  1  1056.664 0.7859135 0.01975356  0.9923     1.0000    
# 3   HHP vs Remission  1  1523.782 1.0814167 0.02510170  0.1855     0.5565 

# Visualize dispersions
png('figures/BetaDiversity/Aitchison/SSD_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(f.disper1,main = "Centroids and Dispersion based on Aitchison Distance", col=colorset1$Flare_Color)
dev.off()

png('figures/BetaDiversity/Aitchison/SSD_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(f.disper1,xlab="By Flare", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance", col=colorset1$Flare_Color)
dev.off()


## NOW just Remission vs Flare: run the same comparison but remove HHP samples so it's just Remission vs Flare...

# first create metadata.v2 that does not include HHP samples
no.hhp.metadata.v2<-metadata.v2[!grepl("HHP", metadata.v2$Flare),]
rownames(f.clr) %in% rownames(no.hhp.metadata.v2) #f.clr was used to make the distance matrix f.euc_dist

# use no-HHP-metadata.v2 to create new f.clr df that does not have those HHP samples
f.clr.no.hhp=f.clr[rownames(no.hhp.metadata.v2),] ## reorder metadata.v2 to match order of CLR data
# ^^ this will drop the HHPs from f.clr in new df we make

# first by compare dispersions by sampling date
f.disper2<-betadisper((vegdist(f.clr.no.hhp,method="euclidean")), no.hhp.metadata.v2$Flare)
f.disper2

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(f.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#             Remission Flare
# Remission           0.233
# Flare       0.21674      

anova(f.disper2) # p = 0.2167 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across flare groups
# ANOVA adjusted p-value
aov.beta.p2<-anova(f.disper2)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p2,method="bonferroni",n=length(aov.beta.p2))

TukeyHSD(f.disper2) # tells us which category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj
#   Flare-Remission -4.653369 -12.14977 2.84303 0.2167391

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova2<-adonis2(f.clr.no.hhp ~ Flare,data=no.hhp.metadata.v2,method = "euclidean",by="terms",permutations= 10000)
pnova2 # p-value = 0.9917

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

f.clr.no.hhp.dist = (vegdist(f.clr.no.hhp, "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod2<-pairwise.adonis(f.clr.no.hhp.dist,no.hhp.metadata.v2$Flare, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod2
#                 pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1 Flare vs Remission  1  1056.664 0.7859135 0.01975356  0.9898     0.9898 

# Visualize dispersions
png('figures/BetaDiversity/Aitchison/SSD_pcoa_betadispersion_site.png',width = 700, height = 600, res=100)
plot(f.disper2,main = "Centroids and Dispersion based on Aitchison Distance")
dev.off()

png('figures/BetaDiversity/Aitchison/SSD_boxplot_centroid_distance_site.png',width = 700, height = 600, res=100)
boxplot(f.disper2,xlab="By Flare", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance")
dev.off()

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

# check if rownames of metadata.v2 are in same order as the distance matrix
rownames(meta.all.scaled) %in% rownames(as.matrix(f.euc_dist))

anosim1<-anosim(x = f.euc_dist, grouping = meta.all.scaled$Flare, permutations = 9999)
anosim1
# ANOSIM statistic R: 0.0847; Significance: 0.03741
# ^ sites are slightly dissimilar, result is significant...cannot accept null hypothesis
plot(anosim1)

anosim2<-anosim(x = f.euc_dist, grouping = interaction(meta.all.scaled$Flare,meta.all.scaled$CollectionYear), permutations = 9999)
# Null hypothesis: There is no difference between the means of two or more groups of (ranked) dissimilarities.
# An R value close to "1.0" suggests dissimilarity between groups while an R value close to "0" suggests an even distribution of high and low ranks within and between groups.
# The P-value is the proportion of permutations that resulted in a value of R as large or larger than that calculated using the actual grouping factor

plot(anosim2)

#### Species Contributions to Dissimilarity with SIMPER ####
# NOTES: (some from here https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/simper/)
# When there are multiple sample units, there are also multiple dissimilarities (i.e., pairwise combinations) to consider – each species can contribute to the dissimilarity between each pair of sample units.
## If our purpose is to quantify the contribution of each species to the differences between two groups, we have to consider all of these dissimilarities.
# The average dissimilarity for each pairwise combination can be calculated directly via the vegan::meandist() function

meandist(dist = f.euc_dist, grouping = meta.all.scaled$Flare)
# The distances shown here are a symmetric square matrix where each value is the average distance between two sample units.
# Values on the diagonal are the average distances between observations in the same group;
# all other values are the average distances between an observation in one group and an observation in another group.

# now for SIMPER analysis...
# Similarity percentage (SIMPER) partitions the dissimilarity matrix for every pair of sample units, and then calculates the average contribution of each species to the difference between the sample units.
# These contributions are relativized so that the average contributions of all species sum to 1.
# Statistical significance of these contributions is assessed by permuting the group identities.
# Note that the title is misleading: a high value for SIMPER means that a species has a high contribution to the difference between the two groups, NOT that it has high similarity between them!

simper1<-simper(f.clr,
       meta.all.scaled$Flare,
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

subset.simper<-function(simper_object){
  # comp_names = list of the comparison names from SIMPER output (vegan::simper())
  ## e.g. $BDC-WI, $PD-WI, etc
  # var_subsets = list of dataframes subsetted by column$element from original dataframe;
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata.v2 using vector of metadata.v2$Flare names
  simper.sum<-summary(simper_object, ordered = TRUE,
                           digits = max(3,getOption("digits") - 3)) # create table of SIMPER outputs
  comp_names<-c(names(simper.sum)) # create vector of comparison names in simper object
  for(i in seq_along(comp_names)){
    # print(comp_names[i]) -- comp_names[i] = each element in comp_names
    # print(simper_object[i]) -- simper_object[i] = each comparison in simper object
    df<-as.data.frame(simper.sum[i])
    df$ASV_ID<-rownames(df)
    df2<-merge(df,bac.ASV_tax,by="ASV_ID")
    #df2.order<-df2[order(df2[,grepl(".p",names(df2))],decreasing=FALSE)]
    #print(df)
    assign(paste0(comp_names[i],"_SIMPER.results"), df, envir = .GlobalEnv)
    assign(paste0(comp_names[i],"_SIMPER_taxaIDs"), df2, envir = .GlobalEnv)

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

head(simper1$BDC_DP$cusum)
head(sort(simper1$BDC_DP$p,decreasing=FALSE))

head(simper1$BDC_PD$cusum)
head(simper1$BDC_WI$cusum)

head(simper1$DP_PD$cusum)
head(simper1$DP_WI$cusum)

head(simper1$PD_WI$cusum)


#### Using Shapiro-Wilk test for Checking Normality of PCoA Axes ####
shapiro.test(f.pcoa.vectors$Axis.1) # what is the p-value?
# p-value = 0.5816
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(f.pcoa.vectors$Axis.1, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.githuf.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(f.pcoa.vectors$Axis.1, pch = 1, frame = FALSE)
qqline(f.pcoa.vectors$Axis.1, col = "red", lwd = 2)

#### Linear Models with Surface Type Frequencies & PC Axes ####
head(f.pcoa.meta)

SurfTypFreq[,3:12]

# create dfs of only surface type freq data and only the pcoa axes of interest
STFs_only<-SurfTypFreq[,3:12]
head(STFs_only)
pcoa.axes<-f.pcoa.vectors[,-(ncol(f.pcoa.vectors))]
head(pcoa.axes)

dim(STFs_only) # confirming that both data frames have the same # of rows
dim(pcoa.axes)

rownames(STFs_only) # check rownames to see if they are in the same order in both data frames
rownames(pcoa.axes)

# reorder data frames so they are in the same order by row (SampleID)
STFs_only=STFs_only[rownames(pcoa.axes),] ## reorder metadata.v2 to match order of CLR data

rownames(STFs_only) # check rownames to see if they are in the same order in both data frames after reordering
rownames(pcoa.axes)

# glm_<- vector('list', ncol(pcoa.axes) * ncol(STFs_only)) # create empty list where the GLM output is stored
# results_<- vector('list', ncol(pcoa.axes) * ncol(STFs_only)) # create an empty list where the GLM summaries are stored
# sig.results<-vector('list', ncol(pcoa.axes) * ncol(STFs_only))
# mdlnum <- 1 # counting our model numbers for indexes purposes in the loop
#
# # use a loop to run a bunch of GLMs
# ## pcoa.axes[i] is dependent variable (y), STFs_only[j] is independent variable (x) in GLM
# for (i in 1:ncol(pcoa.axes)){ # for each column in pcoa.axes
#   for (j in 1:ncol(STFs_only)){ # for each column in STFs_only
#     glm_[[mdlnum]] <-glm(pcoa.axes[,i]~STFs_only[,j], family=gaussian) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
#     results_[[mdlnum]] <-summary(glm_[[mdlnum]]) # save results of glm into list called results
#     names(results_)[mdlnum]<-paste(names(pcoa.axes)[i],"~",names(STFs_only)[j]) # rename list element to contain the name of the columns used in the model
#     mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)
#
#   }
# }
#
# ## pcoa.axes[i] is dependent variable (y), STFs_only[j] is independent variable (x) in GLM
# for (i in 1:ncol(pcoa.axes)){ # for each column in pcoa.axes
#   for (j in 1:ncol(STFs_only)){ # for each column in STFs_only
#     glm_[[mdlnum]] <-glm(pcoa.axes[,i]~STFs_only[,j], family=gaussian) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
#     results_[[mdlnum]] <-summary(glm_[[mdlnum]]) # save results of glm into list called results
#     names(results_)[mdlnum]<-paste(names(pcoa.axes)[i],"~",names(STFs_only)[j]) # rename list element to contain the name of the columns used in the model
#
#     ifelse(coef(results_[[mdlnum]])[,4] < 0.05, sig.results[[mdlnum]]<-results_[[mdlnum]], "Not Sig")
#     names(sig.results)[mdlnum]<-paste(names(pcoa.axes)[i],"~",names(STFs_only)[j])
#     mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)
#
#   }
# }
# sig.results[sapply(sig.results, is.null)] <- NULL

multi.univar.glm.fxn<-function(dep.var.df,indep.var.df,distro){
  # create empty lists to store stuff & model number (mdlnum) to keep track of models each iteration of loop in fxn
  glm_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create empty list where the GLM output is stored
  results_<- vector('list', ncol(dep.var.df) * ncol(indep.var.df)) # create an empty list where the GLM summaries are stored
  sig.results<-vector('list', ncol(dep.var.df) * ncol(indep.var.df))
  mdlnum <- 1 # counting our model numbers for indexes purposes in the loop

  # run the nested loop that generates GLMs from each data frame
  ## dep.var.df[i] is dependent variable (y), indep.var.df[j] is independent variable (x) in GLM
  for (i in 1:ncol(dep.var.df)){ # for each column in dep.var.df
    for (j in 1:ncol(indep.var.df)){ # for each column in indep.var.df
      glm_[[mdlnum]] <-glm(dep.var.df[,i]~indep.var.df[,j], family=distro) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
      results_[[mdlnum]] <-summary(glm_[[mdlnum]]) # save results of glm into list called results
      names(results_)[mdlnum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j]) # rename list element to contain the name of the columns used in the model

      # save only significant GLMs to another list called sig.results
      ## if p-value < 0.05, save to sig.results list
      ifelse(coef(results_[[mdlnum]])[,4] < 0.05, sig.results[[mdlnum]]<-results_[[mdlnum]], "Not Sig")
      names(sig.results)[mdlnum]<-paste(names(dep.var.df)[i],"~",names(indep.var.df)[j])
      mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)

    }
  }

  # drop all NULL elements from sig.results list so it only includes significant GLMs
  sig.results[sapply(sig.results, is.null)] <- NULL

  # assign lists to global env so they are saved there are function ends
  assign("results.glms", results_,envir = .GlobalEnv)
  assign("sig.results.glms", sig.results,envir = .GlobalEnv)

}

multi.univar.glm.fxn(pcoa.axes,STFs_only,gaussian) # test the function!

#### Save Everything ####
save.image("data/SSea Dust_BetaDiv_Data.Rdata")

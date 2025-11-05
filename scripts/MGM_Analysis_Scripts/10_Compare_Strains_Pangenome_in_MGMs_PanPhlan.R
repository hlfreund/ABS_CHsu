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

#### Import PanPhlan Results ####

# panphlan ecoli coverage results file
abs.ecoli.cov <- as.data.frame(read.delim('data/Metagenomes/Revisions_3.3.2025/PanPhlan_Results/ABS_panphlan_result_profile_ecoli_coverage.tsv', sep="\t"))
head(abs.ecoli.cov)
names(abs.ecoli.cov)<-gsub("_L00[0-9]_panphlan_ecoli.tsv","",names(abs.ecoli.cov)) # remove seq file info from sample ID
names(abs.ecoli.cov)<-gsub("_panphlan_ecoli.tsv","",names(abs.ecoli.cov)) # remove seq file info from sample ID

abs.ecoli.cov[1:4,1:4]
colnames(abs.ecoli.cov)[which(names(abs.ecoli.cov) == "X")] <- "Reference"
# drop A022_10_S12 sample - excluded from SRA upload
abs.ecoli.cov<-abs.ecoli.cov[,which(!names(abs.ecoli.cov)=="A022_10_S12")]
rownames(abs.ecoli.cov)<-abs.ecoli.cov$Reference

#### CLR Transform All Comp Data ####
rownames(abs.ecoli.cov)
abs.ecoli.cov[1:4,1:4]

abs.ecoli.cov.table<-as.data.frame(t(abs.ecoli.cov[,-1]))
abs.ecoli.cov.table[1:4,1:4]

# CLR transformation of ASV table
# df must have rownames are SampleIDs, columns are taxa IDs for vegan functions below
ecoli.clr<-decostand(abs.ecoli.cov.table,method = "clr", pseudocount = 1) #CLR transformation
ecoli.clr[1:4,1:4]
ecoli.clr$SampleID<-rownames(ecoli.clr)

ecoli.clr.melt<-melt(ecoli.clr,by="SampleID")
head(ecoli.clr.melt)
colnames(ecoli.clr)[which(names(ecoli.clr) == "variable")] <- "Reference"
colnames(ecoli.clr)[which(names(ecoli.clr) == "value")] <- "CLR_Coverage"

# merge CLR transformed data with metadata
ecoli.clr.meta<-merge(ecoli.clr,abs.all.basic.meta,by="SampleID")
# sanity checks: note that $SampleID and $col do not all contain same sample IDs due to nature of dataframe and sample to sample comparisons

ecoli.clr.meta$Condition<-factor(ecoli.clr.meta$Condition,levels=c("HC","Flare","Remission","HHP"))

#### Ecoli Diversity ####

# check rownames of CLR transformed ASV data & abs.all.basic.meta
rownames(ecoli.clr) %in% rownames(abs.all.basic.meta)
abs.all.basic.meta=abs.all.basic.meta[rownames(ecoli.clr),] ## reorder abs.all.basic.meta to match order of CLR data

# calculate our Euclidean distance matrix using CLR data
ecoli.euc_dist <- dist(ecoli.clr, method = "euclidean")
# ecoli.euc_dist2 <- vegdist(ecoli.clr, method = "euclidean") # did this to compare the two matrices to see ifthey're identical
# fviz_dist(ecoli.euc_dist, gradient = list(low = "blue", mid = "white", high = "red"))

# creating our hierarcical clustering dendrogram
ecoli.euc_clust <- hclust(ecoli.euc_dist, method="ward.D2")

# let's make it a little nicer...
ecoli.euc_dend <- as.dendrogram(ecoli.euc_clust, hang=0.02)
ecoli.dend_cols <- as.character(abs.all.basic.meta$Group_Color[order.dendrogram(ecoli.euc_dend)])
labels_colors(ecoli.euc_dend) <- ecoli.dend_cols

par(mar=c(1,1,1,1))
plot(ecoli.euc_dend, ylab="CLR Euclidean Distance", cex = 0.1,horiz=TRUE) + title(main = "Fungi Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
legend("topleft",legend = unique(abs.all.basic.meta$Condition),cex=.4,col = unique(abs.all.basic.meta$Group_Color),pch = 1, bty = "n")
dev.off()

# let's try another hierarchical cluster dendrogram with eclust()
# ecoli.clr.hc<-eclust(ecoli.clr, "hclust",nboot = 501) # use eclust() for "enhanced" hierarchical clustering; more here: https://www.datanovia.com/en/blog/cluster-analysis-in-r-simplified-and-enhanced/
#fviz_dend(ecoli.clr.hc, rect = TRUE) # plot the dendrogam

# PCOA w/ Euclidean distance matrix (of CLR data) aka Aitchison distance
ecoli.pcoa <- pcoa(ecoli.euc_dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("data/Amplicon/ABS_ITS2_CLR_EucDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
ecoli.pcoa$values

# extract principal coordinates
ecoli.pcoa.vectors<-data.frame(ecoli.pcoa$vectors)
ecoli.pcoa.vectors$SampleID<-rownames(ecoli.pcoa$vectors)

# merge pcoa coordinates w/ abs.all.basic.meta
ecoli.pcoa.meta<-merge(ecoli.pcoa.vectors, abs.all.basic.meta, by.x="SampleID", by.y="SampleID")

head(ecoli.pcoa.meta)
rownames(ecoli.pcoa.meta)<-ecoli.pcoa.meta$SampleID

head(ecoli.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 57.47%, PC2 = 17.37%

# Use princomp() to run PCA of Euclidean distances of CLR data (i.e., still Aitchison distances) so we can get to species contributions!
## create PCA of CLR-transformed counts (this will use Euclidean distance and gives the same results as the PCoA; you can compare based on the PC axes %)
ecoli.pca<-prcomp(ecoli.clr[,!colnames(ecoli.clr) %in% "SampleID"])
summary(ecoli.pca) # proportion of variance gives us variance explained by each axis; should match your pcoa results if you're using Aitchison distance
ecoli.pca$rotation # gives the contributions of each ASV to the variance/similarity of samples across the Euclidean space

# visualize PCA with contribution variables and the samples
## color scale corresponds to variable (here ASVs) contributions to PC axes' variation
# fviz_pca_var(ecoli.pca,
#              col.var = "contrib", # Color by contributions to the PC
#              gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"))+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,angle=45),panel.background = element_rect(fill = "white", colour = "white"),
#         plot.background = element_rect(fill = "white"))
# ggsave(filename = "figures/BetaDiversity/ABS_ITS2_CLR_PCA_ContributingASVs_biplot.png", 
#        width = 1800,height = 1200,units = "px",
#        dpi = 200,create.dir = TRUE)

# PCA biplot that shows contributing variables (ASVs) and samples in different colors
# fviz_pca_biplot(ecoli.pca,
#                 col.var = "contrib", # Variables color
#                 gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
#                 col.ind = "#696969") + theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#                                              axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,angle=45),panel.background = element_rect(fill = "white", colour = "white"),
#                                              plot.background = element_rect(fill = "white"))
# ggsave(filename = "figures/BetaDiversity/ABS_ITS2_CLR_PCA_ContributingASVs_Samples_biplot.png", 
#        width = 1800,height = 1200,units = "px",
#        dpi = 200,create.dir = TRUE)

# view the contributions of the variables; aka visualizing ecoli.pca$rotations in a meaningful way
# ## ASV contributions for PC1
# fviz_contrib(ecoli.pca,choice="var",fill="lightslateblue",color=NA,top=50, axes = 1)+theme_minimal() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,angle=45),panel.background = element_rect(fill = "white", colour = "white"),
#         plot.background = element_rect(fill = "white"))+
#   xlab("Top 50 Contributing ASVs - PC Axis 1") + ylab("Contributions (%)")
# ggsave(filename = "figures/BetaDiversity/ABS_ITS2_CLR_PCA_Top50ContributingASVs_PC1.png", 
#        width = 1800,height = 1200,units = "px",
#        dpi = 200,create.dir = TRUE)
# 
# ## ASV contributions for PC2
# fviz_contrib(ecoli.pca,choice="var",fill="lightslateblue",color=NA,top=50, axes = 2)+theme_minimal() +
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),
#         axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,angle=45),panel.background = element_rect(fill = "white", colour = "white"),
#         plot.background = element_rect(fill = "white"))+
#   xlab("Top 50 Contributing ASVs - PC Axis 2") + ylab("Contributions (%)")
# ggsave(filename = "figures/BetaDiversity/ABS_ITS2_CLR_PCA_Top50ContributingASVs_PC2.png", 
#        width = 1800,height = 1200,units = "px",
#        dpi = 200,create.dir = TRUE)

# K-means Clustering - will show us the groups that also appear in our PCoA
## more on K-means clustering: http://www.sthda.com/english/wiki/factoextra-r-package-easy-multivariate-data-analyses-and-elegant-visualization
## more on K-means clustering: https://uc-r.githuecoli.io/kmeans_clustering#:~:text=The%20gap%20statistic%20compares%20the,simulations%20of%20the%20sampling%20process.

## first we are going to play around with different values for k (aka # of clusters) and plot them to see what seems best
# k2 <- kmeans(ecoli.clr, centers = 2, nstart = 25)
# k3 <- kmeans(ecoli.clr, centers = 3, nstart = 25)
# k4 <- kmeans(ecoli.clr, centers = 4, nstart = 25)
# k5 <- kmeans(ecoli.clr, centers = 5, nstart = 25)
# 
# # plots to compare different values for k (aka # clusters)
# p1 <- fviz_cluster(k2, geom = "point",  data = ecoli.clr) + ggtitle("k = 2")
# p2 <- fviz_cluster(k3, geom = "point",  data = ecoli.clr) + ggtitle("k = 3")
# p3 <- fviz_cluster(k4, geom = "point",  data = ecoli.clr) + ggtitle("k = 4")
# p4 <- fviz_cluster(k5, geom = "point",  data = ecoli.clr) + ggtitle("k = 5")
# 
# # next we can calculate the K-means clustering and get back the ideal # of clusters
# ecoli.clr.km <- eclust(ecoli.clr, "kmeans", hc_metric="euclid",nboot = 1000)
# ecoli.clr.km
# 
# # Visualize kmeans clustering with output form eclust()
# # use repel = TRUE to avoid overplotting
# fviz_cluster(ecoli.clr.km, ecoli.clr, ellipse = TRUE,
#              ellipse.level = 0.95, ellipse.alpha = 0.2,ellipse.type = "convex",outlier.color = "black")
# 
# # Gap statistic plot to determine the ideal # of clusters k
# # The gap statistic compares the total intracluster variation for different values of k with their expected values under null reference distribution of the data (i.e. a distribution with no obvious clustering)
# fviz_gap_stat(ecoli.clr.km$gap_stat) # shows 1 cluster
# fviz_nbclust(ecoli.clr, kmeans, method = "gap_stat") # also shows 1 clusters
# fviz_nbclust(ecoli.clr, kmeans, method = "silhouette") # shows 2 clusters
# 
# # Optimal number of clusters using gap statistics = 1
# ecoli.clr.km$nbclust

#### Homogeneity of Variance & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(abs.all.basic.meta)
head(ecoli.clr)
rownames(abs.all.basic.meta) %in% rownames(ecoli.clr) #ecoli.clr was used to make the distance matrix ecoli.euc_dist

# first by compare dispersions by group
ecoli.disper1<-betadisper((vegdist(ecoli.clr[,!colnames(ecoli.clr) %in% "SampleID"],method="euclidean")), abs.all.basic.meta$Condition)
ecoli.disper1

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(ecoli.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#                Flare         HC        HHP Remission
# Flare                1.0000e-03 2.0000e-03     0.010
# HC        4.3332e-10            2.7000e-02     0.001
# HHP       2.3254e-04 2.5270e-02                0.105
# Remission 4.8400e-03 4.8680e-06 1.0768e-01          

anova(ecoli.disper1) # p = 4.003e-10 *** --> reject the Null H, spatial medians (a measure of dispersion) are significantly difference across flare groups
# ANOVA adjusted p-value
aov.beta.p1<-anova(ecoli.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(ecoli.disper1) # tells us which category's dispersion MEANS are significantly different than each other
#                     diff       lwr       upr     p adj
# HC-Flare        -104.13402 -141.21711 -67.05092 0.0000000
# HHP-Flare        -81.17182 -122.36074 -39.98289 0.0000063
# Remission-Flare  -56.30262  -95.81941 -16.78584 0.0017435
# HHP-HC            22.96220  -14.53449  60.45889 0.3854967
# Remission-HC      47.83139   12.17960  83.48319 0.0036418
# Remission-HHP     24.86919  -15.03598  64.77436 0.3696890

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(ecoli.clr[,!colnames(ecoli.clr) %in% "SampleID"] ~ Condition,data=abs.all.basic.meta,method = "euclidean",by="terms",permutations= 10000)
pnova1 # p-value = 9.999e-05
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 1

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

ecoli.clr.dist = (vegdist(ecoli.clr[,!colnames(ecoli.clr) %in% "SampleID"], "euclidean", na.rm = TRUE)) #distance matrix using Bray's dissimilarity index for trait distribution (traits of interest only)
pair.mod1<-pairwise.adonis(ecoli.clr.dist,abs.all.basic.meta$Condition, p.adjust.m='bonferroni',perm=9999) # shows us variation for each sample to see which ones are different
pair.mod1
#                 pairs Df SumsOfSqs  F.Model        R2     p.value     p.adjusted sig
# 1 Flare vs Remission  1  42430.60  3.328240 0.05426929  0.0216     0.1296    
# 2       Flare vs HHP  1  71383.54  5.761841 0.09805413  0.0015     0.0090   *
# 3        Flare vs HC  1 127981.23 15.530447 0.18372607  0.0001     0.0006  **
# 4   Remission vs HHP  1  17569.99  3.139034 0.05219629  0.0351     0.2106    
# 5    Remission vs HC  1  39939.42 12.619803 0.14739350  0.0001     0.0006  **
# 6          HHP vs HC  1   5745.42  2.635735 0.03731447  0.0393     0.2358    

# Visualize dispersions
par(mar=c(1,1,1,1))
plot(ecoli.disper1,main = "Centroids and Dispersion based on Aitchison Distance")
dev.off()

par(mar=c(1,1,1,1))
boxplot(ecoli.disper1,xlab="By Flare", main = "Distance to Centroid by Category", sub="Based on Aitchison Distance")
dev.off()


#### Visualize PCoAs - All Data ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints

# create PCoA ggplot fig
pval<-"p-value = 0.0003" # pvalue from PERMANOVA below that I already ran comparing Flare groups
head(ecoli.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
ecoli.pcoa.meta$Condition<-factor(ecoli.pcoa.meta$Condition,levels=c("HC","Flare","Remission","HHP"))

pcoa1<-ggplot(ecoli.pcoa.meta, aes(x=Axis.1, y=Axis.2,color=Condition))+ annotate("text", x = 40, y = 205, label = pval,size=3.5) +geom_point(size=3)+theme_bw()+
  labs(color="Flare Group")+theme_classic()+stat_ellipse(geom="polygon",aes(x=Axis.1, y=Axis.2,fill=Condition),linetype=1,alpha = 0.25)+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),legend.title = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1),legend.text = element_text(size=11))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Flare Group",values=unique(ecoli.pcoa.meta$Group_Color[order(ecoli.pcoa.meta$Condition)])) +
  scale_fill_manual(name ="Flare Group",values=unique(ecoli.pcoa.meta$Group_Color[order(ecoli.pcoa.meta$Condition)])) +
  xlab("PC1 [57.43%]") + ylab("PC2 [17.39%]")

ggsave(pcoa1,filename = "figures/BetaDiversity/ABS_ITS2_CLR_FlareGroup_Ellipses_Pval_PCOA1.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa1a<-ggplot(ecoli.pcoa.meta, aes(x=Axis.1, y=Axis.2,color=Condition))+ annotate("text", x = 200, y = 205, label = pval,size=5) +geom_point(size=4)+theme_bw()+
  labs(color="Flare Group")+theme_classic()+
  theme(axis.title.x = element_text(size=20),axis.title.y = element_text(size=20),legend.title = element_text(size=18),
        axis.text = element_text(size=17),axis.text.x = element_text(vjust=1),legend.text = element_text(size=17))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Condition",values=unique(ecoli.pcoa.meta$Group_Color[order(ecoli.pcoa.meta$Condition)])) +
  scale_fill_manual(name ="Condition",values=unique(ecoli.pcoa.meta$Group_Color[order(ecoli.pcoa.meta$Condition)])) +
  xlab("PC1 [57.43%]") + ylab("PC2 [17.39%]")

ggsave(pcoa1a,filename = "figures/PanPhlan/ABS_PanPhlan_EcoliOnly_CLR_FlareGroup_Pval_PCOA1.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)


# 3D PCoA

pltly.all.a<-plot_ly(ecoli.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~Flare, colors = c(unique(ecoli.pcoa.meta$Group_Color[order(ecoli.pcoa.meta$Condition)])),
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

save_image(pltly.all.a, "figures/BetaDiversity/ABS_ITS2_CLR_Flare_3D_Aitchison_PCOA1.png",width=1200,height=1000)
save_image(pltly.all.a, "figures/BetaDiversity/ABS_ITS2_CLR_Flare_3D_Aitchison_PCOA2.png",width=1400,height=1100)


# save 3D plot as an HTml
saveWidget(widget = pltly.all.a, #the plotly object,
           file = "figures/BetaDiversity/ABS_ITS2_CLR_Flare_3D_Aitchison_PCOA1.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

# ## Loop to Generate Heat Map for Each PC Axis
# 
# pc.plot.list<-list() # create empty list for each plot to be stored in
# pc.axes<-names(ecoli.pcoa.dts)[grepl("Axis",names(ecoli.pcoa.dts))] # pull out names of columns in df that contain "Axis" in name
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
#   pc.heatmap=hm.fxn(ecoli.pcoa.dts, ecoli.pcoa.dts$Condition, ecoli.pcoa.dts$Condition, ecoli.pcoa.dts[,i])
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
#   #ggsave(a, file = paste0("figures/BetaDiversity/ABS_ITS2.PCoA_heatmap_", f_var,".png"), device = png, width = 15, height = 15, units = "cm")
#
#   return(a)
# }
# # loop with heatmap function to create heatmap
# 
# 
# s.t.pcoa<-melt(ecoli.pcoa.dts[,-1], by=c("Flare","SampDate"))
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

#ggsave(sulecoli.hm1a,filename = "figures/MGM_Figs/FxnDiv/Sulfur/Sulfur_KOFxns_MGMs_SampID_by_Function_heatmap.png", width=18, height=13, dpi=600,create.dir = TRUE)

#### Session Info ####
sessionInfo() # see loaded packages

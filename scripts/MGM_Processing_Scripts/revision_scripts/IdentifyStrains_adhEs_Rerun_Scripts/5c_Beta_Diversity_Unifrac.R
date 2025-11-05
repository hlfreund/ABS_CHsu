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
  library(shades)
  #library(ALDEx2)
  library(rstatix)
  #library(decontam)
  #library(ggvegan)
  library(microbiome)
  library(pairwiseAdonis)
  library(corrplot)
  library(fst)
  library(NatParksPalettes)
  library(plotly)
  library(htmlwidgets)
  library(rbiom)
  library(Maaslin2)
  library(car) # car is imported w/ rstatix so you do not need to load library individually?
  library(factoextra)
})

## NOTE FOR THIS SCRIPT:
# In order to generate Unifrac distances, you need to have created a Multiple Sequence Alignment & built a Phylogenetic Tree
## if you do not have this yet, go back and run 5a_Prep_Data_for_MSA_PhylogeneticTree.R, then 5b_MultipleSequenceAlignment_and_PhylogeneticTree.sh
## to run this script you need the FastTree output: *.bacteria.FastTree.tre

#### Load Global Env to Import Count/ASV Tables ####
load("data/DADA2_Results/FT_16SV3V4_DataReady.Rdata") # save global env to Rdata file

bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(metadata)

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
colorset1 # group colors
# for your reference:
# WT="#D4D4D4",MUC2="#E9BE41",KO="#90BFF9"
# WT="#D4D4D4",MUC2="#E9BE41",KO="#90BFF9"

#### List of Important R Objects for Reference ####

bac.ASV_counts[1:5,1:5] # original ASV counts object, singletons removed, used to create ASV table
bac.ASV_table[1:5,1:5] # this is the updated ASV table where singletons have been removed
colnames(bac.ASV_table) ## bac.ASV_table has sample names as rows, ASVs as column names

head(metadata) # updated metadata that includes GroupName group colors
bac.ASV_meta[1:5,] # combined metadata + ASV counts

head(bac.ASV_tax) # updated taxa table with singleton ASVs removed
head(bac.ASV_all.tax) # updated taxa & ASV counts, not melted (columns include individual samples)
head(bac.AllTax.melt) # melted combined ASVs + taxonomy table -- for relative abundance calcs
## ^ has melted column of SampleID

head(bac.ALL.dat) # ALL data combined & melted
# ^ metadata, ASVs + counts, & taxonomy data

head(new.asv_fasta) # udpated FASTA file to be used for Unifrac distances of 16S V3V4
## ^ excludes contaminant and/or singleton ASVs! updated fasta file from the output from DADA2

# #### Calculate Relative Abundance of ASVs ####
# bac.ASV_table[1:4,1:4]
# which(colnames(bac.ASV_table) %in% "SampleID")
# 
# b.ASV_RelAb<-data.frame(decostand(bac.ASV_table[,-c(3184)], method="total", MARGIN=1, na.rm=TRUE))
# # relative abundance of taxa data where everything is divided by col total (b/c Margin=1 meaning rows == SAMPLES in this case)
# rowSums(b.ASV_RelAb) # sanity check to make sure the transformation worked!
# 
# b.RelAb_mat.t<-as.matrix(t(b.ASV_RelAb))
# b.RelAb_mat.t[1:4,1:4]
# 
# #### Arcsin(SqRt()) Transformation of All Pathway Relative Abundance Data ####
# 
# # now we are going to use the Arcsin(Square Root()) Transformation to transform the pathway relative abundances
# ## sources for this transformation: LLorens-Rico et al 2021
# ## Arcsin(sqrt()) transform is ideal for proportions or percentages data (i.e., what we have)
# ## NOTE: even if your data is proportional, it must be between 0-1 for the arcsin(sqrt()) transformation to work!
# 
# b.RelAb.asinsqrt<-asin(sqrt(b.ASV_RelAb)) # Arcsin(sqrt()) transformation of the pathway rel abs
# b.RelAb.asinsqrt[1:4,1:4]
# 
# b.RA.asinsqrt_mat.t<-as.matrix(t(b.RelAb.asinsqrt))
# b.RA.asinsqrt_mat.t[1:4,1:4]
# 
#### Rarefaction of Raw ASV Counts ####
# this section explores different ways of utilizing rarefaction and rarefaction functions in vegan

# in vegan ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs)
min.rar<-min(rowSums(bac.ASV_table[,-c(3184)])) ## seeing min sum of OTUs so we can see what min is for rarefaction
min.rar

bac.ASV.rar<-vegan::rrarefy(bac.ASV_table[,-c(3184)],min.rar) ## be cognizant of min for rarefaction
bac.ASV.rar 

# create transposed matrix for creating unifrac distances later
b.rar_mat.t<-as.matrix(t(bac.ASV.rar))
b.rar_mat.t[1:4,1:4]

#### Import Phylogenetic Tree File ####

phylo.tree<-read_tree("data/FT_16SV3V4_MSA_ssualign.bacteria.bacteria.FastTree.tre")
# ^ this object is class "phylo"
phylo.tree$tip.label[length(phylo.tree$tip.label)] # shows us the last tip in the tree, could use as an outgroup
# ASV 1632 will serve as outgroup for rooting tree to get Unifrac distances
phylo.tree$tip.label[length(phylo.tree$tip.label)-1] # shows us the second to last tip in the tree, could use as an outgroup

phylo.tree$tip.label[1] # first tip in tree

# use the ASV outgroup to root your phylogenetic tree -- need this for calculating Unifrac distance
phylo.tree.rooted<-root(phylo.tree,outgroup="ASV 1632", resolve.root = TRUE)
is.rooted(phylo.tree.rooted) # should be TRUE

phylo.tree.rooted$tip.label<-gsub("ASV ","ASV_",phylo.tree.rooted$tip.label)
phylo.tree.rooted$tip.label

#cophenetic(phylo.tree)
#dist.nodes(phylo.tree)

#### Drop ASVs from Table that are NOT in Phylo Tree ####

# first from transformed ASV tables
b.rar_mat.t[1:4,1:4]
phylo.tree.rooted$tip.label

# subset rows that have rownames found in the phylogenetic tree tip labels
b.rar_mat.t<-b.rar_mat.t[rownames(b.rar_mat.t) %in% phylo.tree.rooted$tip.label,]

# check dimensions and if they're the same...
## # of rows in b.rar_mat.t should match # of tips in phylo.tree.rooted
dim(as.data.frame(phylo.tree.rooted$tip.label))
dim(b.rar_mat.t)

# then from raw ASV table which we will use for Unweighted Unifrac Distance
bac.ASV_table[1:4,1:4]
bac.ASV_mat.t<-as.matrix(t(bac.ASV_table[,-c(3184)])) # transpose ASV table and turn into a matrix
bac.ASV_mat.t<-bac.ASV_mat.t[rownames(bac.ASV_mat.t) %in% phylo.tree.rooted$tip.label,] # drop ASVs not in tree

#### Calculate Weighted Unifrac Distance ####

# using rarefied data as input for weighted unifrac distances
b.wunifrac.dist<-bdiv_distmat(b.rar_mat.t, bdiv="UniFrac", weighted = TRUE, tree = phylo.tree.rooted, transform="none")
# weighted Unifrac considers the ASV abundances, non-weighted unifrac looks at presence/absence data

class(b.wunifrac.dist) #check if this is now class dist, and it is!

#### Beta Diversity - Weighted Unifrac ####
metadata=metadata[colnames(b.rar_mat.t),] ## reorder metadata to match order of Rarefied data used to generate Unifrac distances

# creating our hierarcical clustering dendrogram
b.wunifrac_clust <- hclust(b.wunifrac.dist, method="ward.D2")

# let's make it a little nicer...
b.wunifrac_dend <- as.dendrogram(b.wunifrac_clust, hang=0.02)
b.dend_cols <- as.character(metadata$GroupName_Color[order.dendrogram(b.wunifrac_dend)])
labels_colors(b.wunifrac_dend) <- b.dend_cols

colorset1 # color dendrogram by group
#png(filename="figures/BetaDiversity/WeightedUnifrac/FT_16SV3V4_Rarefied_WeightedUnifracDist_Dendrogram1.png",width = 7, height = 7, units = "in",res = 800)
par(mar=c(1,1,1,1))
plot(b.wunifrac_dend, ylab="Rarefied Weighted Unifrac Distance", cex = 0.1,horiz=TRUE) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topleft",legend = colorset8$GroupName,cex=.8,col = colorset8$GroupName_Color,pch = 15, bty = "n")
#dev.off()

# PCOA w/ Weighted Unifrac distance matrix (of Rarefied data)
b.wunifrac.pcoa <- pcoa(b.wunifrac.dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("FT_16SV3V4_Rarefied_WeightedUnifracDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.wunifrac.pcoa$values

# extract principal coordinates
b.wunifrac.pcoa.vectors<-data.frame(b.wunifrac.pcoa$vectors)
b.wunifrac.pcoa.vectors$SampleID<-rownames(b.wunifrac.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.wunifrac.pcoa.meta<-merge(b.wunifrac.pcoa.vectors, metadata, by.x="SampleID", by.y="SampleID")

head(b.wunifrac.pcoa.meta)
rownames(b.wunifrac.pcoa.meta)<-b.wunifrac.pcoa.meta$SampleID

head(b.wunifrac.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 27.70%, PC2 = 16.73%
# save.image("data/FT_16SV3V4_WeightedUnifracDist_Ready.Rdata")

#### Visualize PCoAs - Weighted Unifrac ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints

# create PCoA ggplot fig

## quick permanova to get p value for fig
pnova1<-adonis2(b.wunifrac.dist ~ GroupName,data=metadata,by="terms",permutations=10000)
pnova1 # p-value = 0.9082

pcoa1<-ggplot(b.wunifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=GroupName),size=4)+theme_bw()+
  labs(color="Sample Date")+theme_classic()+ 
  theme(axis.title.x = element_text(size=23),axis.title.y = element_text(size=23),legend.title.align=0.5, legend.title = element_text(size=23),axis.text = element_text(size=20),axis.text.x = element_text(vjust=1),legend.text = element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Group",values=unique(b.wunifrac.pcoa.meta$GroupName_Color[order(b.wunifrac.pcoa.meta$GroupName)])) +
  xlab("PC1 [27.70%]") + ylab("PC2 [16.73%]") + annotate("text", x=0.06,y=0.5,label="PERMANOVA, p = 0.91")

ggsave(pcoa1,filename = "figures/BetaDiversity/WeightedUnifrac/FT_16SV3V4_Rarefied_GroupName_WeightedUnifrac_PCOA1.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa2<-ggplot(b.wunifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2))+
  geom_text(label=b.wunifrac.pcoa.meta$SampleID,color=b.wunifrac.pcoa.meta$GroupName_Color,position="identity") +
  theme_bw()+
  labs(title="PCoA: Bacteria/Archaea",subtitle="Using Weighted Unifrac Distance of Rarefied Counts")+theme_classic()+ theme(axis.title.x = element_text(size=23),axis.title.y = element_text(size=23),legend.title.align=0.5, legend.title = element_text(size=23),axis.text = element_text(size=20),axis.text.x = element_text(vjust=1),legend.text = element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("PC1 [27.70%]") + ylab("PC2 [16.73%]")

ggsave(pcoa2,filename = "figures/BetaDiversity/WeightedUnifrac/FT_16SV3V4_Rarefied_GroupName_Labeled_WeightedUnifrac_PCOA2.png", width=14, height=10, dpi=600, create.dir=TRUE)

# 3D PCoA

pltly.all.a<-plot_ly(b.wunifrac.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~GroupName, colors = c("#D4D4D4","#E9BE41","#90BFF9"),
                     symbol=~GroupName,symbols = c("square", "circle","diamond")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 27.70%'),
                      yaxis = list(title = 'PC2 16.73%'),
                      zaxis = list(title = 'PC3 14.66%')))

saveWidget(widget = pltly.all.a, #the plotly object,
           file = "figures/BetaDiversity/WeightedUnifrac/FT_16SV3V4_WeightedUnifrac_GroupName_3D_PCOA1.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)

#### Calculate Unweighted Unifrac Distance ####
bac.ASV_mat.t[1:4,1:4]

b.unifrac.dist<-bdiv_distmat(bac.ASV_mat.t, bdiv="UniFrac", weighted = FALSE, tree = phylo.tree.rooted, transform="none")
# weighted Unifrac considers the ASV abundances, non-weighted unifrac looks at presence/absence data

class(b.unifrac.dist) #check if this is now class dist, and it is!

#### Beta Diversity - Unweighted Unifrac ####
metadata=metadata[colnames(bac.ASV_mat.t),] ## reorder metadata to match order of Rarefied data used to generate Unifrac distances

# creating our hierarcical clustering dendrogram
b.unifrac_clust <- hclust(b.unifrac.dist, method="ward.D2")

# let's make it a little nicer...
b.unifrac_dend <- as.dendrogram(b.unifrac_clust, hang=0.02)
b.dend_cols <- as.character(metadata$GroupName_Color[order.dendrogram(b.unifrac_dend)])
labels_colors(b.unifrac_dend) <- b.dend_cols

colorset1 # color dendrogram by collection date
#png(filename="figures/BetaDiversity/UnweightedUnifrac/FT_16SV3V4_Rarefied_WeightedUnifracDist_Dendrogram1.png",width = 7, height = 7, units = "in",res = 800)
par(mar=c(1,1,1,1))
plot(b.unifrac_dend, ylab="Unweighted Unifrac Distance", cex = 0.1,horiz=TRUE) + title(main = "Bacteria/Archaea Clustering Dendrogram", cex.main = 1, font.main= 1, cex.sub = 0.8, font.sub = 2)
#legend("topleft",legend = colorset8$GroupName,cex=.8,col = colorset8$GroupName_Color,pch = 15, bty = "n")
#dev.off()

# PCOA w/ Weighted Unifrac distance matrix (of Rarefied data)
b.unifrac.pcoa <- pcoa(b.unifrac.dist) # pcoa of euclidean distance matrix = PCA of euclidean distance matrix
#save.image("FT_16SV3V4_Rarefied_WeightedUnifracDist_Ready.Rdata")

# The proportion of variances explained is in its element values$Relative_eig
b.unifrac.pcoa$values

# extract principal coordinates
b.unifrac.pcoa.vectors<-data.frame(b.unifrac.pcoa$vectors)
b.unifrac.pcoa.vectors$SampleID<-rownames(b.unifrac.pcoa$vectors)

# merge pcoa coordinates w/ metadata
b.unifrac.pcoa.meta<-merge(b.unifrac.pcoa.vectors, metadata, by.x="SampleID", by.y="SampleID")

head(b.unifrac.pcoa.meta)
rownames(b.unifrac.pcoa.meta)<-b.unifrac.pcoa.meta$SampleID

head(b.unifrac.pcoa$values) # pull out Relative (Relative_eig) variation % to add to axes labels
# PC1 = 16.90%, PC2 = 10.47%, PC3 = 8.26%

#save.image("data/FT_16SV3V4_UnweightedUnifracDist_Ready.Rdata")

#### Visualize PCoAs - Unweighted Unifrac ####

## *** all figures that end in _PCOA1 come from the same single PCoA
#data is just subsetted to understand how points are related to each other w/in & across timepoints

# create PCoA ggplot fig
# quick permanova to get p-value
pnova2<-adonis2(b.unifrac.dist ~ GroupName,data=metadata,by="terms",permutations=10000)
pnova2 # p-value = 0.5536

pcoa1<-ggplot(b.unifrac.pcoa.meta, aes(x=Axis.1, y=Axis.2)) +geom_point(aes(color=GroupName),size=4)+theme_bw()+
  labs(color="Group Name")+theme_classic()+ 
  theme(axis.title.x = element_text(size=23),axis.title.y = element_text(size=23),legend.title.align=0.5, legend.title = element_text(size=23),axis.text = element_text(size=20),axis.text.x = element_text(vjust=1),legend.text = element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(name ="Group",values=unique(b.unifrac.pcoa.meta$GroupName_Color[order(b.unifrac.pcoa.meta$GroupName)])) +
  xlab("PC1 [16.90%]") + ylab("PC2 [10.47%]") + annotate("text", x=0.3,y=0.21,label="PERMANOVA, p = 0.55")

ggsave(pcoa1,filename = "figures/BetaDiversity/UnweightedUnifrac/FT_16SV3V4_GroupName_UnweightedUnifrac_PCOA1.png", 
       width = 1800,height = 1200,units = "px",
       dpi = 200,create.dir = TRUE)

pcoa2<-ggplot(b.unifrac.pcoa.meta, aes(x=Axis.2, y=Axis.3))+theme_bw()+
  geom_text(label=b.unifrac.pcoa.meta$SampleID,color=b.unifrac.pcoa.meta$GroupName_Color,position="identity") +
  labs(title="PCoA: Bacteria/Archaea",subtitle="Using Unweighted Unifrac Distance")+theme_classic()+ 
  theme(axis.title.x = element_text(size=23),axis.title.y = element_text(size=23),legend.title.align=0.5, legend.title = element_text(size=23),axis.text = element_text(size=20),axis.text.x = element_text(vjust=1),legend.text = element_text(size=20))+
  guides(shape = guide_legend(override.aes = list(size = 5)))+
  xlab("PC1 [16.90%]") + ylab("PC2 [10.47%]")

ggsave(pcoa2,filename = "figures/BetaDiversity/UnweightedUnifrac/FT_16SV3V4_GroupName_UnweightedUnifrac_PCOA2.png", width=14, height=10, dpi=600, create.dir=TRUE)

# 3D PCoA

pltly.all.a<-plot_ly(b.unifrac.pcoa.meta, x=~Axis.1,y=~Axis.2,z=~Axis.3, color = ~GroupName, colors = c("#D4D4D4","#E9BE41","#90BFF9"),
                     symbol=~GroupName,symbols = c("square", "circle","diamond")) %>%
  layout(scene = list(xaxis = list(title = 'PC1 16.90%'),
                      yaxis = list(title = 'PC2 10.47%'),
                      zaxis = list(title = 'PC3 8.26%')))

saveWidget(widget = pltly.all.a, #the plotly object,
           file = "figures/BetaDiversity/UnweightedUnifrac/FT_16SV3V4_UnweightedUnifrac_GroupName_3D_PCOA1.html", #the path & file name
           selfcontained = TRUE #creates a single html file
)


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

head(metadata)
head(bac.ASV.rar)
rownames(metadata) %in% rownames(as.matrix(b.wunifrac.dist)) #b.clr was used to make the distance matrix b.wunifrac.dist

# first by compare dispersions by sampling date
b.disper1<-betadisper(b.wunifrac.dist, metadata$GroupName)
b.disper1
b.disper1$distances

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper1, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#         WT    MUC2    KO
# WT           0.90900 0.771
# MUC2 0.90002         0.787
# KO   0.76050 0.75740      

anova(b.disper1) # p = 0.9346 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sites
# ANOVA adjusted p-value
aov.beta.p1<-anova(b.disper1)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(b.disper1) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova1<-adonis2(b.wunifrac.dist ~ GroupName,data=metadata,by="terms",permutations=10000)
pnova1 # p-value = 0.8831
p.adjust(pnova1$`Pr(>F)`,method="bonferroni",n=length(pnova1$`Pr(>F)`)) # adjusted pval
# 1

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod1<-pairwise.adonis(b.wunifrac.dist,metadata$GroupName, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod1
#       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 MUC2 vs KO  1 0.08406535 0.5709208 0.05400861   0.957          1    
# 2 MUC2 vs WT  1 0.11430912 0.7111178 0.07322718   0.841          1    
# 3   KO vs WT  1 0.10870028 0.8575492 0.07232095   0.616          1 
# Visualize dispersions
plot(b.disper1,main = "Centroids and Dispersion based on Weighted Unifrac Distance", labels=FALSE,col=colorset1$GroupName_Color)
dev.off()

par(mar=c(1,1,1,1))
boxplot(b.disper1,xlab="By Group Name", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset1$GroupName_Color)
dev.off()

# below we make the same plot as above but with scatter plot instead of boxplot

b.disper1.distances<-data.frame(Disper.Dist=b.disper1$distances)
b.disper1.distances$SampleID<-rownames(b.disper1.distances)
b.dipser1.dist.meta<-merge(metadata,b.disper1.distances,by="SampleID")

ggplot(b.dipser1.dist.meta,aes(x=GroupName, y=Disper.Dist, col=GroupName))+
  geom_jitter(aes(color=factor(GroupName)), size=5, width=0.15, height=0) +
  theme(axis.title.x = element_text(size=23),axis.title.y = element_text(size=23),
        axis.text = element_text(size=20),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=23),legend.text = element_text(size=20),plot.title = element_text(size=15)) +
  stat_summary(fun.data=median_hilow,color="black",size=1,fun.args=list(conf.int=0.5),show.legend=FALSE) + ylab("Distance to the Centroid")


#### Homogeneity of Variance of Unweighted Unifrac & PERMANOVA tests - Composition by Groups ####
## betadisper to look at homogeneity of group dispersions (aka variance) when considering multiple variables
# multivariate analogue to Levene's test of homogeneity of variances
# program finds spatial median or centroid of the group, & compare distances of group to centroid/spatial median via ANOVA

#While PERMANOVA tests differences in group means (analogous to MANOVA),
## a related test called PERMDISP can be used to evaluate homogeneity of group dispersion
#(analogous to Levene's test for equal variances). The vegan function for this test is “betadisper”:
## * need a distance matrix!

head(metadata)
head(bac.ASV_table)
rownames(metadata) %in% rownames(as.matrix(b.unifrac.dist)) #b.clr was used to make the distance matrix b.unifrac.dist

# first by compare dispersions by sampling date
b.disper2<-betadisper(b.unifrac.dist, metadata$GroupName)
b.disper2
b.disper2$distances

## Significant differences in homogeneities can be tested using either parametric or permutational tests,
##and parametric post hoc contrasts can also be investigated:

permutest(b.disper2, pairwise=TRUE) # compare dispersions to each other via permutation test to see significant differences in dispersion by pairwise comparisons
#Pairwise comparisons:
#  (Observed p-value below diagonal, permuted p-value above diagonal)
#         WT    MUC2    KO
# WT           0.90900 0.771
# MUC2 0.90002         0.787
# KO   0.76050 0.75740      

anova(b.disper2) # p = 0.8353 --> accept the Null H, spatial medians (a measure of dispersion) are NOT significantly difference across sites
# ANOVA adjusted p-value
aov.beta.p1<-anova(b.disper2)[["Pr(>F)"]] # get p values from ANOVA
p.adjust(aov.beta.p1,method="bonferroni",n=length(aov.beta.p1))

TukeyHSD(b.disper2) # tells us which Sample Dates/category's dispersion MEANS are significantly different than each other
#         diff       lwr       upr     p adj

# If PERMANOVA is significant but betadisper() IS NOT, then you can infer that there is only a location effect.
# If both tests are significant, then there is a dispersion effect for sure and there might also be (not always) a location effect.
# Dispersion effect means the actual spread of the data points is influencing the significant differences, not the actual data itself

pnova2<-adonis2(b.unifrac.dist ~ GroupName,data=metadata,by="terms",permutations=10000)
pnova2 # p-value = 0.5544
p.adjust(pnova2$`Pr(>F)`,method="bonferroni",n=length(pnova2$`Pr(>F)`)) # adjusted pval
# 1

##one issue with adonis is that it doesn't do multiple comparisons *******
# tells us that something is different, but what is different? Which sample/plot/location?
## our four provinces differ, but do all of them differ,or just one?

##random person on the internet to the rescue!
#install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")

pair.mod2<-pairwise.adonis(b.unifrac.dist,metadata$GroupName, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
pair.mod2
#       pairs Df SumsOfSqs   F.Model         R2 p.value p.adjusted sig
# 1 MUC2 vs KO  1 0.2823718 0.8053136 0.07452941   0.868      1.000    
# 2 MUC2 vs WT  1 0.3378840 0.9784181 0.09805343   0.505      1.000    
# 3   KO vs WT  1 0.3786314 1.0818067 0.08954014   0.279      0.837  

# Visualize dispersions
par(mar=c(1,1,1,1))
plot(b.disper2,main = "Centroids and Dispersion based on Weighted Unifrac Distance", labels=FALSE,col=colorset1$GroupName_Color)
dev.off()

par(mar=c(1,1,1,1))
boxplot(b.disper2,xlab="By Group Name", main = "Distance to Centroid by Category", sub="Based on Weighted Unifrac Distance", col=colorset1$GroupName_Color)
dev.off()

# below we make the same plot as above but with scatter plot instead of boxplot

b.disper2.distances<-data.frame(Disper.Dist=b.disper2$distances)
b.disper2.distances$SampleID<-rownames(b.disper2.distances)
b.dipser2.dist.meta<-merge(metadata,b.disper2.distances,by="SampleID")

ggplot(b.dipser2.dist.meta,aes(x=GroupName, y=Disper.Dist, col=GroupName))+
  geom_jitter(aes(color=factor(GroupName)), size=5, width=0.15, height=0) +
  theme(axis.title.x = element_text(size=23),axis.title.y = element_text(size=23),
        axis.text = element_text(size=20),axis.text.x = element_text(vjust = 1, hjust=1),legend.title.align=0.5, legend.title = element_text(size=23),legend.text = element_text(size=20),plot.title = element_text(size=15)) +
  stat_summary(fun.data=median_hilow,color="black",size=1,fun.args=list(conf.int=0.5),show.legend=FALSE) + ylab("Distance to the Centroid")


#### Species Contributions to Dissimilarity with SIMPER ####
# NOTES: (some from here https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/simper/)
# When there are multiple sample units, there are also multiple dissimilarities (i.e., pairwise combinations) to consider
# each species can contribute to the dissimilarity between each pair of sample units.
## If our purpose is to quantify the contribution of each species to the differences between two groups, we have to consider all of these dissimilarities.
# The average dissimilarity for each pairwise combination can be calculated directly via the vegan::meandist() function

meandist(dist = b.wunifrac.dist, grouping = metadata$GroupName)
# The distances shown here are a symmetric square matrix where each value is the average distance between two sample units.
# Values on the diagonal are the average distances between observations in the same group;
# all other values are the average distances between an observation in one group and an observation in another group.

# now for SIMPER analysis...
# Similarity percentage (SIMPER) partitions the dissimilarity matrix for every pair of sample units, and then calculates the average contribution of each species to the difference between the sample units.
# These contributions are relativized so that the average contributions of all species sum to 1.
# Statistical significance of these contributions is assessed by permuting the group identities.
# Note that the title is misleading: a high value for SIMPER means that a species has a high contribution to the difference between the two groups, NOT that it has high similarity between them!

simper1<-simper(b.clr,
                metadata$GroupName,
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
  ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$GroupName names
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
                interaction(metadata$GroupName,metadata$CollectionYear),
                permutations = 999
)


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

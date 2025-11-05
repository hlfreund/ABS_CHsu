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
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/ITS2_DADA2/ABS_ITS2_DataReady.Rdata") # save global env to Rdata file

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

head(metadata.v2) # updated metadata.v2 that includes Flare group colors and does not have Donor sample!
its2.ASV_meta[1:5,1:10] # combined metadata.v2 + ASV counts

head(its2.ASV_tax.clean) # updated taxa table with singleton ASVs removed
head(its2.AllTax.melt) # melted ASVs + taxonomy table -- for relative abundance calcs

head(its2.ALL.dat) # ALL data combined
# ^ metadata.v2, ASVs + counts, & taxonomy data

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
f.phyla_melt_meta.v2$Phylum<-factor(f.phyla_melt_meta.v2$Phylum,levels=phyla.list,ordered=TRUE)
f.phyla_melt_meta.v2<-f.phyla_melt_meta.v2[order(f.phyla_melt_meta.v2$Phylum,f.phyla_melt_meta.v2$Count),]

#### Using Shapiro-Wilk test for Normality - Rarefied Data ####
# use list of fungal phyla from most to least abundant to go through Shapiro-Wilks & aov() tests
phyla.list<-names(sort(colSums(f.phyla_RelAb[,-ncol(f.phyla_RelAb)]),decreasing=TRUE)) # see most abundant phyla
phyla.list

## only checking the two most abundant phyla - both were non-normal so we will use KruskalWallis, Wilcoxon tests
shapiro.test(f.phyla_RA_meta.v2$Ascomycota) # what is the p-value?
# W = 0.36004, p-value = 5.747e-15
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(f.phyla_RA_meta.v2$Ascomycota, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(f.phyla_RA_meta.v2$Ascomycota, pch = 1, frame = FALSE)
qqline(f.phyla_RA_meta.v2$Ascomycota, col = "red", lwd = 2)

shapiro.test(f.phyla_RA_meta.v2$Basidiomycota) # what is the p-value?
# W = 0.81983, p-value = 0.0002437
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(f.phyla_RA_meta.v2$Basidiomycota, col="blue")

# visualize Q-Q plot for species richness
qqnorm(f.phyla_RA_meta.v2$Basidiomycota, pch = 1, frame = FALSE) # with outliars
qqline(f.phyla_RA_meta.v2$Basidiomycota, col = "red", lwd = 2)

#### Compare Variance by Phylum ####
phyla.list
# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

# first, Ascomycota
# Kruskal-Wallis test is an ANOVA for non-normal data
fit1<-kruskal.test(Ascomycota ~ Flare, data=f.phyla_RA_meta.v2)
fit1
# data:  Ascomycota by Flare
# Kruskal-Wallis chi-squared = 1.9764, df = 2, p-value = 0.3722

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(f.phyla_RA_meta.v2, Ascomycota ~ Flare, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Ascomycota ~ Flare, data = f.phyla_RA_meta.v2)
# Fligner-Killeen:med chi-squared = 2.0941, df = 2, p-value = 0.351
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Ascomycota ~ Flare, data=f.phyla_RA_meta.v2, method="wilcox.test",p.adjust.method = "bonferroni")
# .y.        group1    group2        p p.adj p.format p.signif method  
# <chr>      <chr>     <chr>     <dbl> <dbl> <chr>    <chr>    <chr>   
#   1 Ascomycota HHP       Remission 0.851  1    0.85     ns       Wilcoxon
# 2 Ascomycota HHP       Flare     0.277  0.83 0.28     ns       Wilcoxon
# 3 Ascomycota Remission Flare     0.187  0.56 0.19     ns       Wilcoxon

compare_means(Ascomycota ~ Flare, data=f.phyla_RA_meta.v2, method="kruskal.test",p.adjust.method = "bonferroni")
# .y.            p p.adj p.format p.signif method        
# <chr>      <dbl> <dbl> <chr>    <chr>    <chr>         
#   1 Ascomycota 0.372  0.37 0.37     ns       Kruskal-Wallis

# next, Basidiomycota
# Kruskal-Wallis test is an ANOVA for non-normal data
fit2<-kruskal.test(Basidiomycota ~ Flare, data=f.phyla_RA_meta.v2)
fit2

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(f.phyla_RA_meta.v2, Basidiomycota ~ Flare, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(Basidiomycota ~ Flare, data = f.phyla_RA_meta.v2)
# Fligner-Killeen:med chi-squared = 1.6657, df = 2, p-value = 0.4348
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(Basidiomycota ~ Flare, data=f.phyla_RA_meta.v2, method="wilcox.test",p.adjust.method = "bonferroni")
compare_means(Basidiomycota ~ Flare, data=f.phyla_RA_meta.v2, method="kruskal.test",p.adjust.method = "bonferroni")

#### Visualize Phylum by Group ####

mycomp <- list(c("Remission", "Flare"), c("HHP", "Flare"), c("HHP", "Remission"))

phy.bxplt1<-ggboxplot(
  f.phyla_melt_meta.v2, x="Flare", y="Count", 
  add = c("jitter"), 
  facet.by = "Phylum", nrow = 1,
  fill = "Flare", palette = c("#007561", "#ff6361", "#ffa600"))+
  stat_compare_means(aes(group = Flare),label.y = 1.4) +
  stat_compare_means(comparisons = mycomp)+
  scale_y_continuous(name = "Relative Abundance", labels = scales::percent)+
  rremove("x.ticks")+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(colour = "black", size=12),
        axis.text=element_text(colour="black", size=12),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))

ggsave(phy.bxplt1,
  filename="figures/RelativeAbundance/Phylum/ABS_ITS2_FungalPhyla_variance_by_flare_boxplots.png",
  width = 5000,
  height = 2200,
  units = "px",
  dpi = 200,
) 

# drop fungal phyla that do not significantly change between the groups (i.e., No Wilcoxon tests reported in last fig)
exclude.phy<-c("Basidiobolomycota", "Aphelidiomycota", "Monoblepharomycota")
f.phyla_melt_meta.v3<-f.phyla_melt_meta.v2[!grepl(paste(exclude.phy,collapse="|"),f.phyla_melt_meta.v2$Phylum),]

phy.bxplt2<-ggboxplot(
  f.phyla_melt_meta.v3, x="Flare", y="Count", 
  add = c("jitter"), 
  facet.by = "Phylum", nrow = 1,
  fill = "Flare", palette = c("#007561", "#ff6361", "#ffa600"))+
  stat_compare_means(aes(group = Flare),label.y = 1.4) +
  stat_compare_means(comparisons = mycomp)+
  scale_y_continuous(name = "Relative Abundance", labels = scales::percent)+
  rremove("x.ticks")+
  theme(axis.title.x = element_blank(),
        legend.position = "none",
        text = element_text(colour = "black", size=12),
        axis.text=element_text(colour="black", size=12),
        axis.text.x = element_text(angle=45,vjust=1,hjust=1))

ggsave(phy.bxplt2,
       filename="figures/RelativeAbundance/Phylum/ABS_ITS2_FungalPhyla_variance_by_flare_NoLowPhyla_boxplots.png",
       width = 4000,
       height = 2200,
       units = "px",
       dpi = 200,
) 

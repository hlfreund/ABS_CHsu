#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
#setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaDust")
suppressPackageStartupMessages({ # load packages quietly
  library(devtools)
  library(phyloseq)
  library(ggplot2)
  library(vegan)
  library(ggpubr)
  library(lme4)
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
  library(nationalparkcolors)
  library(shades)
  library(ALDEx2)
  library(rstatix)
  library(decontam)
  library(ggvegan)
  library(microbiome)
  library(pairwiseAdonis)
  library(corrplot)
  library(lmtest)
  #library(car) # car is imported w/ rstatix so you do not need to load library individually?
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file
#load("data/SSeaDust_AlphaDiv_Data_Rarefied.Rdata") # includes Shannon entropy, Shannon div, and species richness of rarefied counts

#save.image("data/Env_Seqs_All/env.seq_analysis.Rdata") # save global env to Rdata file
b.dust.all[1:6,1:6]
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(meta.all.scaled)

## DO NOT RUN THIS LINE, THIS IS YOUR COLOR REFERENCE!!!!
#(WI="#ef781c",DP="#03045e",BDC="#059c3f")

#### Rarefaction Curves & Species Accumulation Curves ####
# bacteria/archaea

# Species Accumulation Curve
sc2<-specaccum(bac.ASV_round.table[,-1],"random")
plot(sc2, ci.type="poly", col="darkgreen", lwd=2, ci.lty=0, ci.col="lightgreen")
boxplot(sc2, col="yellow", add=TRUE, pch=20)

# Prep for Rarefaction Curve
rowSums(bac.ASV_round.table[,-1]) # total # ASVs per sample, excluding SampleID from calculation
sort(colSums(bac.ASV_round.table[,-1])) # counts per ASV
sort(rowSums(bac.ASV_round.table[,-1])) # ASVs per sample

# Create Rarefaction curve
png('figures/AlphaDiversity/SSD_16S_rarecurve.png')
rarecurve(as.matrix(bac.ASV_round.table[,-1]),col=meta.all.scaled$SampDate_Color, step=1000, label=F,ylab="ASVs")
# to show sampel labels per curve, change label=T
dev.off()

#### Good's Coverage ####

# Good's coverage is defined as 1 - (F1/N) where F1 is the number of singleton OTUs and N is the sum of counts for all OTUs.
# If a sample has a Good's coverage == . 98, this means that 2% of reads in that sample are from OTUs that appear only once in that sample.
# NOTE: Good's Coverage is not the best estimate of coverage because DADA2 removes singletons during process - so we are basing this on singeltons that remain (??)
bac.ASV.melt<-melt(bac.ASV_round.table,by="SampleID")
colnames(bac.ASV.melt)[which(colnames(bac.ASV.melt) == "variable")] <- "ASV_ID"
colnames(bac.ASV.melt)[which(colnames(bac.ASV.melt) == "value")] <- "Count"

# calculate # of singletons per sample
sings<-data.frame(Singletons=rowSums(bac.ASV_round.table[,-1]==1),SampleID=bac.ASV_round.table$SampleID) # how many ASVs have a count of only 1
# no singletons?

# find total counts per sample
totseq<-data.frame(TotalSeqs=rowSums(bac.ASV_round.table[,-1]),SampleID=bac.ASV_round.table$SampleID) # total # of counts per sample

# Calculate Good's Coverage
good.cov<-data.frame(GoodsCov=(100*(1-(sings$Singletons/totseq$TotalSeqs))),SampleID=sings$SampleID)

# Merge Total Seqs and Goods DFs together for plotting
good.df<-merge(good.cov,totseq,by="SampleID")

ggplot(good.df, aes(x=TotalSeqs,y=GoodsCov))+geom_point()+xlab("Total Seqs per Sample")+ylab("Good's Coverage (%)")

#### Repeated Rarefaction of Raw Counts & Averaging Shannon Diversity ####
# this section explores different ways of utilizing rarefaction and rarefaction functions in vegan

# in vegan ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs)
min.rar<-min(rowSums(bac.ASV_round.table[,-1])) ## seeing min sum of OTUs so we can see what min is for rarefaction
min.rar

bac.ASV.rar1<-rarefy(bac.ASV_round.table[,-1],min.rar) ## be cognizant of min for rarefaction
bac.ASV.rar1 # rarefy gives the expected species richness in random subsamples of size sample from the community

bac.ASV.rar<-rrarefy(bac.ASV_round.table[,-1],min.rar)
# rrarefy generates one randomly rarefied community data frame or vector of given sample size
# The random rarefaction is made without replacement so that the variance of rarefied communities is rather related to rarefaction proportion than to the size of the sample
rowSums(bac.ASV.rar) # should all be the same after rarefaction

# Calculating Average Shannon Diversity w/ repeated rarefaction with rrarefy()
sdiv.rarefy<-function(df,min){
  # rownames must be samples, ASVs must by columns
  # min = min(rowSums(ASV table)) aka minimum sum of ASVs in each sample
  rare<-rrarefy(df,min)
  # rrarefy generates one randomly rarefied community data frame or vector of given sample size
  # The random rarefaction is made without replacement so that the variance of rarefied communities is rather related to rarefaction proportion than to the size of the sample

  s.ent<-vegan::diversity(rare, index="shannon")
  # calculate Shannon entropy from rarefied data

  sdiv<- exp(s.ent) # Shannon Diversity aka Hill number 1
  # calculate Shannon diversity from rarefied data

  return(sdiv)
}

# vv contains average Shannon Diversity calculations after rarefaction and Shan div calculations 100 times
ave.sdiv <- data.frame(AveShanDiv=rowMeans(data.frame(lapply(as.list(1:100), function(x) sdiv.rarefy(bac.ASV_round.table[,-1], min.rar)))))
# sdiv.rarefy function --> Rrarefy ASV warnings() table with given min, calculate Shannon diversity
# lapply(as.list(1:100), function(x) sdiv.rarefy(bac.ASV_round.table[,-1], min.rar))) --> creates list by running sdiv.rarefy() function 100 times & storing results in list each time
# data.frame(lapply(...)) --> saves lapply() output as data frame, not list
# rowMeans(df(lapply(...))) --> get rowMeans aka sample means (aka the average) of Shannon diversity after 100 calculations of Shan div.

# compare our results with just one iteration of Shannon Div calculation
## calculate Shannon entropy from rarefied data
bac.S.Ent<-vegan::diversity(bac.ASV.rar, index="shannon")

## calculations Shannon diversity from entropy
bac.S.Div<- exp(bac.S.Ent) # Shannon Diversity aka Hill number 1

# how are our results verses the vegan::diversity() function's results?
as.data.frame(bac.S.Div)
ave.sdiv
# Results from bac.S.Div are similar to ave.sdiv results which is great!

# do the same thing with Shannon entropy just so you have those data
# calculate the average Shannon entropy of rarefied data (taking average of 100 calculations of Shannon entropy)
ave.s.ent<-data.frame(AveShanEnt=rowMeans(data.frame(lapply(as.list(1:100), function(x) vegan::diversity(bac.ASV.rar, index="shannon")))))

bac.ASV.probs<-drarefy(bac.ASV_round.table[,-1],min.rar)
# drarefy returns probabilities that species occur in a rarefied community of size sample


#### Rarefaction of Raw Counts & Averaging Species Richness ####
# this section explores different ways of utilizing rarefaction and rarefaction functions in vegan

# in vegan ROWS need to be SITES/samples; COLUMNS are SPECIES (OTUs, ASVs)
min.rar<-min(rowSums(bac.ASV_round.table[,-1])) ## seeing min sum of OTUs so we can see what min is for rarefaction
min.rar

bac.ASV.rar1<-rarefy(bac.ASV_round.table[,-1],min.rar) ## be cognizant of min for rarefaction
bac.ASV.rar1 # rarefy gives the expected species richness in random subsamples of size sample from the community

bac.ASV.rar<-rrarefy(bac.ASV_round.table[,-1],min.rar)
# rrarefy generates one randomly rarefied community data frame or vector of given sample size
# The random rarefaction is made without replacement so that the variance of rarefied communities is rather related to rarefaction proportion than to the size of the sample
rowSums(bac.ASV.rar) # should all be the same after rarefaction

# Calculating Average Shannon Diversity w/ repeated rarefaction with rrarefy()
spech.rich.rarefy<-function(df,min){
  # rownames must be samples, ASVs must by columns
  # min = min(rowSums(ASV table)) aka minimum sum of ASVs in each sample
  rare<-rrarefy(df,min)
  # rrarefy generates one randomly rarefied community data frame or vector of given sample size
  # The random rarefaction is made without replacement so that the variance of rarefied communities is rather related to rarefaction proportion than to the size of the sample

  s.rich<-specnumber(rare)

  return(s.rich)
}

# vv contains average Shannon Diversity calculations after rarefaction and Shan div calculations 100 times
ave.spec.rich <- data.frame(AveSpecRich=rowMeans(data.frame(lapply(as.list(1:100), function(x) spech.rich.rarefy(bac.ASV_round.table[,-1], min.rar)))))
# spec.rich.rarefy function --> Rrarefy ASV warnings() table with given min, calculate Shannon diversity
# lapply(as.list(1:100), function(x) sdiv.rarefy(bac.ASV_round.table[,-1], min.rar))) --> creates list by running spec.rich.rarefy() function 100 times & storing results in list each time
# data.frame(lapply(...)) --> saves lapply() output as data frame, not list
# rowMeans(df(lapply(...))) --> get rowMeans aka sample means (aka the average) of species richness after 100 calculations of species richness

ave.spec.rich

# Calculate species richness (number of species per sample) from rarefied data to compare
specnumber(bac.ASV.rar)
specnumber(bac.ASV_round.table[,-1])

# results are similar!

bac.ASV.probs<-drarefy(bac.ASV_round.table[,-1],min.rar)
# drarefy returns probabilities that species occur in a rarefied community of size sample


# #### Variance Stabilizing Transformation VST of Raw (?) counts ####
# #Prepare Contig Feature Count Data for Normalization w/ DESeq2
# # make sure count data & mgm_meta are in the same order
# bac.ASV_round.table[1:5,1:5]
#
# #bac.ASV_matrix<-as.matrix(bac.ASV_round.table[,!names(bac.ASV_round.table) %in% c("SampleID")]) # convert count table into matrix & exclude column called SampleID
# bac.ASV_matrix2<-t(as.matrix(bac.ASV_round.table[,-1])) # transpose matrix so ASVs are rows, samples are columns
# # ^ will be used in DESeq2 functions
# rownames(meta.all.scaled) %in% colnames(bac.ASV_matrix2) # check if rownames in meta.all.scaled (SampleID) match column names in count data
# dim(meta.all.scaled)
# dim(bac.ASV_matrix2)
#
# # Reorder matrix to match order of meta.all.scaled
# bac.ASV_matrix2=bac.ASV_matrix2[,rownames(meta.all.scaled)] ## reorder ASV matrix by column name to match order of rownames in meta.all.scaled
# colnames(bac.ASV_matrix2) # sanity check that this reordering worked
# rownames(meta.all.scaled)
#
# # create the DESeq DataSet object for DGE analysis
# # DESeq2 needs whole # data, so need raw read counts, NOT coverage for these tables...questionable
# # b_dds has the scaled coverage that was calculated by dividing reads from featureCounts by gene lengths
# b_dds<-DESeqDataSetFromMatrix(countData=round(bac.ASV_matrix2), colData = meta.all.scaled, design = ~ 1)
#
# # design = ~ 1 means no design
# head(counts(b_dds))
# colSums(counts(b_dds)) %>% barplot
#
# # Estimate size factor - aka normalization factor, divide all read counts by each size factor per sample
# b_dds <- estimateSizeFactors(b_dds,type="ratio")
# ## To calculate size factor in DESeq2:
# # calculates geometric mean of each gene in each sample and uses this as a pseudoreference
# # calculates ratio of each sample by dividing each gene count by it's pseudoreference in each sample
# # The median value of all ratios for a given sample is taken as the normalization factor (size factor)
# # The differentially expressed genes should not affect the median value
# # median of ratios method makes the assumption that not ALL genes are differentially expressed; therefore, the normalization factors should account for sequencing depth and RNA composition of the sample
# ## (large outlier genes will not represent the median ratio values)
# # more here: https://hbctraining.github.io/DGE_workshop/lessons/02_DGE_count_normalization.html
#
# sizeFactors(b_dds)
#
# plot(sizeFactors(b_dds), colSums(counts(b_dds)))
#
# ### Variance Stabilizing Transformation
# # you should be able to use matrix or DESeq2 object for this next function, but matrix was not working?
# b_vst1 <- varianceStabilizingTransformation(b_dds) # add pseudocount
# assay(b_vst1) #see output of VST
#
# b.vst<-assay(b_vst1)
#
#### Compare Transformed ASVs vs Raw Counts ####

bac.ASV_round.table$SampleID = factor(bac.ASV_round.table$SampleID, levels=unique(bac.ASV_round.table$SampleID[order(rowSums(bac.ASV_round.table[,-1]))]), ordered=TRUE)

raw.tot.counts<-ggplot(data=bac.ASV_round.table, aes(x=bac.ASV_round.table$SampleID, y=rowSums(bac.ASV_round.table[,-1]))) +
  geom_bar(stat="identity",colour="black",fill="dodgerblue")+theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="Total ASVs per Sample",subtitle="Based on Raw ASV Counts")+ylab("Total ASVs")+xlab("SampleID")

ggsave(raw.tot.counts,filename = "figures/AlphaDiversity/SSD_16S_Total_ASVs_per_Sample_barplot.png", width=13, height=10, dpi=600,create.dir = TRUE)

# Calculate rarefied ASVs per Sample

ggplot(data=as.data.frame(bac.ASV.rar), aes(x=rownames(bac.ASV.rar), y=rowSums(bac.ASV.rar))) +
  geom_bar(stat="identity",colour="black",fill="dodgerblue")+theme(axis.text.x = element_text(angle = 45, hjust=1)) +
  labs(title="Total Rarefied ASVs per Sample",subtitle="Based on Rarefied ASV Counts")+ylab("Total Rarefied ASVs")+xlab("SampleID")

# # Calculate VST ASVs per Sample
# total_vst_asvs<-data.frame(ASV_Total=colSums(b.vst),meta.all.scaled)
# total_vst_asvs$SampleID = factor(total_vst_asvs$SampleID, levels=unique(total_vst_asvs$SampleID[order(total_vst_asvs$ASV_Total)]), ordered=TRUE)
#
# ggplot(data=total_vst_asvs, aes(x=SampleID, y=ASV_Total,fill=Sample_Type)) +
#   geom_bar(stat="identity",colour="black")+theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

## average ASV per sample month & depth
# aggregate(bac.ASV_all$Count, list(bac.ASV_all$SampleMonth), FUN=mean)
# aggregate(bac.ASV_all$Count, list(bac.ASV_all$Depth_m), FUN=mean)

#### Alpha Diversity & Species Richness - Rarefied Data ####

# for code on how average Shannon entropy & Shannon diversity were calculated, please see the Rarefaction of Raw Counts section
head(ave.s.ent)
head(ave.sdiv)

# create data frame with Shannon entropy and Shannon diversity values
div_16s.rar<-data.frame(Bac_Shannon_Entropy=ave.s.ent,AveShanDiv=ave.sdiv)
class(div_16s.rar)
head(div_16s.rar)

div_16s.rar$SampleID<-rownames(div_16s.rar)
head(div_16s.rar)

ave.spec.rich$SampleID<-rownames(ave.spec.rich)

# merge richness and diversity dataframes together
d.r_16s.rar<-merge(div_16s.rar, ave.spec.rich, by.x="SampleID", by.y="SampleID")

# merge w/ meta.all.scaled
bac.div.metadat.rar <- merge(d.r_16s.rar,meta.all.scaled, by.x="SampleID", by.y="SampleID")
head(bac.div.metadat.rar)
class(bac.div.metadat.rar) # want data frame

unique(bac.div.metadat.rar$SampDate) # see how many elements there are in the Group variable
unique(bac.div.metadat.rar$Site) # see how many elements there are in the Group variable

# drop the outliers
#bac.div.metadat.rar<-bac.div.metadat.rar[bac.div.metadat.rar$AveShanDiv<300 & bac.div.metadat.rar$AveSpecRich>100,]

# create numeric variable for depth to be used for models later

# Find highest/lowest values of Shannon div per sample date
max(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="WI"]) # max div WI
min(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="WI"]) # min div WI
bac.div.metadat.rar[bac.div.metadat.rar$Site=="WI",]

max(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="DP"]) # max div DP
min(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="DP"]) # min div DP
bac.div.metadat.rar[bac.div.metadat.rar$Site=="DP",]

max(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="BDC"]) # max div BDC
min(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="BDC"]) # min div BDC
bac.div.metadat.rar[bac.div.metadat.rar$Site=="BDC",]

max(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="PD"]) # max div BDC
min(bac.div.metadat.rar$AveShanDiv[bac.div.metadat.rar$Site=="PD"]) # min div BDC
bac.div.metadat.rar[bac.div.metadat.rar$Site=="PD",]

# Find highest/lowest values of Species richness per sample date
max(bac.div.metadat.rar$AveSpecRich[bac.div.metadat.rar$Site=="WI"]) # max sr WI
min(bac.div.metadat.rar$AveSpecRich[bac.div.metadat.rar$Site=="WI"]) # min sr WI
bac.div.metadat.rar[bac.div.metadat.rar$Site=="WI",]

max(bac.div.metadat.rar$AveSpecRich[bac.div.metadat.rar$Site=="DP"]) # max sr DP
min(bac.div.metadat.rar$AveSpecRich[bac.div.metadat.rar$Site=="DP"]) # min sr DP
bac.div.metadat.rar[bac.div.metadat.rar$Site=="DP",]

max(bac.div.metadat.rar$AveSpecRich[bac.div.metadat.rar$Site=="BDC"]) # max sr BDC
min(bac.div.metadat.rar$AveSpecRich[bac.div.metadat.rar$Site=="BDC"]) # min sr BDC
bac.div.metadat.rar[bac.div.metadat.rar$Site=="BDC",]

max(bac.div.metadat.rar$AveSpecRich[bac.div.metadat.rar$Site=="PD"]) # max sr BDC
min(bac.div.metadat.rar$AveSpecRich[bac.div.metadat.rar$Site=="PD"]) # min sr BDC
bac.div.metadat.rar[bac.div.metadat.rar$Site=="PD",]

# save diversity data
save.image("data/SSeaDust_AlphaDiv_Data_Rarefied.Rdata")

#### Using Shapiro-Wilk test for Normality - Rarefied Data ####
shapiro.test(bac.div.metadat.rar$AveShanDiv) # what is the p-value?
# W = 0.80807, p-value = 0.0001486
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat.rar$AveShanDiv, col="blue") # with outliars

# visualize Q-Q plot for alpha div
# The Q-Q plot, or quantile-quantile plot, is a graphical tool to help us assess if a set of data plausibly came from some theoretical distribution such as a normal or exponential.
# For example, if we run a statistical analysis that assumes our residuals are normally distributed, we can use a normal Q-Q plot to check that assumption
# more on Q-Q plots here: https://data.library.virginia.edu/understanding-q-q-plots/
# more here too: https://grodri.github.io/glms/notes/c2s9#:~:text=8%20The%20Q%2DQ%20Plot,versus%20quantiles%20of%20a%20distribution.
qqnorm(bac.div.metadat.rar$AveShanDiv, pch = 1, frame = FALSE)
qqline(bac.div.metadat.rar$AveShanDiv, col = "red", lwd = 2)

shapiro.test(bac.div.metadat.rar$AveSpecRich) # what is the p-value?
# W = 0.81983, p-value = 0.0002437
# p > 0.05 states distribution of data are not significantly different from normal distribution
# p < 0.05 means that data is significantly different from a normal distribution
hist(bac.div.metadat.rar$AveSpecRich, col="blue")

# visualize Q-Q plot for species richness
qqnorm(bac.div.metadat.rar$AveSpecRich, pch = 1, frame = FALSE) # with outliars
qqline(bac.div.metadat.rar$AveSpecRich, col = "red", lwd = 2)

# #### Compare Means of Shannon Diversity - Rarefied Data ####
# head(bac.div.metadat.rar)
#
# t.test(AveShanDiv ~ Site, data=bac.div.metadat.rar)
# site_list<-unique(bac.div.metadat.rar$Site)
#
# site_div_subsets<-lapply(site_list, function(x) {subset(bac.div.metadat.rar, Site==x)}) # create list of separated metadata by Site
#
# # set up the function and run this to store it in our Global environment
# df_specific.subset<-function(var_vec,var_subsets){
#   # var_vec = vector of variable elements from specific categorical variable;
#   ## e.g. vector of names from Site categorical variable (metadata sites)
#   # var_subsets = list of dataframes subsetted by column$element from original dataframe;
#   ## e.g. list of dataframes (each df = element of list) subsetted from metadata using vector of metadata$Site names
#   for(i in seq_along(var_vec)){
#     # print(var_vec[i]) -- var_vec[i] = each element in var_vec
#     # print(var_subsets[[i]]) -- var_subsets[[i]] = each sub
#     df<-paste(var_vec[i])
#     #print(df)
#     assign(df, var_subsets[[i]], envir = .GlobalEnv)
#     print(paste("Dataframe", var_vec[i] ,"done"))
#
#   }
#
# } # create function that takes each subset in the list and makes a dataframe out of it saved in our Global Env
#
# # run the function
# df_specific.subset(site_list, site_div_subsets) # used metadata + rarefied div values, species richness values attached
#
# multi.t.test.fxn<-function(df1,df2){
#   # create empty lists to store stuff & model number (mdlnum) to keep track of models each iteration of loop in fxn
#   t.tests<- vector('list', ncol(df1) * ncol(df2)) # create empty list where the GLM output is stored
#   sig.ttest.results<-vector('list', ncol(df1) * ncol(df2))
#   mdlnum <- 1 # counting our model numbers for indexes purposes in the loop
#
#   # run the nested loop that generates GLMs from each data frame
#   ## df1[i] is dependent variable (y), df2[j] is independent variable (x) in GLM
#   for (i in 1:ncol(df1)){ # for each column in df1
#     for (j in 1:ncol(df2)){ # for each column in df2
#       t.tests[[mdlnum]] <-t.test(df1[,i],df2[,j]) # run the GLM with the gaussian distribution, where df1[i] is your dependent variable and df2[j] is your independent variable
#       #results_[[mdlnum]] <-summary(t.tests[[mdlnum]]) # save results of glm into list called results
#       names(t.tests)[mdlnum]<-paste(names(df1)[i],",",names(df2)[j]) # rename list element to contain the name of the columns used in the model
#
#       # save only significant GLMs to another list called sig.ttest.results
#       ## if p-value < 0.05, save to sig.ttest.results list
#       ifelse(t.tests[[mdlnum]]$p.value < 0.05, sig.ttest.results[[mdlnum]]<-t.tests[[mdlnum]], "Not Sig")
#       # structure of ifelse(): ifelse(the test or condition, what to do if it satisfies condition, otherwise if condition isn't satisfied do this)
#       names(sig.ttest.results)[mdlnum]<-paste(names(df1)[i],",",names(df2)[j])
#       mdlnum <- mdlnum + 1 # add 1 to modelnumber so we keep track of # of models (for indexing purposes in list)
#
#     }
#   }
#
#   # drop all NULL elements from sig.ttest.results list so it only includes significant GLMs
#   sig.ttest.results[sapply(sig.ttest.results, is.null)] <- NULL
#
#   # assign lists to global env so they are saved there are function ends
#   assign("results.t.tests", t.tests,envir = .GlobalEnv)
#   assign("sig.ttest.results.t.tests", sig.ttest.results,envir = .GlobalEnv)
#
# }
#
# # First run individual t-tests by time point
# # WI vs DP
# t.test.WI.DP<-t.test(WI.divmeta$AveShanDiv,DP.divmeta$AveShanDiv)
# t.test.WI.DP$p.value
#
# # WI vs BDC
# t.test.WI.BDC<-t.test(WI.divmeta$AveShanDiv,BDC.divmeta$AveShanDiv)
# t.test.WI.BDC$p.value
#
# # WI vs PD
# t.test.WI.PD<-t.test(WI.divmeta$AveShanDiv,PD.divmeta$AveShanDiv)
# t.test.WI.PD$p.value
#
# # DP vs BDC
# t.test.DP.BDC<-t.test(DP.divmeta$AveShanDiv,BDC.divmeta$AveShanDiv)
# t.test.DP.BDC$p.value
#
# # Combine the p-values
# ttest.Div.pvals<-c(t.test.WI.DP$p.value, t.test.WI.BDC$p.value, t.test.DP.BDC$p.value)
#
# # Adjust the p-values based on the # of comparisons you did
# p.adjust(ttest.Div.pvals, method="bonferroni",n=3)
# # ^ matches findings on the figure, which uses the t_test function from the rstatix package (see geom_pwc())
#
# #### Compare Means of Species Richness - Rarefied Data ####
# head(bac.div.metadat.rar)
# #WI.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$SampDate=="WI",]
# #DP.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$SampDate=="DP",]
# #BDC.divmeta<-bac.div.metadat.rar[bac.div.metadat.rar$SampDate=="BDC",]
#
# # First run individual wilcox-tests by time point
# # WI vs DP
# w.test.WI.DP<-wilcox.test(WI.divmeta$AveSpecRich,DP.divmeta$AveSpecRich)
# w.test.WI.DP$p.value
#
# # WI vs BDC
# w.test.WI.BDC<-wilcox.test(WI.divmeta$AveSpecRich,BDC.divmeta$AveSpecRich)
# w.test.WI.BDC$p.value
#
# # DP vs BDC
# w.test.DP.BDC<-wilcox.test(DP.divmeta$AveSpecRich,BDC.divmeta$AveSpecRich)
# w.test.DP.BDC$p.value
#
# # Combine the p-values
# wilc.SR.pvals<-c(w.test.WI.DP$p.value, w.test.WI.BDC$p.value, w.test.DP.BDC$p.value)
#
# # Adjust the p-values based on the # of comparisons you did
# p.adjust(wilc.SR.pvals, method="bonferroni",n=3)
# # ^ matches findings on the figure, which uses the wilcox_test function from the rstatix package (see geom_pwc())
#
#### Visualize Alpha Diversity & Species Richness - from Rarefied Data ####

## Shannon Diversity by Site & Sample Date
bac.a.div.rar<-ggplot(bac.div.metadat.rar, aes(x=Site, y=AveShanDiv)) +geom_jitter(aes(color=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Site",values=unique(bac.div.metadat.rar$Site_Color[order(bac.div.metadat.rar$Site)]),labels=c("PD","BDC","DP","WI")) +
  geom_boxplot(fill=NA, outlier.color=NA)+scale_x_discrete(labels=c("PD","BDC","DP","WI"))+theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Shannon Diversity by Site", subtitle="Using Rarefied Counts", x="Site", y="Shannon Diversity", color="Depth (m)")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
  #geom_pwc(method = "t_test", label = "p.adj.format",p.adjust.method = "bonferroni")

ggsave(bac.a.div.rar,filename = "figures/AlphaDiversity/RarefiedCounts/SSD_16S_rarefied_alpha_diversity_site_boxplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

bac.a.div.rarB<-ggplot(bac.div.metadat.rar, aes(x=Site, y=AveShanDiv)) +geom_jitter(aes(color=SampDate,shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(bac.div.metadat.rar$SampDate_Color[order(bac.div.metadat.rar$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA)+theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Shannon Diversity by Collection Date & Site", subtitle="Using Rarefied Counts", x="Site", y="Shannon Diversity", color="Collection Date")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_shape_manual(values = c(0,1,16,15),labels=unique(bac.div.metadat.rar$Site[order(bac.div.metadat.rar$Site)])) + scale_x_discrete(labels=unique(bac.div.metadat.rar$Site[order(bac.div.metadat.rar$Site)]))
ggsave(bac.a.div.rarB,filename = "figures/AlphaDiversity/RarefiedCounts/SSD_16S_rarefied_alpha_diversity_sampledate_site_boxplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

bac.a.div.rarC<-ggplot(bac.div.metadat.rar, aes(x=SampDate, y=AveShanDiv)) +geom_jitter(aes(color=SampDate,shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(bac.div.metadat.rar$SampDate_Color[order(bac.div.metadat.rar$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Shannon Diversity by Collection Date & Site", subtitle="Using Rarefied Counts", x="Site", y="Shannon Diversity", color="Collection Date")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,size=10,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_shape_manual(values = c(0,1,16,15),labels=unique(bac.div.metadat.rar$Site[order(bac.div.metadat.rar$Site)]))

ggsave(bac.a.div.rarC,filename = "figures/AlphaDiversity/RarefiedCounts/SSD_16S_rarefied_alpha_diversity_sampledate_boxplot.png", width=15, height=10, dpi=600,create.dir = TRUE)


## Species Richness by Site & Sample Date
bac.a.sr.rar<-ggplot(bac.div.metadat.rar, aes(x=Site, y=AveSpecRich)) +geom_jitter(aes(color=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Site",values=unique(bac.div.metadat.rar$Site_Color[order(bac.div.metadat.rar$Site)]),labels=c("PD","BDC","DP","WI")) +
  geom_boxplot(fill=NA, outlier.color=NA)+scale_x_discrete(labels=c("PD","BDC","DP","WI"))+theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Species Richness by Site", subtitle="Using Rarefied Counts", x="Site", y="Species Richness", color="Depth (m)")+theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))
#geom_pwc(method = "t_test", label = "p.adj.format",p.adjust.method = "bonferroni")

ggsave(bac.a.sr.rar,filename = "figures/SpeciesRichness/RarefiedCounts/SSD_16S_rarefied_species_richness_site_boxplot.png", width=13, height=10, dpi=600,create.dir = TRUE)

bac.a.sr.rarB<-ggplot(bac.div.metadat.rar, aes(x=Site, y=AveSpecRich)) +geom_jitter(aes(color=SampDate,shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(bac.div.metadat.rar$SampDate_Color[order(bac.div.metadat.rar$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  geom_boxplot(fill=NA, outlier.color=NA)+theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Species Richness by Collection Date & Site", subtitle="Using Rarefied Counts", x="Site", y="Species Richness", color="Collection Date")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_shape_manual(values = c(0,1,16,15),labels=unique(bac.div.metadat.rar$Site[order(bac.div.metadat.rar$Site)])) + scale_x_discrete(labels=unique(bac.div.metadat.rar$Site[order(bac.div.metadat.rar$Site)]))

ggsave(bac.a.sr.rarB,filename = "figures/SpeciesRichness/RarefiedCounts/SSD_16S_rarefied_species_richness_sampledate_site_boxplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

bac.a.sr.rarC<-ggplot(bac.div.metadat.rar, aes(x=SampDate, y=AveSpecRich)) +geom_jitter(aes(color=SampDate,shape=Site), size=3, width=0.15, height=0) +
  scale_color_manual(name ="Collection Date",values=unique(bac.div.metadat.rar$SampDate_Color[order(bac.div.metadat.rar$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +theme_bw()+theme_classic()+
  labs(title = "Dust Bacterial Species Richness by Collection Date & Site", subtitle="Using Rarefied Counts", x="Site", y="Species Richness", color="Collection Date")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,size=10,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  scale_shape_manual(values = c(0,1,16,15),labels=unique(bac.div.metadat.rar$Site[order(bac.div.metadat.rar$Site)]))

ggsave(bac.a.sr.rarC,filename = "figures/SpeciesRichness/RarefiedCounts/SSD_16S_rarefied_species_richness_sampledate_boxplot.png", width=15, height=10, dpi=600,create.dir = TRUE)

#### Compare Variance w/ Shannon Diversity ####
# shannon diversity is NOT normally distributed; used rarefied counts to calculate SR
# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

# Kruskal-Wallis test is an ANOVA for non-normal data
fit1<-kruskal.test(AveShanDiv ~ Site, data=bac.div.metadat.rar)

fit1

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(bac.div.metadat.rar, AveShanDiv ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(AveShanDiv ~ Site, data = bac.div.metadat.rar)
# Fligner-Killeen:med chi-squared = 0.48165, df = 3, p-value = 0.9229
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(AveShanDiv ~ Site, data=bac.div.metadat.rar, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(AveShanDiv ~ Site, data=bac.div.metadat.rar, method="kruskal.test",p.adjust.method = "bonferroni")

# Kruskal-Wallis test is an ANOVA for non-normal data
fit2<-kruskal.test(AveShanDiv ~ interaction(Site,CollectionYear), data=bac.div.metadat.rar)

fit2
# Kruskal-Wallis chi-squared = 3.4803, df = 7, p-value = 0.8373

# Kruskal-Wallis test is an ANOVA for non-normal data
fit2<-kruskal.test(AveShanDiv ~ interaction(Site,CollectionYear), data=bac.div.metadat.rar)

fit2
# Kruskal-Wallis chi-squared = 3.4803, df = 7, p-value = 0.8373

#### Compare Variance w/ Species Richness ####
# species richness is NOT normally distributed; used rarefied counts to calculate SR
# use the following statisitcal tests for variance comparisons
## Kruskal: are variances significantly different between groups
## Dunn test: which groups' variances are significant different from one another
## Fligner test: is variance homogenous aka equal across samples?

# Kruskal-Wallis test is an ANOVA for non-normal data
fit3<-kruskal.test(AveSpecRich ~ Site, data=bac.div.metadat.rar)

fit3
# Kruskal-Wallis chi-squared = 1.4757, df = 3, p-value = 0.6879

# Instead of using Tukey test, we can use Dunn's test to see which groups significantly vary if Kruskal-Wallis test is significant
# ANOVA + Tukey for normally distributed data, Kruskal-Wallis + Dunn's test for non-normal data
rstatix::dunn_test(bac.div.metadat.rar, AveSpecRich ~ Site, p.adjust.method = "bonferroni", detailed = TRUE)

## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
### Fligner's test is a Levene's test for data that are not normally distributed
### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
## more here: https://www.geeksforgeeks.org/fligner-killeen-test-in-r-programming/
fligner.test(AveSpecRich ~ Site, data = bac.div.metadat.rar)
# Fligner-Killeen:med chi-squared = 0.67288, df = 3, p-value = 0.8796
# Which shows that the data do NOT deviate significantly from homogeneity.

compare_means(AveSpecRich ~ Site, data=bac.div.metadat.rar, method="wilcox.test",p.adjust.method = "bonferroni")

compare_means(AveSpecRich ~ Site, data=bac.div.metadat.rar, method="kruskal.test",p.adjust.method = "bonferroni")

# Kruskal-Wallis test is an ANOVA for non-normal data
fit4<-kruskal.test(AveSpecRich ~ interaction(Site,CollectionYear), data=bac.div.metadat.rar)

fit4


#### Linear Models with Surface Type Frequencies & Diversity,Richness ####
head(bac.div.metadat.rar)
head(meta.all.scaled)
head(d.r_16s.rar)

# let's fit distributions to Shan Div and Species Richness to see which family distribution would be best
# code & notes from here: https://pages.vassar.edu/jutouchon/files/2018/08/Vassar_R.GLMs_2018.pdf
# first, ShanDiv
fit.div1<-fitdistr(bac.div.metadat.rar$AveShanDiv,"normal")
fit.div2<-fitdistr(bac.div.metadat.rar$AveShanDiv,"lognormal")
fit.div3<-fitdistr(round(bac.div.metadat.rar$AveShanDiv),"Poisson")
fit.div4<-fitdistr(round(bac.div.metadat.rar$AveShanDiv),"negative binomial")
AIC(fit.div1,fit.div2,fit.div3,fit.div4)
# fit.div4 has the lowest AIC value, so negative binomial is probably the way to go

# first, ShanDiv
fit.sr1<-fitdistr(bac.div.metadat.rar$AveSpecRich,"normal")
fit.sr2<-fitdistr(bac.div.metadat.rar$AveSpecRich,"lognormal")
fit.sr3<-fitdistr(round(bac.div.metadat.rar$AveSpecRich),"Poisson")
fit.sr4<-fitdistr(round(bac.div.metadat.rar$AveSpecRich),"negative binomial")
AIC(fit.sr1,fit.sr2,fit.sr3,fit.sr4)
# fit.sr4 has the lowest AIC value, so negative binomial is probably the way to go

# just look at everything at once in step-wise fashion
# first shannon div
# step1<-step(glm.nb(formula = round(AveShanDiv) ~ ., data=bac.div.metadat.rar[,c(3,9:10,13:15,38:47)]))
# #                 Estimate Std. Error t value Pr(>|t|)
# summary(step1)
# plot(step1)
#
# summary(glm.nb(formula = AveShanDiv ~ Shrub , data=bac.div.metadat.rar))
#
# # then species richness
# step2<-step(glm.nb(formula = round(AveSpecRich) ~ ., data=bac.div.metadat.rar[,c(4,9:10,13:15,38:47)]))
# #                 Estimate Std. Error t value Pr(>|t|)
# summary(step2)
# plot(step2)
#
# summary(glm.nb(formula = AveSpecRich ~ ave.wind_direction+ave.air_temp , data=bac.div.metadat.rar))
#
# # create dfs of only STF + climate data and only the pcoa axes of interest
# Clim_only<-meta.all.scaled[,c(4,7:8,11:12,37:46)]
#
# head(Clim_only)
#
# # view df of just diversity data
# head(d.r_16s.rar)
# rownames(d.r_16s.rar)<-d.r_16s.rar$SampleID
#
# dim(Clim_only) # confirming that both data frames have the same # of rows
# dim(d.r_16s.rar)
#
# rownames(Clim_only) # check rownames to see if they are in the same order in both data frames
# rownames(d.r_16s.rar)
#
# # reorder data frames so they are in the same order by row (SampleID)
# Clim_only=Clim_only[rownames(d.r_16s.rar),] ## reorder metadata to match order of CLR data
#
# rownames(Clim_only) # check rownames to see if they are in the same order in both data frames after reordering
# rownames(d.r_16s.rar)
#
# glmtest<-glm.nb(d.r_16s.rar$AveShanDiv~Clim_only$ave.air_temp)

#### Correlations ####

# quick visualization of comparing all the variables to each other
## cannot include all variables of interest because figure margins become too big
par(mar=c(1,1,1,1))
pairs(bac.div.metadat.rar[, c(3:4,8,11:12,14:16,41:50)]) # average shannon div first
# maybe Herbaceoous x Shannon Div?
pairs(bac.div.metadat.rar[, c(3:4,8,11:12,14:16,41:50)]) # species richness next

# Visualize with a corrplot
cor_div.sr.stf <- cor(bac.div.metadat.rar[, c(3:4,8,11:12,14:16,41:50)], method='pearson')
cor_div.sr.stf

corrplot.mixed(cor_div.sr.stf, tl.pos='lt', tl.cex=0.4, sig.level = 0.05, insig="label_sig", number.cex=0.8,
               diag='l',cl.ratio = 0.2, tl.srt = 45)


# now let's run the correlations
head(meta.all.scaled)
head(bac.div.metadat.rar)

# create dfs of only surface type freq data, only climate data, & only the alpha div & species richness respectively
STF_Clim_Only<-meta.all.scaled[,c(4,7:8,10:12,38:47)]
head(STF_Clim_Only)

div.rich<-bac.div.metadat.rar[,c(1,3:4)]
head(div.rich)
rownames(div.rich)<-div.rich$SampleID
head(div.rich)

dim(STF_Clim_Only) # confirming that both data frames have the same # of rows
dim(div.rich)

rownames(STF_Clim_Only) # check rownames to see if they are in the same order in both data frames
rownames(div.rich)

# reorder data frames so they are in the same order by row (SampleID)
STF_Clim_Only=STF_Clim_Only[rownames(div.rich),] ## reorder metadata to match order of CLR data

rownames(STF_Clim_Only) # check rownames to see if they are in the same order in both data frames after reordering
rownames(div.rich)

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
  assign("env.div.rich.results.corrs", results_,envir = .GlobalEnv)
  assign("env.div.rich.sig.results.corrs", sig.results,envir = .GlobalEnv)

}

multi.univar.corr.fxn(div.rich[,-1],STF_Clim_Only) # test the function!

# ggplot(bac.div.metadat.rar, aes(x=ave.wind_direction, y=AveShanDiv)) +geom_jitter(aes(color=SampDate,shape=Site), size=3, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Date",values=unique(bac.div.metadat.rar$SampDate_Color[order(bac.div.metadat.rar$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +theme_bw()+theme_classic()+
#   labs(title = "Dust Bacterial Shannon Diversity by Collection Date & Site", subtitle="Using Rarefied Counts & Scaled Environmental Data", x="Ave Wind Direction", y="Shannon Diversity", color="Collection Date")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   scale_shape_manual(values = c(0,1,16,15),labels=unique(bac.div.metadat.rar$Site[order(bac.div.metadat.rar$Site)]))
#
# ggplot(bac.div.metadat.rar, aes(x=ave.wind_direction, y=AveSpecRich)) +geom_jitter(aes(color=SampDate,shape=Site), size=3, width=0.15, height=0) +
#   scale_color_manual(name ="Collection Date",values=unique(bac.div.metadat.rar$SampDate_Color[order(bac.div.metadat.rar$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +theme_bw()+theme_classic()+
#   labs(title = "Dust Bacterial Shannon Diversity by Collection Date & Site", subtitle="Using Rarefied Counts & Scaled Environmental Data", x="Ave Wind Direction", y="Species Richness", color="Collection Date")+
#   theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1,hjust=1,size=10),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
#   scale_shape_manual(values = c(0,1,16,15),labels=unique(bac.div.metadat.rar$Site[order(bac.div.metadat.rar$Site)]))
#

#### Save Everything ####
save.image("data/SSeaDust_AlphaDiv_Data_Rarefied.Rdata")


## Old ANOVA exmaple below ####
# # shannon diversity is normally distributed; used rarefied counts to calculate ShanDiv
# # use the following statisitcal tests for variance comparisons
# ## ANOVA: are variances significantly different between groups
# ## Tukey test: which groups' variances are significant different from one another
# ## Levene's test: is variance homogenous aka equal across samples?
#
# head(bac.div.metadat.rar)
#
# fit1<-anova(AveShanDiv ~ Site, data=bac.div.metadat.rar)
# # ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
# #pairwise.adonis(bac.div.metadat.rar$AveShanDiv, bac.div.metadat.rar$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
#
# fit1
# # Kruskal-Wallis chi-squared = 1.4631, df = 3, p-value = 0.6908
#
# p.adjust(summary(fit1)[[1]][["Pr(>F)"]][1],method="bonferroni")
#
# # Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
# Tuk1<-TukeyHSD(fit1)
# Tuk1$Site
# #                             diff        lwr      upr      p adj
# #              diff        lwr       upr     p adj
# # BDC-PD -18.178964 -103.47268  67.11475 0.9347493
# # DP-PD    3.784434  -81.50928  89.07815 0.9993253
# # WI-PD  -25.564970 -110.85869  59.72875 0.8412084
# # DP-BDC  21.963398  -63.33032 107.25712 0.8920056
# # WI-BDC  -7.386006  -92.67972  77.90771 0.9950723
# # WI-DP  -29.349404 -114.64312  55.94431 0.7788528
#
# # fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=bac.div.metadat.rar)
# # summary(fit.0)
# # TukeyHSD(fit.0)
#
# # Levene's test - test for homogeneity of variance
# ## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
# ## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# ## t test assumes that variances are the same, so Levene's test needs to be non significant
# car::leveneTest(AveShanDiv ~ Site, data = bac.div.metadat.rar, center=mean)
# # Levene's Test for Homogeneity of Variance (center = mean)
# #       Df F value Pr(>F)
# # group  3  0.8558 0.4774
# #       24
# # ^ p value is > 0.05 so we cannot reject the Null hypothesis -- variances are equal
#
# ## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
# ### Fligner's test is a Levene's test for data that are not normally distributed
# ### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
# ## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# # fligner.test(AveShanDiv ~ SampDate, data = bac.div.metadat.rar)
# # Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# # Which shows that the data do not deviate significantly from homogeneity.
#
# compare_means(AveShanDiv ~ Site, data=bac.div.metadat.rar, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input
#
#
# fit2<-aov(AveShanDiv ~ Site*CollectionYear, data=bac.div.metadat.rar)
# # ANOVA is basically a regression but w/ categorical variables more info here https://www.statology.org/anova-vs-regression/
# #pairwise.adonis(bac.div.metadat.rar$AveShanDiv, bac.div.metadat.rar$SampDate, p.adjust.m='bonferroni') # shows us variation for each sample to see which ones are different
#
# summary(fit2)
# #Df           Sum Sq Mean Sq    F value   Pr(>F)
# #Site         3   4194    1398   0.418  0.742
# #Residuals   24  80303    3346
#
# p.adjust(summary(fit2)[[1]][["Pr(>F)"]][1],method="bonferroni")
#
# # Tukey test - tells us which groups are significantly different from each other (more here: https://www.r-bloggers.com/2013/06/anova-and-tukeys-test-on-r/)
# Tuk2<-TukeyHSD(fit2)
# Tuk2
#
# # fit.0<-aov(DustComplexity ~ as.factor(Elevation), data=bac.div.metadat.rar)
# # summary(fit.0)
# # TukeyHSD(fit.0)
#
# # Levene's test - test for homogeneity of variance
# ## Levene’s test is an inferential statistic used to check if the variances of a variable obtained for two or more groups are equal or not when data comes from a non-normal distribution
# ## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# ## t test assumes that variances are the same, so Levene's test needs to be non significant
# car::leveneTest(AveShanDiv ~ Site, data = bac.div.metadat.rar, center=mean)
# # Levene's Test for Homogeneity of Variance (center = mean)
# #       Df F value Pr(>F)
# # group  3  0.8558 0.4774
# #       24
# # ^ p value is > 0.05 so we cannot reject the Null hypothesis -- variances are equal
#
# ## The Fligner-Killeen test is a non-parametric test for homogeneity of group variances based on ranks. It is useful when the data are non-normally distributed or when problems related to outliers in the dataset cannot be resolved.
# ### Fligner's test is a Levene's test for data that are not normally distributed
# ### It is also one of the many tests for homogeneity of variances which is most robust against departures from normality.
# ## Null hypothesis: all populations variances are equal; Alt Hypothesis: at least 1 sample has different variance (aka variances are NOT equal across samples)
# # fligner.test(AveShanDiv ~ SampDate, data = bac.div.metadat.rar)
# # Fligner-Killeen:med chi-squared = 1.9504, df = 2, p-value = 0.3771
# # Which shows that the data do not deviate significantly from homogeneity.
#
# compare_means(AveShanDiv ~ Site, data=bac.div.metadat.rar, method="anova",p.adjust.method = "bonferroni") # won't take as.factor(Elevation) as input



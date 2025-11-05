#### Set WD & Load Libraries ####
getwd() # use setwd("path/to/files") if you are not in the right directory
setwd("/Volumes/HLF_SSD/Aronson_Lab_Data/Salton_Sea/SaltonSeaDust")
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
  library(ggvegan)
  library(microbiome)
})

#### Load Global Env to Import Count/ASV Tables ####
load("data/Amplicon/SSDust_16S.V3V4_W23_Data_Ready.Rdata") # save global env to Rdata file
#load("data/Amplicon/SSD_16S_CLR_EucDist_Ready.Rdata")

head(b.dust.all)
bac.ASV_table[1:4,1:4]
bac.ASV_table[(nrow(bac.ASV_table)-4):(nrow(bac.ASV_table)),(ncol(bac.ASV_table)-4):(ncol(bac.ASV_table))] # last 4 rows & cols
head(dust_meta)
head(meta.all.scaled)

#### Create Centered Log-Ratio Table from ASV table ####
bac.ASV_table[1:4,1:4]
b.clr<-decostand(bac.ASV_table[,-1],method = "clr", pseudocount = 1) #CLR transformation
b.clr[1:4,1:4]
# df must have rownames are SampleIDs, columns are ASV IDs for vegan functions below

rownames(b.clr) %in% rownames(meta.all.scaled) # check if metadata and micro comp data have matching rownames

#### Look at Surface Type Frequencies by Sample for Modeling Later ####

# melt down data so that we can make a stacked barplot
STF.melt<-melt(SurfTypFreq[,c(1,3:12,15:16)],by=c("Site","SampleID","STF_Date"))
colnames(STF.melt)[which(names(STF.melt) == "variable")] <- "SurfaceType"
colnames(STF.melt)[which(names(STF.melt) == "value")] <- "Frequency"
head(STF.melt)

# create palette
colorset9 = as.data.frame(t(data.frame("BarrenLand"="peachpuff2","CropLand"="gold1","Developed"="gray","Forest"="darkgreen","Herbaceous"="limegreen",
                   "Mexico"="red1","OpenWater"="mediumblue","Others"="black","SaltonSea"="darkturquoise","Shrub"="saddlebrown")))
colorset9$SurfaceType<-rownames(colorset9)
colnames(colorset9)[which(names(colorset9) == "V1")] <- "ST_Color"
colorset9

# merge color palette and STF data together
STF.melt<-merge(STF.melt, colorset9, by="SurfaceType")
head(STF.melt)
STF.melt$SampDate_Color <- as.character(STF.melt$ST_Color)
head(SurfTypFreq)

# create factors for organizing data in plot
STF.melt$Site<-factor(STF.melt$Site,levels=c("PD","BDC","DP","WI"))
unique(STF.melt$STF_Date)
STF.melt$SampleID = factor(STF.melt$SampleID, levels=unique(STF.melt$SampleID[order(STF.melt$Site,STF.melt$STF_Date)]), ordered=TRUE)

# plot time
just.stfs1<-ggplot(STF.melt, aes(x=SampleID, y=Frequency, fill=SurfaceType))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Surface Type Frequencies in Salton Sea Dust", x="SampleID", y="Frequency", subtitle="",fill="Surface Type")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name ="Surface Type",values=unique(STF.melt$ST_Color[order(STF.melt$SurfaceType)]),labels=c("Barren Land","Crop Land","Developed","Forest","Herbaceous","Mexico","Open Water","Others","Salton Sea","Shrub"))
ggsave(just.stfs1,filename = "figures/SurfaceTypeFrequencies/SSD_SurfaceTypeFrequencys_barplot.png", width=10, height=10, dpi=600)

just.stfs2<-ggplot(STF.melt, aes(x=SampleID, y=Frequency, fill=SurfaceType))+geom_bar(stat="identity",colour="black")+scale_x_discrete()+theme_classic()+
  labs(title = "Surface Type Frequencies in Salton Sea Dust", x="SampleID", y="Frequency", subtitle="",fill="Surface Type")+
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(hjust=1,angle=45),legend.title.align=0.5, legend.title = element_text(size=13),legend.text = element_text(size=11),plot.title = element_text(size=15))+
  guides(fill=guide_legend(ncol=1))+
  scale_fill_manual(name ="Surface Type",values=unique(STF.melt$ST_Color[order(STF.melt$SurfaceType)]),labels=c("Barren Land","Crop Land","Developed","Forest","Herbaceous","Mexico","Open Water","Others","Salton Sea","Shrub")) +
  facet_wrap(vars(Site), scales = "free")
ggsave(just.stfs2,filename = "figures/SurfaceTypeFrequencies/SSD_SurfaceTypeFrequencys_bySite_barplot.png", width=10, height=10, dpi=600)

#### Separate All Data by Time points ####
# create metadata df that will contain scaled chemical data
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

# matching data with user defined function -- here is the function, must run to store function in Global env
# note: metadata must have rownames for this to work!
match_dat<-function(compdata, subset_metadata){
  ### * comp data and metadata need to have row names - rownames should be Sample IDs
  subset_comp_data = pullrow<-(is.element(row.names(compdata), row.names(subset_metadata)))
  subset_comp_data=compdata[pullrow,]
  return(subset_comp_data)
}

# loop through list containing each site's metadata and use match_dat to pair with CLR data
for (i in seq_along(site_subsets)){
  print(site_subsets[[i]]) # shows what is in each element within list
  #print(names(site_subsets[i]))
  new.clr.df<-match_dat(b.clr,site_subsets[[i]])
  assign(paste0("b.clr_",names(site_subsets[i])), new.clr.df,envir = .GlobalEnv)
}
# names(site_subsets[i]) --> gives us name of each element in list

# did the function work the way we wanted it to?

b.clr_WI[1:4,1:4]
rownames(WI) %in% rownames(b.clr_WI) # hopefully all of the rownames match, aka will get output of TRUE

#### Check Count Data Relationship w/ Env Variables (w/ DCA) ####
## remember, CCA assumes that our species have a unimodal relationship with our variables.
### unimodal = one maximum, think upsidedown bellcurve or something
## RDA assumes a linear relationship
## check the assumption w/ DCA
# ^ more on DCA here: https://ordination.okstate.edu/DCA.htm

# ALL data
# add pseudocount so row sums are > 0
b.clr.pseudo<-b.clr+1
b.dca = decorana(b.clr.pseudo)

#plot(b.dca) # may take too long to load, do not run unless you have to
b.dca #DCA1 axis length = 0.47172; use RDA
## The length of first DCA axis:
## > 4 indicates heterogeneous dataset on which unimodal methods should be used (CCA),
##  < 3 indicates homogeneous dataset for which linear methods are suitable (RDA)
## between 3 and 4 both linear and unimodal methods are OK.

# BY MONTH

b.clr_WI.pseudo<-b.clr_WI+1
b.WI.dca = decorana(b.clr_WI.pseudo)
b.WI.dca #DCA1 axis length = 0.54543; use RDA

b.clr_DP.pseudo<-b.clr_DP+1
b.DP.dca = decorana(b.clr_DP.pseudo)
b.DP.dca #DCA1 axis length = 0.64485; use RDA

b.clr_BDC.pseudo<-b.clr_BDC+1
b.BDC.dca = decorana(b.clr_BDC.pseudo)
b.BDC.dca #DCA1 axis length = 0.46958; use RDA

b.clr_PD.pseudo<-b.clr_PD+1
b.PD.dca = decorana(b.clr_PD.pseudo)
b.PD.dca #DCA1 axis length = 0.52948; use RDA

#### RDA w/ All Data ####

rownames(meta.all.scaled) %in% rownames(b.clr) # check order of DFs
head(meta.all.scaled)

rda.all.0<-rda(b.clr ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.0
summary(rda.all.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.0) # -0.03243668 -- bad model!
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.0, permutations = how(nperm=999)) # p = 0.699

## we can also do a permutation test by RDA axis
#anova(rda.all.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.0)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 12.795263              8.980873              5.426993              6.155984             43.324469             49.271714
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others
# 12.288065             21.023371             14.459425             22.459255              7.222826             30.982534
# SaltonSea                 Shrub
# 19.061052                    NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.a = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,9:10,35:44)]),
                     scope=formula(rda.all.0),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ ave.wind direction, with Developed + Forest near sig

rda.all.a$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
# + ave.wind_direction  1 257.35 2.3999  0.002 **
#   + Developed           1 257.70 1.5215  0.042 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.a2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,9:10,35:44)]),
                        scope=formula(rda.all.0),
                        permutations = how(nperm=999))
rda.all.a2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# Dropping OpenWater and Forest due to largest p values in ordistep
rda.all.1<-rda(b.clr ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Mexico+Others+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.1
summary(rda.all.1)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all.1) # -0.0164286 -- bad model!
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.1, permutations = how(nperm=999)) # p = 0.611

## we can also do a permutation test by RDA axis
#anova(rda.all.1, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.all.1, by = "terms", permutations = how(nperm=999)) ### by variables
# ave wind speed near sig
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.1)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 4.705177              8.873171              4.300262              6.138801           1567.007750             61.555889
# Developed            Herbaceous                Mexico                Others             SaltonSea                 Shrub
# 1156.140790             31.666730             21.855612             46.656702            480.368028           2809.609836

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.b = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,9:10,35:37,39:40,42:44)]),
                     scope=formula(rda.all.1),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ ave.wind direction is sig, Developed is near sig

rda.all.b$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
# ave.wind_direction  1 257.35 2.3999  0.007 **

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.b2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,9:10,35:37,39:40,42:44)]),
                        scope=formula(rda.all.1),
                        permutations = how(nperm=999))
rda.all.b2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# drop Others, CropLand due to high p values in ordistep
rda.all.2<-rda(b.clr ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+Developed+Herbaceous+Mexico+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.2
summary(rda.all.2)

# how much variation does our model explain?
RsquareAdj(rda.all.2) # 0.008195629
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.2, permutations = how(nperm=999)) # p = 0.611

anova(rda.all.2, by = "terms", permutations = how(nperm=999)) ### by variables
# ave wind speed near sig
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.2)

# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand             Developed
# 3.447116              8.455555              2.816973              3.570053            210.867239            131.743267
# Herbaceous                Mexico             SaltonSea                 Shrub
# 15.943352             16.744903             61.141458            273.409823
## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.c = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,9:10,35,37,39:40,43:44)]),
                     scope=formula(rda.all.2),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ ave.wind direction is sig, Developed is near sig

rda.all.c$anova # see significance of individual terms in model
#                               Df    AIC      F Pr(>F)
# + ave.wind_direction  1 257.35 2.3999  0.004 **
# + Developed           1 257.70 1.5215  0.044 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.c2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,9:10,35,37,39:40,42:44)]),
                        scope=formula(rda.all.2),
                        permutations = how(nperm=999))
rda.all.c2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# dropping Mexico
rda.all.3<-rda(b.clr ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+Developed+Herbaceous+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.3
summary(rda.all.3)

# how much variation does our model explain?
RsquareAdj(rda.all.3) # 0.02010821
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.3, permutations = how(nperm=999)) # p = 0.611

anova(rda.all.3, by = "terms", permutations = how(nperm=999)) ### by variables
#                       Df Variance      F Pr(>F)
# ave.air_temp           1    493.7 1.4122  0.085 .
# ave.wind_speed         1    515.4 1.4745  0.078 .
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.3)

# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand             Developed
# 3.014973              7.675913              2.688341              3.564667             57.532573             41.092763
# Herbaceous             SaltonSea                 Shrub
# 8.800896             39.177838             82.458257

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.d = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,9:10,35,37,39,43:44)]),
                     scope=formula(rda.all.3),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ ave.wind direction is sig
rda.all.d$anova # see significance of individual terms in model
#                       Df    AIC      F Pr(>F)
# + ave.wind_direction  1 257.35 2.3999  0.005 **

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.d2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,9:10,35,37,39,43:44)]),
                        scope=formula(rda.all.3),
                        permutations = how(nperm=999))
rda.all.d2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# drop ave relative humidity due to large p values and low R^2 variation
rda.all.4<-rda(b.clr ~ ave.air_temp+ave.wind_speed+ave.wind_direction+BarrenLand+Developed+Herbaceous+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.4
summary(rda.all.4)

# how much variation does our model explain?
RsquareAdj(rda.all.4) # 0.01867898
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.4, permutations = how(nperm=999)) # p = 0.611

anova(rda.all.4, by = "terms", permutations = how(nperm=999)) ### by variables
#                       Df Variance      F Pr(>F)
# ave.air_temp        1    493.7 1.4101  0.089 .
# ave.wind_speed      1    515.4 1.4723  0.084 .
# ave.wind_direction  1    557.7 1.5932  0.063 .

## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.4)

# ave.air_temp     ave.wind_speed ave.wind_direction         BarrenLand          Developed         Herbaceous          SaltonSea
# 2.164729           7.450533           3.563559          54.243260          37.356764           8.795223          35.796114
# Shrub
# 64.205386

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.d = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,10,35,37,39,43:44)]),
                     scope=formula(rda.all.4),
                     direction = "forward",
                     permutations = how(nperm=999))
# b.clr ~ ave.wind direction + Developed is sig
rda.all.d$anova # see significance of individual terms in model
#                       Df    AIC      F Pr(>F)
#+ ave.wind_direction  1 257.35 2.3999  0.002 **
#+ Developed           1 257.70 1.5215  0.033 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.d2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5:6,10,35,37,39,43:44)]),
                        scope=formula(rda.all.4),
                        permutations = how(nperm=999))
rda.all.d2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# drop ave wind speed due to large p values and low R^2 variation
rda.all.5<-rda(b.clr ~ ave.air_temp+ave.wind_direction+BarrenLand+Developed+Herbaceous+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.5
summary(rda.all.5)

# how much variation does our model explain?
RsquareAdj(rda.all.5) # 0.02767908
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.5, permutations = how(nperm=999)) # p = 0.611

anova(rda.all.5, by = "terms", permutations = how(nperm=999)) ### by variables
#                       Df Variance      F Pr(>F)
# ave.air_temp        1    493.7 1.4232  0.094 .
# ave.wind_direction  1    757.5 2.1838  0.012 *

## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.5)
# ave.air_temp ave.wind_direction         BarrenLand          Developed         Herbaceous          SaltonSea              Shrub
# 2.164728           2.421237          45.974211          33.534243           5.190045          35.708238          63.679575

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.e = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,35,37,39,43:44)]),
                     scope=formula(rda.all.5),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.all.e$anova # see significance of individual terms in model
#                       Df    AIC      F Pr(>F)
# + ave.wind_direction  1 257.35 2.3999  0.006 **

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.e2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,35,37,39,43:44)]),
                        scope=formula(rda.all.5),
                        permutations = how(nperm=999))
rda.all.e2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# drop Herbaceous due to large p values and low R^2 variation
rda.all.6<-rda(b.clr ~ ave.air_temp+ave.wind_direction+BarrenLand+Developed+SaltonSea+Shrub,data=meta.all.scaled)

# check summary of RDA
rda.all.6
summary(rda.all.6)

# how much variation does our model explain?
RsquareAdj(rda.all.6) # 0.0364673
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.6, permutations = how(nperm=999)) # p = 0.611

anova(rda.all.6, by = "terms", permutations = how(nperm=999)) ### by variables
#                       Df Variance      F Pr(>F)
# ave.air_temp        1    493.7 1.4362  0.087 .
# ave.wind_direction  1    757.5 2.2037  0.012 *
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.6)
#ave.air_temp ave.wind_direction         BarrenLand          Developed          SaltonSea              Shrub
# 2.108749           2.411153          26.168186          22.005662          20.748237          38.909089

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.f = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,35,37,43:44)]),
                     scope=formula(rda.all.6),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.all.f$anova # see significance of individual terms in model
#                       Df    AIC      F Pr(>F)
# + ave.wind_direction  1 257.35 2.3999  0.006 **
# + Developed           1 257.70 1.5215  0.041 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.f2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,35,37,43:44)]),
                        scope=formula(rda.all.6),
                        permutations = how(nperm=999))
rda.all.f2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# nothing sig

# drop Shrub due to large p values and low R^2 variation
rda.all.7<-rda(b.clr ~ ave.air_temp+ave.wind_direction+BarrenLand+Developed+SaltonSea,data=meta.all.scaled)

# check summary of RDA
rda.all.7
summary(rda.all.7)

# how much variation does our model explain?
RsquareAdj(rda.all.7) # 0.05903603
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.7, permutations = how(nperm=999)) # p =

anova(rda.all.7, by = "terms", permutations = how(nperm=999)) ### by variables
#                       Df Variance      F Pr(>F)
# ave.air_temp        1    493.7 1.4706  0.076 .
# ave.wind_direction  1    757.5 2.2566  0.009 **
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.7)
#ave.air_temp ave.wind_direction         BarrenLand          Developed          SaltonSea
# 2.041044           1.879083           7.606135           4.054603           5.574042

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.g = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,35,37,43)]),
                     scope=formula(rda.all.7),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.all.g$anova # see significance of individual terms in model
#                       Df    AIC      F Pr(>F)
# + ave.wind_direction  1 257.35 2.3999  0.005 **

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.g2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,35,37,43)]),
                        scope=formula(rda.all.7),
                        permutations = how(nperm=999))
rda.all.g2$anova # see significance of individual terms in model
#                       R2.adj Df    AIC      F Pr(>F)
# + ave.wind_direction 0.049292  1 257.35 2.3999  0.009 **

# drop SaltonSea due to large p values and low R^2 variation
rda.all.8<-rda(b.clr ~ ave.air_temp+ave.wind_direction+BarrenLand+Developed,data=meta.all.scaled)

# check summary of RDA
rda.all.8
summary(rda.all.8)

# how much variation does our model explain?
RsquareAdj(rda.all.8) # 0.0710518
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.8, permutations = how(nperm=999)) # p =  0.008 **

anova(rda.all.8, by = "terms", permutations = how(nperm=999)) ### by variables
#                       Df Variance      F Pr(>F)
# ave.air_temp        1    493.7 1.4896  0.078 .
# ave.wind_direction  1    757.5 2.2858  0.003 **
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.8)
# ave.air_temp ave.wind_direction         BarrenLand          Developed
# 1.792784           1.818044           3.103824           3.985710

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.h = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,35,37)]),
                     scope=formula(rda.all.8),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.all.h$anova # see significance of individual terms in model
#                       Df    AIC      F Pr(>F)
# + ave.wind_direction  1 257.35 2.3999  0.004 **
# + Developed           1 257.70 1.5215  0.047 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.h2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,35,37)]),
                        scope=formula(rda.all.8),
                        permutations = how(nperm=999))
rda.all.h2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# + ave.wind_direction 0.049292  1 257.35 2.3999  0.005 **
# + Developed          0.067986  1 257.70 1.5215  0.044 *

# drop BarrenLand due to large p values and low R^2 variation
rda.all.9<-rda(b.clr ~ ave.air_temp+ave.wind_direction+Developed,data=meta.all.scaled)

# check summary of RDA
rda.all.9
summary(rda.all.9)

# how much variation does our model explain?
RsquareAdj(rda.all.9) # 0.07560572
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.9, permutations = how(nperm=999)) # p =  0.003 **

anova(rda.all.9, by = "terms", permutations = how(nperm=999)) ### by variables
#                       Df Variance      F Pr(>F)
# ave.air_temp        1    493.7 1.4970  0.064 .
# ave.wind_direction  1    757.5 2.2971  0.004 **
# Developed           1    466.4 1.4143  0.089 .
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.9)
# ave.air_temp ave.wind_direction          Developed
# 1.161867           1.295476           1.319659

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.i = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,37)]),
                     scope=formula(rda.all.9),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.all.i$anova # see significance of individual terms in model
#                       Df    AIC      F Pr(>F)
# + ave.wind_direction  1 257.35 2.3999  0.006 **
# + Developed           1 257.70 1.5215  0.041 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.i2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(5,10,37)]),
                        scope=formula(rda.all.9),
                        permutations = how(nperm=999))
rda.all.i2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)

# drop ave air temp due due to p value compared to other vars
rda.all.10<-rda(b.clr ~ ave.wind_direction+Developed,data=meta.all.scaled)

# check summary of RDA
rda.all.10
summary(rda.all.10)

# how much variation does our model explain?
RsquareAdj(rda.all.10) # 0.06798582
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.all.10, permutations = how(nperm=999)) # p =  **0.003

anova(rda.all.10, by = "terms", permutations = how(nperm=999)) ### by variables
#                       Df Variance      F Pr(>F)
# ave.wind_direction  1    813.9 2.4480  0.005 **
# Developed           1    505.9 1.5215  0.055 .
## this will help us interpret our RDA and we can see some variable are not significant

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.all.10)
# ave.wind_direction          Developed
# 1.175934           1.175934

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/

head(meta.all.scaled[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
# more info on ordistep & ordiR2step here: https://www.davidzeleny.net/anadat-r/doku.php/en:forward_sel_examples
rda.all.j = ordistep(rda(b.clr ~ 1, data = meta.all.scaled[,c(10,37)]),
                     scope=formula(rda.all.10),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.all.j$anova # see significance of individual terms in model
#                       Df    AIC      F Pr(>F)
# ave.wind_direction  1 257.35 2.3999   0.01 **

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.all.j2 = ordiR2step(rda(b.clr ~ 1, data = meta.all.scaled[,c(10,37)]),
                        scope=formula(rda.all.10),
                        permutations = how(nperm=999))
rda.all.j2$anova # see significance of individual terms in model
#                               R2.adj Df    AIC      F Pr(>F)
# + ave.wind_direction 0.049292  1 257.35 2.3999  0.006 **

#### RDA - WI ####

rownames(WI) %in% rownames(b.clr_WI) # check order of DFs
head(WI)

# included all Surface type frequencies and it overfit the model so pulling some out...
rda.WI.0<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=WI)

# check summary of RDA
rda.WI.0
summary(rda.WI.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.WI.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.WI.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.WI.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.WI.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.WI.0)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 18.91092              35.48104              24.87321              31.29241              11.32342              21.88877
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others
# NA                    NA                    NA                    NA                    NA                    NA
# SaltonSea                 Shrub
# NA                    NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.WI.a = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9:10,35:44)]),
                         scope=formula(rda.WI.0),
                         direction = "forward",
                         permutations = how(nperm=999))
rda.WI.a$anova # see significance of individual terms in model
# + OpenWater  1 63.384 2.0387  0.048 *

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.WI.a2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9:10,35:44)]),
                            scope=formula(rda.WI.0),
                            permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.WI.a, permutations = how(nperm=999)) #not significant

# Dropping Mexico & CropLand because because ordistep p values were really high
rda.WI.1<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+Developed+Herbaceous+OpenWater+Others+SaltonSea+Shrub,data=WI)
summary(rda.WI.1)
RsquareAdj(rda.WI.1) # how much variation is explained by our model?
anova(rda.WI.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.1)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand             Developed
# 167.20211              49.65380              19.15610              59.68793             103.67779              14.52182
# Herbaceous             OpenWater                Others             SaltonSea                 Shrub
# NA                    NA                    NA                    NA                    NA

head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.b1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9:10,35,37,39,41:44)]),
                          scope=formula(rda.WI.1),
                          direction = "forward",
                          permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# + OpenWater              1 63.384 2.0387  0.051 .
# + ave.air_temp           1 63.902 1.5367  0.098 .
rda.WI.b1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.b2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9:10,35,37,39,41:44)]),
                            scope=formula(rda.WI.1),
                            permutations = how(nperm=999))
# too many variables

rda.WI.b2$anova
# check best fit model based on above results
anova(rda.WI.b1, permutations = how(nperm=999)) #

anova(rda.WI.0, rda.WI.1) #

# Dropping OpenWater because its negligible across sites
rda.WI.2<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+Developed+Herbaceous+Others+SaltonSea+Shrub,data=WI)
summary(rda.WI.2)
RsquareAdj(rda.WI.2) # how much variation is explained by our model?
anova(rda.WI.2, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.2)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand             Developed
# 167.20211              49.65380              19.15610              59.68793             103.67779              14.52182
# Herbaceous             OpenWater                Others             SaltonSea                 Shrub
# NA                    NA                    NA                    NA                    NA

head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.c1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9:10,35,37,39,42:44)]),
                     scope=formula(rda.WI.2),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# ae air temp near sig
rda.WI.c1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.c2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9:10,35,37,39,42:44)]),
                       scope=formula(rda.WI.2),
                       permutations = how(nperm=999))
# too many variables

# Dropping ave wind direction and Others because of highest p values in ordistep fxn above
rda.WI.3<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+BarrenLand+Developed+Herbaceous+SaltonSea+Shrub,data=WI)
summary(rda.WI.3)
RsquareAdj(rda.WI.3) # how much variation is explained by our model?
anova(rda.WI.3, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.3)
# ave.air_temp        ave.wind_speed ave.relative_humidity            BarrenLand             Developed            Herbaceous
# 247.95167             175.33588             209.39141             667.46893              72.90821            2025.67906
# SaltonSea                 Shrub
# NA                    NA

head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.d1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,35,37,39,43:44)]),
                     scope=formula(rda.WI.3),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# ave air temp near sig
rda.WI.d1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.d2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,35,37,39,43:44)]),
                       scope=formula(rda.WI.3),
                       permutations = how(nperm=999))
# too many variables

# drop Herbaceous due to highest p value in ordistep
rda.WI.4<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+BarrenLand+Developed+SaltonSea+Shrub,data=WI)
summary(rda.WI.4)
RsquareAdj(rda.WI.4) # how much variation is explained by our model?
anova(rda.WI.4, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.4)
# ave.air_temp        ave.wind_speed ave.relative_humidity            BarrenLand             Developed             SaltonSea
# 18.553062             15.247675             11.772372             16.601723              1.935937             54.181397
# Shrub
# NA

head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.e1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,35,37,43:44)]),
                     scope=formula(rda.WI.4),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# ave air temp near sig
rda.WI.e1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.e2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,35,37,43:44)]),
                       scope=formula(rda.WI.4),
                       permutations = how(nperm=999))
# too many variables

# dropping Salton Sea - too high p value in ordistep
rda.WI.5<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+BarrenLand+Developed+Shrub,data=WI)
summary(rda.WI.5)
RsquareAdj(rda.WI.5) # how much variation is explained by our model?
anova(rda.WI.5, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.5)
# ave.air_temp        ave.wind_speed ave.relative_humidity            BarrenLand             Developed                 Shrub
# 286.39590              67.79934              79.06143             125.81347              14.17166             822.16423

head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.f1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,35,37,44)]),
                     scope=formula(rda.WI.5),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# ave air temp near sig
rda.WI.f1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.f2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,35,37,44)]),
                       scope=formula(rda.WI.5),
                       permutations = how(nperm=999))
# too many variables

rda.WI.6<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+BarrenLand+Developed,data=WI)
summary(rda.WI.6)
RsquareAdj(rda.WI.6) # how much variation is explained by our model? -0.02300137
anova(rda.WI.6, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.6)
# ave.air_temp        ave.wind_speed ave.relative_humidity            BarrenLand             Developed
# 4.915098              1.973293              2.676292              8.365654              1.930336

head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.g1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,35,37)]),
                     scope=formula(rda.WI.6),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# ave air temp near sig
rda.WI.g1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.g2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,35,37)]),
                       scope=formula(rda.WI.6),
                       permutations = how(nperm=999))
# too many variables

# dropped BarrenLand due to highest p value in ordistep
rda.WI.7<-rda(b.clr_WI ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+Developed,data=WI)
summary(rda.WI.7)
RsquareAdj(rda.WI.7) # how much variation is explained by our model? 0.05204929
anova(rda.WI.7, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.7)
# ave.air_temp        ave.wind_speed ave.relative_humidity             Developed
# 1.605232              1.934618              2.078610              1.166318

head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.h1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,37)]),
                     scope=formula(rda.WI.7),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# ave air temp near sig
rda.WI.h1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.h2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5:6,9,37)]),
                       scope=formula(rda.WI.7),
                       permutations = how(nperm=999))

# dropped ave.wind_speed due to largest p value in ordistep
rda.WI.8<-rda(b.clr_WI ~ ave.air_temp+ave.relative_humidity+Developed,data=WI)
summary(rda.WI.8)
RsquareAdj(rda.WI.8) # how much variation is explained by our model? 0.1918306
anova(rda.WI.8, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.8)
# ave.air_temp ave.relative_humidity             Developed
# 1.538701              1.492563              1.091951

head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.i1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5,9,37)]),
                     scope=formula(rda.WI.8),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# ave air temp near sig
rda.WI.i1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.i2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5,9,37)]),
                       scope=formula(rda.WI.8),
                       permutations = how(nperm=999))
#

# drop ave relative humidity due to highest p value in ordistep
rda.WI.9<-rda(b.clr_WI ~ ave.air_temp+Developed,data=WI)
summary(rda.WI.9)
RsquareAdj(rda.WI.9) # how much variation is explained by our model? 0.2063033
anova(rda.WI.9, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.WI.9)
# ave.air_temp    Developed
# 1.037627     1.037627

head(WI[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.WI.j1 = ordistep(rda(b.clr_WI ~ 1, data = WI[,c(5,37)]),
                     scope=formula(rda.WI.9),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# ave air temp near sig
rda.WI.j1$anova
#
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.WI.j2 = ordiR2step(rda(b.clr_WI ~ 1, data = WI[,c(5,37)]),
                       scope=formula(rda.WI.9),
                       permutations = how(nperm=999))
#

# best model for WI
rda.WI.9$call

#### RDA - DP ####

rownames(DP) %in% rownames(b.clr_DP) # check order of DFs
head(DP)

# included all Surface type frequencies and it overfit the model so pulling some out...
## dropping forest, open waters, others, and shrub first (after some experimenting)
rda.DP.0<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=DP)

# check summary of RDA
rda.DP.0
summary(rda.DP.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.DP.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.DP.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.DP.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.DP.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.DP.0)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 24003.065             14260.183              1110.023             86550.715            303149.675              1202.145
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others
# NA                    NA                    NA                    NA                    NA                    NA
# SaltonSea                 Shrub
# NA                    NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.DP.a = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:6,9:10,35:44)]),
                    scope=formula(rda.DP.0),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.DP.a$anova # see significance of individual terms in model
# Df    AIC      F Pr(>F)
# + BarrenLand  1 65.028 2.4284  0.001 ***

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.DP.a2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:6,9:10,35:44)]),
                       scope=formula(rda.DP.0),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.DP.a, permutations = how(nperm=999)) #not significant
#anova(rda.DP.a2, permutations = how(nperm=999)) # p =  0.001, significant

# Dropping Open Water & Forest because negligible values across samples
rda.DP.1<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Mexico+Others+SaltonSea+Shrub,data=DP)
summary(rda.DP.1)
RsquareAdj(rda.DP.1) # how much variation is explained by our model?
anova(rda.DP.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.DP.1)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 24003.065             14260.183              1110.023             86550.715            303149.675              1202.145
# Developed            Herbaceous                Mexico                Others             SaltonSea                 Shrub
# NA                    NA                    NA                    NA                    NA                    NA

head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.b1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:6,9:10,35:37,39:40,42:44)]),
                     scope=formula(rda.DP.1),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
#
rda.DP.b1$anova
#             Df    AIC      F Pr(>F)
# + BarrenLand  1 65.028 2.4284  0.002 **

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.b2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:6,9:10,35:37,39:40,42:44)]),
                       scope=formula(rda.DP.1),
                       permutations = how(nperm=999))
# too many variables

rda.DP.b2$anova
# check best fit model based on above results
anova(rda.DP.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.DP.0, rda.DP.1) #

# dropping Others
rda.DP.2<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Mexico+SaltonSea+Shrub,data=DP)
summary(rda.DP.2)
RsquareAdj(rda.DP.2) # how much variation is explained by our model? NA
anova(rda.DP.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.DP.2, by=NULL,permutations = how(nperm=999)) #

vif.cca(rda.DP.2)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 24003.065             14260.183              1110.023             86550.715            303149.675              1202.145
# Developed            Herbaceous                Mexico             SaltonSea                 Shrub
# NA                    NA                    NA                    NA                    NA

# check if ORP & Sulfide are significantly correlated in August, which they are [strong, sig negative corr]
#cor.test(meta.all.scaled[metadata$SampDate=="DP",]$Sulfide_microM, meta.all.scaled[metadata$SampDate=="DP",]$ORP_mV, method="pearson") # ******

head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.c1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:6,9:10,35:37,39:40,43:44)]),
                     scope=formula(rda.DP.2),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.DP.c1$anova
# Df    AIC      F Pr(>F)
# + BarrenLand  1 65.028 2.4284  0.001 ***

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.c2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:6,9:10,35:37,39:40,43:44)]),
                       scope=formula(rda.DP.2),
                       permutations = how(nperm=999))
# too many terms

# check best fit model based on above results
anova(rda.DP.c1, permutations = how(nperm=999)) #

# dropping Mexico
rda.DP.3<-rda(b.clr_DP ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+SaltonSea+Shrub,data=DP)
summary(rda.DP.3)
RsquareAdj(rda.DP.3) # how much variation is explained by our model? NA%
anova(rda.DP.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.DP.3, by=NULL,permutations = how(nperm=999)) # p =

vif.cca(rda.DP.3)
# #  ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 24003.065             14260.183              1110.023             86550.715            303149.675              1202.145
# Developed            Herbaceous             SaltonSea                 Shrub
# NA                    NA                    NA                    NA

head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.d1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5:6,9:10,35:37,39,42:44)]),
                     scope=formula(rda.DP.3),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.DP.d1$anova
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.d2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5:6,9:10,35:37,39,42:44)]),
                       scope=formula(rda.DP.3),
                       permutations = how(nperm=999))
# too many terms

# drop ave.wind_speed due to highest P value in last ordistep function
rda.DP.4<-rda(b.clr_DP ~ ave.air_temp+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+SaltonSea+Shrub,data=DP)
summary(rda.DP.4)
RsquareAdj(rda.DP.4) # how much variation is explained by our model? NA
anova(rda.DP.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
#
anova(rda.DP.4, by=NULL,permutations = how(nperm=999))

vif.cca(rda.DP.4)
# ave.air_temp ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand             SaltonSea
# 24.06550              15.85546              27.90553              47.91774             120.44522              35.67224
# Shrub
# NA
head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.e1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5,9:10,35:36,43:44)]),
                     scope=formula(rda.DP.4),
                     direction = "forward",
                     permutations = how(nperm=999))
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.e2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5,9:10,35:36,43:44)]),
                       scope=formula(rda.DP.4),
                       permutations = how(nperm=999))
# too many terms

# dropped CropLand due to high p value in last ordistep
rda.DP.5<-rda(b.clr_DP ~ ave.air_temp+ave.relative_humidity+ave.wind_direction+BarrenLand+SaltonSea+Shrub,data=DP)
summary(rda.DP.5)
RsquareAdj(rda.DP.5) # how much variation is explained by our model? NA
anova(rda.DP.5, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)

anova(rda.DP.5, by=NULL,permutations = how(nperm=999)) # p

vif.cca(rda.DP.5)
# ave.air_temp ave.relative_humidity    ave.wind_direction            BarrenLand             SaltonSea                 Shrub
# 106.67071              10.14565              29.29001            1094.86853             218.97169            1257.48280

head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.f1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(5,9:10,35,43:44)]),
                     scope=formula(rda.DP.5),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.DP.f1$anova
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.f2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(5,9:10,35,43:44)]),
                       scope=formula(rda.DP.5),
                       permutations = how(nperm=999))
# too many terms

# dropped ave air temp due to high p value in last ordistep
rda.DP.6<-rda(b.clr_DP ~ ave.relative_humidity+ave.wind_direction+BarrenLand+SaltonSea+Shrub,data=DP)
summary(rda.DP.6)
RsquareAdj(rda.DP.6) # how much variation is explained by our model? 0.5862192
anova(rda.DP.6, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.wind_direction     1   3323.5 4.5482  0.014 *
# SaltonSea              1   2159.1 2.9547  0.026 *

anova(rda.DP.6, by=NULL,permutations = how(nperm=999)) # p = 0.049

vif.cca(rda.DP.6)
# ave.relative_humidity    ave.wind_direction            BarrenLand             SaltonSea                 Shrub
# 4.55945              26.37536             187.94045              46.46174             265.83238

head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.g1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(9:10,35,43:44)]),
                     scope=formula(rda.DP.6),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.DP.g1$anova
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.g2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(9:10,35,43:44)]),
                       scope=formula(rda.DP.6),
                       permutations = how(nperm=999))
rda.DP.g2$anova
# R2.adj Df    AIC      F Pr(>F)
# + ave.wind_direction 0.19612  1 64.994 2.4638  0.013 *

# dropped ave relative_humidity due to high p value in last ordistep
rda.DP.7<-rda(b.clr_DP ~ ave.wind_direction+BarrenLand+SaltonSea+Shrub,data=DP)
summary(rda.DP.7)
RsquareAdj(rda.DP.7) # how much variation is explained by our model? 0.3888711
anova(rda.DP.7, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.wind_direction  1   3497.7 3.2409  0.025 *
# BarrenLand          1   1458.9 1.3517  0.293
# SaltonSea           1   1171.9 1.0858  0.418
# Shrub               1   2309.1 2.1395  0.085 .

anova(rda.DP.7, by=NULL,permutations = how(nperm=999)) # p = 0.097

vif.cca(rda.DP.7)
# ave.wind_direction         BarrenLand          SaltonSea              Shrub
# 22.33765           98.58229           39.55891          146.66815

head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.h1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(10,35,43:44)]),
                     scope=formula(rda.DP.7),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.DP.h1$anova
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.h2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(10,35,43:44)]),
                       scope=formula(rda.DP.7),
                       permutations = how(nperm=999))
rda.DP.h2$anova
# R2.adj Df    AIC      F Pr(>F)
# + ave.wind_direction 0.19612  1 64.994 2.4638  0.012 *

# dropped SaltonSea due to high p value in last ordistep
rda.DP.8<-rda(b.clr_DP ~ ave.wind_direction+BarrenLand+Shrub,data=DP)
summary(rda.DP.8)
RsquareAdj(rda.DP.8) # how much variation is explained by our model? 0.1611839
anova(rda.DP.8, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.wind_direction  1   3497.7 2.3612  0.021 *

anova(rda.DP.8, by=NULL,permutations = how(nperm=999)) # p = 0.193

vif.cca(rda.DP.8)
# ave.wind_direction         BarrenLand              Shrub
# 21.81808           57.61915           16.33774

head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.i1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(10,35,44)]),
                     scope=formula(rda.DP.8),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.DP.i1$anova
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.i2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(10,35,44)]),
                       scope=formula(rda.DP.8),
                       permutations = how(nperm=999))
rda.DP.i2$anova
# R2.adj Df    AIC      F Pr(>F)

# dropped Shrub due to high p value in last ordistep
rda.DP.9<-rda(b.clr_DP ~ ave.wind_direction+BarrenLand,data=DP)
summary(rda.DP.9)
RsquareAdj(rda.DP.9) # how much variation is explained by our model? 0.2016672
anova(rda.DP.9, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                             Df Variance      F Pr(>F)
# ave.wind_direction          1   3497.7 2.4809  0.009 **

anova(rda.DP.9, by=NULL,permutations = how(nperm=999)) # p = 0.026

vif.cca(rda.DP.9)
# ave.wind_direction         BarrenLand
# 13.03326           13.03326

head(DP[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.DP.j1 = ordistep(rda(b.clr_DP ~ 1, data = DP[,c(10,35)]),
                     scope=formula(rda.DP.9),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.DP.j1$anova
# BarrenLand
# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.DP.j2 = ordiR2step(rda(b.clr_DP ~ 1, data = DP[,c(10,35)]),
                       scope=formula(rda.DP.9),
                       permutations = how(nperm=999))
rda.DP.j2$anova
# R2.adj Df    AIC      F Pr(>F)

# final model:
rda.DP.9$call

#### RDA - BDC ####

rownames(BDC) %in% rownames(b.clr_BDC) # check order of DFs
head(BDC)

# included all Surface type frequencies and it overfit the model so pulling some out...
## dropping forest, open waters, others, and shrub first (after some experimenting)
rda.BDC.0<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=BDC)

# check summary of RDA
rda.BDC.0
summary(rda.BDC.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.BDC.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.BDC.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.BDC.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.BDC.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.BDC.0)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 248.82571             166.10742            2138.69051             895.18785            1550.43893              45.23175
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others
# NA                    NA                    NA                    NA                    NA                    NA
# SaltonSea                 Shrub
# NA                    NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.BDC.a = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,9:10,35:44)]),
                    scope=formula(rda.BDC.0),
                    direction = "forward",
                    permutations = how(nperm=999))
rda.BDC.a$anova # see significance of individual terms in model
# Df    AIC      F Pr(>F)
# BarrenLand, Developed, ave air temp, and Forest were all near sig

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.BDC.a2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,9:10,35:44)]),
                       scope=formula(rda.BDC.0),
                       permutations = how(nperm=999))
# too many terms

# Dropping SaltonSea and Mexico
rda.BDC.1<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+OpenWater+Others+Shrub,data=BDC)
summary(rda.BDC.1)
RsquareAdj(rda.BDC.1) # how much variation is explained by our model?
anova(rda.BDC.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.BDC.1)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 248.82571             166.10742            2138.69051             895.18785            1550.43893              45.23175
# Developed            Herbaceous                Forest             OpenWater                Others                 Shrub
# NA                    NA                    NA                    NA                    NA                    NA

head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.b1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,9:10,35:39,41:42,44)]),
                     scope=formula(rda.BDC.1),
                     direction = "forward",
                     permutations = how(nperm=999))
#                      Df    AIC      F Pr(>F)
# near sig...
# + BarrenLand             1 62.539 1.7548  0.062 .
# + ave.air_temp           1 62.680 1.6207  0.074 .
# + Developed              1 62.645 1.6538  0.086 .
# + Forest                 1 62.785 1.5223  0.100 .

rda.BDC.b1$anova

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.b2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,9:10,35:39,41:42,44)]),
                       scope=formula(rda.BDC.1),
                       permutations = how(nperm=999))
# too many variables

rda.BDC.b2$anova
# check best fit model based on above results
anova(rda.BDC.b1, permutations = how(nperm=999)) # p =  0.001, significant

anova(rda.BDC.0, rda.BDC.1) #

# dropping Shrub due to highest p val in ordistep
rda.BDC.2<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+OpenWater+Others,data=BDC)
summary(rda.BDC.2)
RsquareAdj(rda.BDC.2) # how much variation is explained by our model? NA
anova(rda.BDC.2, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.BDC.2, by=NULL,permutations = how(nperm=999)) #

vif.cca(rda.BDC.2)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 248.82571             166.10742            2138.69051             895.18785            1550.43893              45.23175
# Developed            Herbaceous                Forest             OpenWater
# NA                    NA                    NA                    NA

head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.c1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,9:10,35:39,41:42)]),
                     scope=formula(rda.BDC.2),
                     direction = "forward",
                     permutations = how(nperm=999))
# + BarrenLand, Developed, and ave air temp near sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.c2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5:6,9:10,35:39,41:42)]),
                       scope=formula(rda.BDC.2),
                       permutations = how(nperm=999))
# too many terms

# drop ave wind speed due to highest p value in ordistep
rda.BDC.3<-rda(b.clr_BDC ~ ave.air_temp+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+OpenWater+Others,data=BDC)
summary(rda.BDC.3)
RsquareAdj(rda.BDC.3) # how much variation is explained by our model? NA
anova(rda.BDC.3, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.BDC.3, by=NULL,permutations = how(nperm=999)) #

vif.cca(rda.BDC.3)
# ave.air_temp ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand             Developed
# 28.716429            133.126108             92.489689            115.088007              3.877069              3.488450
# Herbaceous                Forest             OpenWater                Others
# NA                    NA                    NA                    NA

head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.d1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5,9:10,35:39,41:42)]),
                      scope=formula(rda.BDC.3),
                      direction = "forward",
                      permutations = how(nperm=999))
# + BarrenLand, Developed, ave air temp, and Forest near sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.d2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5,9:10,35:39,41:42)]),
                        scope=formula(rda.BDC.3),
                        permutations = how(nperm=999))
# too many terms

# dropping OpenWater due to ordistep p value being highest
rda.BDC.4<-rda(b.clr_BDC ~ ave.air_temp+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+Others,data=BDC)
summary(rda.BDC.4)
RsquareAdj(rda.BDC.4) # how much variation is explained by our model? NA
anova(rda.BDC.4, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.BDC.4, by=NULL,permutations = how(nperm=999)) #

vif.cca(rda.BDC.4)
# ave.air_temp ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand             Developed
# 28.716429            133.126108             92.489689            115.088007              3.877069              3.488450
# Herbaceous                Forest     Others
# NA                    NA                     NA

head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.e1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5,9:10,35:39,42)]),
                      scope=formula(rda.BDC.4),
                      direction = "forward",
                      permutations = how(nperm=999))
# + BarrenLand, Developed, ave air temp, and Forest near sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.e2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5,9:10,35:39,42)]),
                        scope=formula(rda.BDC.4),
                        permutations = how(nperm=999))
# too many terms

# dropping ave relative humidity + CropLand due to ordistep p values being highest
rda.BDC.5<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_direction+BarrenLand+Developed+Herbaceous+Forest+Others,data=BDC)
summary(rda.BDC.5)
RsquareAdj(rda.BDC.5) # how much variation is explained by our model? NA
anova(rda.BDC.5, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.BDC.5, by=NULL,permutations = how(nperm=999)) #

vif.cca(rda.BDC.5)
# ave.air_temp ave.wind_direction         BarrenLand          Developed         Herbaceous             Forest             Others
# 198.2678           552.3799           579.1409           251.9826            85.2811          2823.2953                 NA

head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.f1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5,10,35,37:39,42)]),
                      scope=formula(rda.BDC.5),
                      direction = "forward",
                      permutations = how(nperm=999))
# + BarrenLand, Developed, ave air temp, and Forest near sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.f2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5,10,35,37:39,42)]),
                        scope=formula(rda.BDC.5),
                        permutations = how(nperm=999))
# too many terms

# dropping Others due to ordistep p values being highest
rda.BDC.6<-rda(b.clr_BDC ~ ave.air_temp+ave.wind_direction+BarrenLand+Developed+Herbaceous+Forest,data=BDC)
summary(rda.BDC.6)
RsquareAdj(rda.BDC.6) # how much variation is explained by our model? NA
anova(rda.BDC.6, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant

anova(rda.BDC.6, by=NULL,permutations = how(nperm=999)) #

vif.cca(rda.BDC.6)
# ave.air_temp ave.wind_direction         BarrenLand          Developed         Herbaceous             Forest
# 198.2678           552.3799           579.1409           251.9826            85.2811          2823.2953

head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.g1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5,10,35,37:39)]),
                      scope=formula(rda.BDC.6),
                      direction = "forward",
                      permutations = how(nperm=999))
# + BarrenLand, Developed, ave air temp, and Forest near sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.g2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5,10,35,37:39)]),
                        scope=formula(rda.BDC.6),
                        permutations = how(nperm=999))
# too many terms

# dropping Herbaceous + ave wind direction due to ordistep p values being highest
rda.BDC.7<-rda(b.clr_BDC ~ ave.air_temp+BarrenLand+Developed+Forest,data=BDC)
summary(rda.BDC.7)
RsquareAdj(rda.BDC.7) # how much variation is explained by our model? 0.1694793
anova(rda.BDC.7, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
# ave.air_temp  1   1653.1 1.7685  0.071 . (near sig)

anova(rda.BDC.7, by=NULL,permutations = how(nperm=999)) # p = 0.211

vif.cca(rda.BDC.7)
# ave.air_temp   BarrenLand    Developed       Forest
# 2.844049     5.871025     2.539155     7.509786

head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.h1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(5,35,37:38)]),
                      scope=formula(rda.BDC.7),
                      direction = "forward",
                      permutations = how(nperm=999))
# + BarrenLand, Developed, and Forest near sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.h2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(5,35,37:38)]),
                        scope=formula(rda.BDC.7),
                        permutations = how(nperm=999))
# BarrenLand is close

# drop ave.air_temp due to highest p value in ordistep
rda.BDC.8<-rda(b.clr_BDC ~ BarrenLand+Developed+Forest,data=BDC)
summary(rda.BDC.8)
RsquareAdj(rda.BDC.8) # how much variation is explained by our model? 0.2037182
anova(rda.BDC.8, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#            Df Variance      F Pr(>F)
# BarrenLand  1   1754.3 1.9575  0.040 *

anova(rda.BDC.8, by=NULL,permutations = how(nperm=999)) # p = 0.084

vif.cca(rda.BDC.8)
# BarrenLand  Developed     Forest
# 5.699074   2.533920   4.758348

head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.i1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(35,37:38)]),
                      scope=formula(rda.BDC.8),
                      direction = "forward",
                      permutations = how(nperm=999))
# + BarrenLand, Developed, near sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.i2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(35,37:38)]),
                        scope=formula(rda.BDC.8),
                        permutations = how(nperm=999))
# BarrenLand near sig

# drop Forest due to highest p value in ordistep
rda.BDC.9<-rda(b.clr_BDC ~ BarrenLand+Developed,data=BDC)
summary(rda.BDC.9)
RsquareAdj(rda.BDC.9) # how much variation is explained by our model? 0.1373921
anova(rda.BDC.9, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#            Df Variance      F Pr(>F)
# BarrenLand  1 1754.3 1.8070  0.044 *

anova(rda.BDC.9, by=NULL,permutations = how(nperm=999)) # p = 0.063

vif.cca(rda.BDC.9)
# BarrenLand  Developed
# 1.214723   1.214723

head(BDC[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.BDC.j1 = ordistep(rda(b.clr_BDC ~ 1, data = BDC[,c(35,37)]),
                      scope=formula(rda.BDC.9),
                      direction = "forward",
                      permutations = how(nperm=999))
# + BarrenLand, Developed, near sig

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.BDC.j2 = ordiR2step(rda(b.clr_BDC ~ 1, data = BDC[,c(35,37)]),
                        scope=formula(rda.BDC.9),
                        permutations = how(nperm=999))
# BarrenLand near sig

# BDC FINAL MODEL:
rda.BDC.9$call

#### RDA - PD ####

rownames(PD) %in% rownames(b.clr_PD) # check order of DFs
head(PD)

# included all Surface type frequencies and it overfit the model so pulling some out...
## dropping forest, open waters, others, and shrub first (after some experimenting)
rda.PD.0<-rda(b.clr_PD ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Forest+Herbaceous+Mexico+OpenWater+Others+SaltonSea+Shrub,data=PD)

# check summary of RDA
rda.PD.0
summary(rda.PD.0)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.PD.0) # NA
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
anova(rda.PD.0, permutations = how(nperm=999))

## we can also do a permutation test by RDA axis
#anova(rda.PD.0, by = "axis", permutations = how(nperm=999)) ### by RDA axis
## or by terms (aka variables)
anova(rda.PD.0, by = "terms", permutations = how(nperm=999)) ### by variables
## this will help us interpret our RDA and we can see some variable are not significant
#                               Df Variance      F Pr(>F)

# Calculating variance inflation factor (VIF) for each predictor variable to check multicolinearity of predictor variables
## VIF helps determien which predictors are too strongly correlated with other predictor variables to explain variation observed
vif.cca(rda.PD.0)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 41.338500             20.671797              9.674336             12.547748              7.979633             19.388377
# Developed                Forest            Herbaceous                Mexico             OpenWater                Others
# NA                    NA                    NA                    NA                    NA                    NA
# SaltonSea                 Shrub
# NA                    NA

## Understanding VIF results...
# A value of 1 indicates there is no correlation between a given predictor variable and any other predictor variables in the model.
# A value between 1 and 5 indicates moderate correlation between a given predictor variable and other predictor variables in the model, but this is often not severe enough to require attention.
# A value greater than 5 indicates potentially severe correlation between a given predictor variable and other predictor variables in the model. In this case, the coefficient estimates and p-values in the regression output are likely unreliable.
# when to ignore high VIF values: https://statisticalhorizons.com/multicollinearity/
head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important (by p values)
rda.PD.a = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(5:6,9:10,35:44)]),
                     scope=formula(rda.PD.0),
                     direction = "forward",
                     permutations = how(nperm=999))
rda.PD.a$anova # see significance of individual terms in model
# Df    AIC      F Pr(>F)

# can also use model seletion to pick most important variables by which increases variation (R^2) the most
rda.PD.a2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(5:6,9:10,35:44)]),
                        scope=formula(rda.PD.0),
                        permutations = how(nperm=999))
# too many terms

# Dropping SaltonSea and Mexico due to lack of contributions
rda.PD.1<-rda(b.clr_PD ~ ave.air_temp+ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+OpenWater+Others+Shrub,data=PD)
summary(rda.PD.1)
RsquareAdj(rda.PD.1) # how much variation is explained by our model?
anova(rda.PD.1, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.1)
# ave.air_temp        ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand
# 41.338500             20.671797              9.674336             12.547748              7.979633             19.388377
# Developed            Herbaceous                Forest             OpenWater                Others                 Shrub
# NA                    NA                    NA                    NA                    NA                    NA

head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.b1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(5:6,9:10,35:39,41:42,44)]),
                      scope=formula(rda.PD.1),
                      direction = "forward",
                      permutations = how(nperm=999))

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.b2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(5:6,9:10,35:39,41:42,44)]),
                        scope=formula(rda.PD.1),
                        permutations = how(nperm=999))
# too many variables

# dropping ave air temp and OpenWater due to very high P values in ordistep
rda.PD.2<-rda(b.clr_PD ~ ave.wind_speed+ave.relative_humidity+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+Others+Shrub,data=PD)
summary(rda.PD.2)
RsquareAdj(rda.PD.2) # how much variation is explained by our model?
anova(rda.PD.2, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.2)
# ave.wind_speed ave.relative_humidity    ave.wind_direction            BarrenLand              CropLand             Developed
# 31.27854              22.05587              41.29358             719.64333              76.76934            1149.82156
# Herbaceous                Forest                Others                 Shrub
# NA                    NA                    NA                    NA

head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.c1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(6,9:10,35:39,42,44)]),
                     scope=formula(rda.PD.2),
                     direction = "forward",
                     permutations = how(nperm=999))

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.c2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(6,9:10,35:39,42,44)]),
                       scope=formula(rda.PD.2),
                       permutations = how(nperm=999))
# too many variables

# dropping ave relative humidity due to very high P values in ordistep
rda.PD.3<-rda(b.clr_PD ~ ave.wind_speed+ave.wind_direction+BarrenLand+CropLand+Developed+Herbaceous+Forest+Others+Shrub,data=PD)
summary(rda.PD.3)
RsquareAdj(rda.PD.3) # how much variation is explained by our model?
anova(rda.PD.3, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.3)
# ave.wind_speed ave.wind_direction         BarrenLand           CropLand          Developed         Herbaceous             Forest
# 10.002152           5.098642         493.689075          21.001346         232.988151         135.562206                 NA
# Others              Shrub
# NA                 NA

head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.d1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(6,10,35:39,42,44)]),
                     scope=formula(rda.PD.3),
                     direction = "forward",
                     permutations = how(nperm=999))

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.d2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(6,10,35:39,42,44)]),
                       scope=formula(rda.PD.3),
                       permutations = how(nperm=999))
# too many variables

# dropping BarrenLand due to very high P values in ordistep
rda.PD.4<-rda(b.clr_PD ~ ave.wind_speed+ave.wind_direction+CropLand+Herbaceous+Forest+Others+Shrub,data=PD)
summary(rda.PD.4)
RsquareAdj(rda.PD.4) # how much variation is explained by our model?
anova(rda.PD.4, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.4)
# ave.wind_speed ave.wind_direction           CropLand         Herbaceous             Forest             Others              Shrub
# 24.387609           2.193757          17.361640          13.118227          44.395684          30.004394                 NA

head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.e1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(6,10,36:39,42,44)]),
                     scope=formula(rda.PD.4),
                     direction = "forward",
                     permutations = how(nperm=999))

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.e2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(6,10,36:39,42,44)]),
                       scope=formula(rda.PD.4),
                       permutations = how(nperm=999))
# too many variables

# dropping Shrub and Herbaceous due to very high P values in ordistep
rda.PD.5<-rda(b.clr_PD ~ ave.wind_speed+ave.wind_direction+CropLand+Developed+Forest+Others,data=PD)
summary(rda.PD.5)
RsquareAdj(rda.PD.5) # how much variation is explained by our model?
anova(rda.PD.5, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.5)
# ave.wind_speed ave.wind_direction           CropLand          Developed             Forest             Others
# 1980.3353           270.5016           604.8987         14616.5694         15165.4746          3724.0446

head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.f1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(6,10,36:38,42)]),
                     scope=formula(rda.PD.5),
                     direction = "forward",
                     permutations = how(nperm=999))

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.f2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(6,10,36:38,42)]),
                       scope=formula(rda.PD.5),
                       permutations = how(nperm=999))
# too many variables

# dropping Developed due to very high P values in ordistep
rda.PD.6<-rda(b.clr_PD ~ ave.wind_speed+ave.wind_direction+CropLand+Forest+Others,data=PD)
summary(rda.PD.6)
RsquareAdj(rda.PD.6) # how much variation is explained by our model? 0.2019816
anova(rda.PD.6, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.6)
# ave.wind_speed ave.wind_direction           CropLand             Forest             Others
# 16.291695           1.733389           8.062481          23.297245          20.576821

head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.g1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(6,10,36,38,42)]),
                     scope=formula(rda.PD.6),
                     direction = "forward",
                     permutations = how(nperm=999))

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.g2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(6,10,36,38,42)]),
                       scope=formula(rda.PD.6),
                       permutations = how(nperm=999))
# ave wind speed is nearest sig but not that sig

# dropping Others due to low contribution
rda.PD.7<-rda(b.clr_PD ~ ave.wind_speed+ave.wind_direction+CropLand+Forest,data=PD)
summary(rda.PD.7)
RsquareAdj(rda.PD.7) # how much variation is explained by our model? 0.2296914
anova(rda.PD.7, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.7)
#

head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.h1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(6,10,36,38)]),
                     scope=formula(rda.PD.7),
                     direction = "forward",
                     permutations = how(nperm=999))

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.h2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(6,10,36,38)]),
                       scope=formula(rda.PD.7),
                       permutations = how(nperm=999))
#

# dropping ave wind direction due to very high p value identified by ANOVA
rda.PD.8<-rda(b.clr_PD ~ ave.wind_speed+CropLand+Forest,data=PD)
summary(rda.PD.8)
RsquareAdj(rda.PD.8) # how much variation is explained by our model? 0.2296914
anova(rda.PD.8, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)
# Df Variance      F Pr(>F)
# ave.wind_speed  1   2405.9 1.6588  0.112
# CropLand        1   2417.5 1.6668  0.103
# Forest          1   2791.7 1.9248  0.083 .

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.8)
# ave.wind_speed       CropLand         Forest
# 11.97063        5.20556       20.49904

head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.i1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(6,36,38)]),
                     scope=formula(rda.PD.8),
                     direction = "forward",
                     permutations = how(nperm=999))

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.i2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(6,36,38)]),
                       scope=formula(rda.PD.8),
                       permutations = how(nperm=999))
#

# dropping Forest due to negatie R^2 values from ordiR2step
rda.PD.9<-rda(b.clr_PD ~ ave.wind_speed+CropLand,data=PD)
summary(rda.PD.9)
RsquareAdj(rda.PD.9) # how much variation is explained by our model? 0.1046263
anova(rda.PD.9, by = "terms", permutations = how(nperm=999)) ### by variables
#                             Df Variance      F Pr(>F)
# Df Variance      F Pr(>F)
# ave.wind_speed  1   2405.9 1.3473  0.190
# CropLand        1   2417.5 1.3538  0.226
# Residual        4   7142.8

## this will help us interpret our RDA and we can see some variable are not significant
vif.cca(rda.PD.9)
# ave.wind_speed       CropLand
# 1.135983       1.135983

head(PD[,c(5:6,9:10,35:44)])
## we can use model selection instead of picking variables we think are important -- based on p values
rda.PD.j1 = ordistep(rda(b.clr_PD ~ 1, data = PD[,c(6,36)]),
                     scope=formula(rda.PD.9),
                     direction = "forward",
                     permutations = how(nperm=999))

# Can also use model selection to pick variables by which ones increase variation (R^2)
rda.PD.j2 = ordiR2step(rda(b.clr_PD ~ 1, data = PD[,c(6,36)]),
                       scope=formula(rda.PD.9),
                       permutations = how(nperm=999))
#

# final RDA for PD:
rda.PD.9$call

#### Final RDAs ####
# RDA by sampling timepoint
head(meta.all.scaled)
head(b.clr)
rownames(b.clr) %in% rownames(meta.all.scaled) # sanity check 1

# all data
#rda.all9$call # best model for all data

rda.all<-rda(b.clr ~ ave.air_temp + ave.wind_direction + Developed,data=meta.all.scaled)
rda.all
summary(rda.all)
RsquareAdj(rda.all) # how much variation is explained by our model? 0.07560572
anova(rda.all, permutations = how(nperm=999)) # p-value = 0.007 **
anova(rda.all, by = "terms", permutations = how(nperm=999))
#                               Df Variance      F Pr(>F)
# ave.air_temp        1    493.7 1.4970  0.054 .
# ave.wind_direction  1    757.5 2.2971  0.009 **
# Developed           1    466.4 1.4143  0.090 .
# Residual           24   7914.6

# p value adjusted for each term
aov.rda.all.terms<-anova(rda.all, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.all.terms$`Pr(>F)`,method="bonferroni",n=length(aov.rda.all.terms)) # adjusted pvalues
# [1] 0.336 0.016 0.356    NA

# p value adjusted for entire model
aov.rda.all<-anova(rda.all, by = NULL, permutations = how(nperm=999))
p.adjust(aov.rda.all$`Pr(>F)`,method="bonferroni",n=length(aov.rda.all)) # adjusted pvalues
# [1] 0.02   NA

# WI
#rda.WI.10$call # best model

rda.WI<-rda(b.clr_WI ~ ave.air_temp+Developed,data=WI)
summary(rda.WI)
RsquareAdj(rda.WI) # how much variation is explained by our model? 0.2063033
anova(rda.WI, permutations = how(nperm=999)) # p-value = 0.122
anova(rda.WI, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
# ave.air_temp  1   1866.4 1.7772  0.094 .
# Developed     1   1871.9 1.7824  0.144
# Residual        4   4774.3
aov.rda.WI<-anova(rda.WI, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.WI$`Pr(>F)`,method="bonferroni",n=length(aov.rda.WI)) # adjusted pvalues
# [1] 0.432 0.596    NA

# DP
#rda.DP.9$call # best model from above

rda.DP<-rda(b.clr_DP ~ ave.wind_direction + BarrenLand,data=DP)
summary(rda.DP)
RsquareAdj(rda.DP) # how much variation is explained by our model? 0.2016672
anova(rda.DP, permutations = how(nperm=999)) # p-value = 0.039
anova(rda.DP, by = "terms", permutations = how(nperm=999))
#                 Df Variance      F Pr(>F)
# ave.wind_direction  1   3497.7 2.4809  0.013 *
#   BarrenLand          1   1458.9 1.0348  0.428
# Residual            4   5639.5

aov.rda.DP<-anova(rda.DP, by = NULL, permutations = how(nperm=999))
p.adjust(aov.rda.DP$`Pr(>F)`,method="bonferroni",n=length(aov.rda.DP)) # adjusted pvalues
# 0.2

aov.rda.DP.terms<-anova(rda.DP, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.DP.terms$`Pr(>F)`,method="bonferroni",n=length(aov.rda.DP.terms)) # adjusted pvalues
# 0.044 1.000    NA

# BDC
#rda.BDC.9$call - best model

rda.BDC<-rda(b.clr_BDC ~ BarrenLand + Developed,data=BDC)
summary(rda.BDC)
RsquareAdj(rda.BDC) # how much variation is explained by our model? 0.1373921
anova(rda.BDC, permutations = how(nperm=999)) # p-value = 0.052
anova(rda.BDC, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
# BarrenLand  1   1754.3 1.8070  0.047 *
#   Developed   1   1115.2 1.1487  0.241
# Residual    4   3883.5

aov.rda.BDC<-anova(rda.BDC, by = NULL, permutations = how(nperm=999))
p.adjust(aov.rda.BDC$`Pr(>F)`,method="bonferroni",n=length(aov.rda.BDC)) # adjusted pvalues
# [1] 0.208    NA
aov.rda.BDC.terms<-anova(rda.BDC, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.BDC.terms$`Pr(>F)`,method="bonferroni",n=length(aov.rda.BDC.terms)) # adjusted pvalues
# 0.196 1.000    NA

# best model : rda.PD.9$call
rda.PD<-rda(b.clr_PD ~ ave.wind_speed+CropLand,data=PD)
summary(rda.PD)
RsquareAdj(rda.PD) # how much variation is explained by our model? 0.1046263
anova(rda.PD, permutations = how(nperm=999)) # p-value = 0.156
anova(rda.PD, by = "terms", permutations = how(nperm=999))
#                           Df Variance      F Pr(>F)
# ave.wind_speed  1   2405.9 1.3473  0.189
# CropLand        1   2417.5 1.3538  0.247
# Residual        4   7142.8

aov.rda.PD<-anova(rda.PD, by = NULL, permutations = how(nperm=999))
p.adjust(aov.rda.PD$`Pr(>F)`,method="bonferroni",n=length(aov.rda.PD)) # adjusted pvalues
#
aov.rda.PD.terms<-anova(rda.PD, by = "terms", permutations = how(nperm=999))
p.adjust(aov.rda.PD.terms$`Pr(>F)`,method="bonferroni",n=length(aov.rda.PD.terms)) # adjusted pvalues
#

# save RDAs as R object
save.image("data/Amplicon/SSD_Amplicon_EnvDriver_RDAsOnly.Rdata")

#### Plot RDA - ALL data ####
#plot(rda.WI) # depending on how many species you have, this step may take a while
plot(rda.all, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.all, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.all)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.all) # 0.07560572 aka 7.56%
## ^^ use this b/c chance correlations can inflate R^2

# we can then test for significance of the model by permutation
# if it is not significant, it doesn't matter how much of the variation is explained
##anova(rda.all, permutations = how(nperm=999)) # p = 0.001, significant

png('figures/EnvDrivers/SSD_AllData_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.all, arrows = TRUE,data = rda.all ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first


# variance partitioning of RDA
rda.all.part<-varpart(b.clr, meta.all.scaled$ave.wind_direction, meta.all.scaled$Developed,meta.all.scaled$ave.air_temp)
rda.all.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSD_AllData_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.all.part,
     Xnames = c("Ave Wind Direction", "Developed STF","Ave Air Temp"), # name the partitions
     bg = c("magenta", "gray","red2"), alpha = 80, # colour the circles
     digits = 3, # only show 2 digits
     cex = 1.5)
dev.off()

rda.sum.all<-summary(rda.all)
rda.sum.all$sites[,1:2]
rda.sum.all$cont #cumulative proportion of variance per axis
# RDA1 = 9.80, RDA2 = 5.30

# create data frame w/ RDA axes for sites
# first check rownames of RDA & metadata, then make df
rownames(rda.sum.all$sites) %in% rownames(meta.all.scaled)
rda.axes.all1<-data.frame(RDA1=rda.sum.all$sites[,1], RDA2=rda.sum.all$sites[,2], SampleID=rownames(rda.sum.all$sites), Site=meta.all.scaled$Site, SampDate=meta.all.scaled$SampDate)

# then merge with metadata to get all category colors!
rda.axes.all<-merge(rda.axes.all1,meta.all.scaled,by=c("SampleID","Site","SampDate"))

# create data frame w/ RDA axes for variables
arrows.all<-data.frame(RDA1=rda.sum.all$biplot[,1], RDA2=rda.sum.all$biplot[,2], Label=rownames(rda.sum.all$biplot))
arrows.all$Label[(arrows.all$Label) == "ave.air_temp"] <- "Ave Air Temp"
arrows.all$Label[(arrows.all$Label) == "ave.wind_direction"] <- "Ave Wind Dir."
arrows.all$Label[(arrows.all$Label) == "Developed"] <- "Developed STF"

rda.sum.all$cont #cumulative proportion of variance per axis
# RDA1 = 9.80, RDA2 = 5.30

rda.plot1<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot2<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate,shape=Site),size=4) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*5.5, yend = RDA2*5.5),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*7, y = RDA2*7, fontface="bold"), size=4)+
  coord_fixed(ratio = 1, xlim = c(-8,8), ylim = c(-8,8)) + theme_classic() +
  scale_color_manual(name ="Collection Date",values=unique(rda.axes.all$SampDate_Color[order(rda.axes.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_manual(values = c(7,10, 15,16)) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [9.80%]") + ylab("RDA2 [5.30%]")

ggsave(rda.plot2,filename = "figures/EnvDrivers/SSD_16S_RDA_AllData.png", width=10, height=10, dpi=600)


rda.plot3<-ggplot(rda.axes.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate,shape=Site),size=5) +
  geom_segment(data = arrows.all,mapping = aes(x = 0, y = 0, xend = RDA1*6, yend = RDA2*6),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.all,aes(label = Label, x = RDA1*8, y = RDA2*8, fontface="bold"), size=5)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.all$SampDate_Color[order(rda.axes.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "August 2021", "September 2021", "December 2021")) +
  scale_shape_discrete(labels=c("WI","DP","BDC"),name="Sample Date") +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [9.80%]") + ylab("RDA2 [5.30%]")

ggsave(rda.plot3,filename = "figures/EnvDrivers/SSD_16S_RDA_AllData_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - WI ####
#plot(rda.WI) # depending on how many species you have, this step may take a while
plot(rda.WI, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.WI, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.WI)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.WI) # 0.2063033
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSD_WI_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.WI, arrows = TRUE,data = rda.WI ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.WI.part<-varpart(b.clr_WI, WI$ave.air_temp, WI$Developed)
rda.WI.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSD_WI_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.WI.part,
     Xnames = c("Ave. Air Temp", "Developed STF"), # name the partitions
     bg = c("dodgerblue", "gray"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.WI<-summary(rda.WI)
rda.sum.WI$sites[,1:2]
rda.sum.WI$cont #cumulative proportion of variance per axis
# RDA1=28.01%, RDA2=19.08%

# create data frame w/ RDA axes for sites
rda.axes.WI<-data.frame(RDA1=rda.sum.WI$sites[,1], RDA2=rda.sum.WI$sites[,2], SampleID=rownames(rda.sum.WI$sites), Site=WI$Site)

# then merge with metadata to get all category colors!
rda.axes.WI.all<-merge(rda.axes.WI,WI,by=c("SampleID","Site"))

# create data frame w/ RDA axes for variables
arrows.WI<-data.frame(RDA1=rda.sum.WI$biplot[,1], RDA2=rda.sum.WI$biplot[,2], Label=rownames(rda.sum.WI$biplot))
arrows.WI$Label[(arrows.WI$Label) == "ave.air_temp"] <- "Ave Air Temp"
arrows.WI$Label[(arrows.WI$Label) == "Developed"] <- "Developed STF"

rda.plot5<-ggplot(rda.axes.WI.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.WI,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.WI,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.WI.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=4) +
  geom_segment(data = arrows.WI,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.WI,aes(label = Label, x = RDA1*9.85, y = RDA2*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.WI.all$SampDate_Color[order(rda.axes.WI.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, WI",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [28.01%]") + ylab("RDA2 [19.08%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/SSD_16S_RDA_WI.png", width=16, height=12, dpi=600)

rda.plot6b<-ggplot(rda.axes.WI.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=5) +
  geom_segment(data = arrows.WI,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.WI,aes(label = Label, x = RDA1*9, y = RDA2*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.WI.all$SampDate_Color[order(rda.axes.WI.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, WI",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [28.01%]") + ylab("RDA2 [19.08%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/SSD_16S_RDA_WI_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - DP ####
#plot(rda.DP) # depending on how many species you have, this step may take a while
plot(rda.DP, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.DP, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.DP)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.DP) # 0.2016672
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSD_DP_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.DP, arrows = TRUE,data = rda.DP ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.DP.part<-varpart(b.clr_DP, DP$ave.wind_direction, DP$BarrenLand)
rda.DP.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSD_DP_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.DP.part,
     Xnames = c("Ave. Wind Dir.", "Barren Land STF"), # name the partitions
     bg = c("magenta", "peachpuff2"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.DP<-summary(rda.DP)
rda.sum.DP$sites[,1:2]
rda.sum.DP$cont #cumulative proportion of variance per axis
# RDA1=33.25%, RDA2=13.53%

# create data frame w/ RDA axes for sites
rda.axes.DP<-data.frame(RDA1=rda.sum.DP$sites[,1], RDA2=rda.sum.DP$sites[,2], SampleID=rownames(rda.sum.DP$sites), Site=DP$Site)

# then merge with metadata to get all category colors!
rda.axes.DP.all<-merge(rda.axes.DP,DP,by=c("SampleID","Site"))

# create data frame w/ RDA axes for variables
arrows.DP<-data.frame(RDA1=rda.sum.DP$biplot[,1], RDA2=rda.sum.DP$biplot[,2], Label=rownames(rda.sum.DP$biplot))
arrows.DP$Label[(arrows.DP$Label) == "ave.wind_direction"] <- "Ave Wind Dir"
arrows.DP$Label[(arrows.DP$Label) == "BarrenLand"] <- "Barren Land STF"

rda.plot5<-ggplot(rda.axes.DP.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.DP,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.DP,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.DP.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=4) +
  geom_segment(data = arrows.DP,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.DP,aes(label = Label, x = RDA1*9.85, y = RDA2*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.DP.all$SampDate_Color[order(rda.axes.DP.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "August 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, DP",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [33.25%]") + ylab("RDA2 [13.53%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/SSD_16S_RDA_DP.png", width=16, height=12, dpi=600)

rda.plot6b<-ggplot(rda.axes.DP.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=5) +
  geom_segment(data = arrows.DP,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.DP,aes(label = Label, x = RDA1*9, y = RDA2*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.DP.all$SampDate_Color[order(rda.axes.DP.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "August 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, DP",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [33.25%]") + ylab("RDA2 [13.53%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/SSD_16S_RDA_DP_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - BDC ####
#plot(rda.BDC) # depending on how many species you have, this step may take a while
plot(rda.BDC, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.BDC, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.BDC)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.BDC) # 0.1373921
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSD_BDC_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.BDC, arrows = TRUE,data = rda.BDC ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.BDC.part<-varpart(b.clr_BDC, BDC$BarrenLand, BDC$Developed)
rda.BDC.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSD_BDC_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.BDC.part,
     Xnames = c("Ave. Air Temp", "Developed STF"), # name the partitions
     bg = c("peachpuff2", "lightgray"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.BDC<-summary(rda.BDC)
rda.sum.BDC$sites[,1:2]
rda.sum.BDC$cont #cumulative proportion of variance per axis
# RDA1=31.19%, RDA2=11.31%

# create data frame w/ RDA axes for sites
rda.axes.BDC<-data.frame(RDA1=rda.sum.BDC$sites[,1], RDA2=rda.sum.BDC$sites[,2], SampleID=rownames(rda.sum.BDC$sites), Site=BDC$Site)

# then merge with metadata to get all category colors!
rda.axes.BDC.all<-merge(rda.axes.BDC,BDC,by=c("SampleID","Site"))

# create data frame w/ RDA axes for variables
arrows.BDC<-data.frame(RDA1=rda.sum.BDC$biplot[,1], RDA2=rda.sum.BDC$biplot[,2], Label=rownames(rda.sum.BDC$biplot))
arrows.BDC$Label[(arrows.BDC$Label) == "BarrenLand"] <- "Barren Land STF"
arrows.BDC$Label[(arrows.BDC$Label) == "Developed"] <- "Developed STF"

rda.plot5<-ggplot(rda.axes.BDC.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.BDC,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.BDC,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.BDC.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=4) +
  geom_segment(data = arrows.BDC,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.BDC,aes(label = Label, x = RDA1*9.85, y = RDA2*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.BDC.all$SampDate_Color[order(rda.axes.BDC.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, BDC",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [31.19%]") + ylab("RDA2 [11.31%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/SSD_16S_RDA_BDC.png", width=16, height=12, dpi=600)

rda.plot6b<-ggplot(rda.axes.BDC.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=5) +
  geom_segment(data = arrows.BDC,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.BDC,aes(label = Label, x = RDA1*9, y = RDA2*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.BDC.all$SampDate_Color[order(rda.axes.BDC.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, BDC",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [31.19%]") + ylab("RDA2 [11.31%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/SSD_16S_RDA_BDC_bigger.png", width=15, height=15, dpi=600)

#### Plot RDA - PD ####
#plot(rda.PD) # depending on how many species you have, this step may take a while
plot(rda.PD, scaling = 1)
## scaling = 1 -> emphasizes relationships among sites
plot(rda.PD, scaling = 2)
## scaling = 2 -> emphasizes relationships among species

# check summary of RDA
summary(rda.PD)

# how much variation does our model explain?
## reminder: R^2 = % of variation in dependent variable explained by model
RsquareAdj(rda.PD) # 0.1046263
## ^^ use this b/c chance correlations can inflate R^2

png('figures/EnvDrivers/SSD_PD_autoplot_rda_example.png',width = 700, height = 600, res=100)
autoplot(rda.PD, arrows = TRUE,data = rda.PD ,layers=c("biplot","sites"),label = FALSE, label.size = 3, shape = FALSE, loadings = TRUE, loadings.colour = 'blue', loadings.label = TRUE, loadings.label.size = 3, scale= 0)+theme_classic()
dev.off()
## FOR AUTOPLOT -> must load packagve ggvegan first

# variance partitioning of RDA
rda.PD.part<-varpart(b.clr_PD, PD$ave.wind_speed, PD$CropLand)
rda.PD.part$part
# plot variance partitioning results
png('figures/EnvDrivers/SSD_PD_RDA_VariancePartitioning.png',width = 900, height = 900, res=100)
plot(rda.PD.part,
     Xnames = c("Ave. Wind Speed", "Crop Land STF"), # name the partitions
     bg = c("blue1", "gold1"), alpha = 80, # colour the circles
     digits = 3, # only show 3 digits
     cex = 1.5)
dev.off()

rda.sum.PD<-summary(rda.PD)
rda.sum.PD$sites[,1:2]
rda.sum.PD$cont #cumulative proportion of variance per axis
# RDA1=20.89%, RDA2=19.42%

# create data frame w/ RDA axes for sites
rda.axes.PD<-data.frame(RDA1=rda.sum.PD$sites[,1], RDA2=rda.sum.PD$sites[,2], SampleID=rownames(rda.sum.PD$sites), Site=PD$Site)

# then merge with metadata to get all category colors!
rda.axes.PD.all<-merge(rda.axes.PD,PD,by=c("SampleID","Site"))

# create data frame w/ RDA axes for variables
arrows.PD<-data.frame(RDA1=rda.sum.PD$biplot[,1], RDA2=rda.sum.PD$biplot[,2], Label=rownames(rda.sum.PD$biplot))
arrows.PD$Label[(arrows.PD$Label) == "ave.wind_speed"] <- "Ave Wind Speed"
arrows.PD$Label[(arrows.PD$Label) == "CropLand"] <- "Crop Land STF"

rda.plot5<-ggplot(rda.axes.PD.all, aes(x = RDA1, y = RDA2)) + geom_point(size=2) +
  geom_segment(data = arrows.PD,mapping = aes(x = 0, y = 0, xend = RDA1, yend = RDA2),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.5,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.PD,aes(label = Label, x = RDA1, y = RDA2, fontface="bold"))+
  coord_fixed() + theme_classic() +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1))

rda.plot6<-ggplot(rda.axes.PD.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=4) +
  geom_segment(data = arrows.PD,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 0.8,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.PD,aes(label = Label, x = RDA1*9.85, y = RDA2*9.5, fontface="bold"), size=4)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.PD.all$SampDate_Color[order(rda.axes.PD.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, PD",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [20.89%]") + ylab("RDA2 [19.42%]")

ggsave(rda.plot6,filename = "figures/EnvDrivers/SSD_16S_RDA_PD.png", width=16, height=12, dpi=600)

rda.plot6b<-ggplot(rda.axes.PD.all, aes(x = RDA1, y = RDA2)) + geom_point(aes(color=SampDate),size=5) +
  geom_segment(data = arrows.PD,mapping = aes(x = 0, y = 0, xend = RDA1*8, yend = RDA2*8),lineend = "round", # See available arrow types in example above
               linejoin = "round",
               size = 1,
               arrow = arrow(length = unit(0.15, "inches")),
               colour = "black") +
  geom_label(data = arrows.PD,aes(label = Label, x = RDA1*9, y = RDA2*9.5, fontface="bold"), size=5)+
  coord_fixed(ratio = 1) + theme_classic() + scale_color_manual(name ="Collection Date",values=unique(rda.axes.PD.all$SampDate_Color[order(rda.axes.PD.all$SampDate)]),labels=c("July 2020", "August 2020", "October 2020","November 2020", "July 2021", "September 2021", "December 2021")) +
  theme(axis.title.x = element_text(size=13),axis.title.y = element_text(size=13),axis.text = element_text(size=11),axis.text.x = element_text(vjust=1)) +
  labs(title="RDA: Bacteria/Archaea Composition in Salton Sea Dust, PD",subtitle="Using Centered-Log Ratio Data",color="Depth (m)") +
  xlab("RDA1 [20.89%]") + ylab("RDA2 [19.42%]")

ggsave(rda.plot6b,filename = "figures/EnvDrivers/SSD_16S_RDA_PD_bigger.png", width=15, height=15, dpi=600)

#### Save Progress ####

save.image("data/SSD_Amplicon_EnvDriver.Rdata")

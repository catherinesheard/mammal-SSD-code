library(MCMCglmm)
library(phangorn)

###################
#Preparing the data
###################

#Reading in the trees, and selecting a smaller subset of them
trees<-read.tree("trees1k.nex")
t100<-trees[1:100]
tree<-t100[[1]]

#Define which families of traits will be used as predictors
predictors<-c(
  "BodyMass",
  "Diet.Plants",
  "MatingSystem.updated",
  "PrecipVar",
  "TempVar"
)

#Read in and prepare the dataset. An order to look at can be specified here;
#if not, the whole mammalia dataset will be used.
setwd("C:/Users/cs16502/Dropbox/LOCKDOWN WOOT/Bristol/Mammals/SSD/Code")
source("Prepare_Data.R")
x<-Prepare_Data(phylo=tree,filename="data_2023.csv",variables=predictors)
x<-subset(x,MatingSystem.updated!="cooperative_breeder") #drop cooperative breeders

#trim out everything from the tree that's not in the dataset
t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$UphamTreeName.full))
#t100 and x should now match



####################
# Categorical Models
####################

#These models will loop through as many trees as are passed to them, so ensure
#the correct number of trees are input
#MAKE SURE TO CHANGE SAVE FILE NAMES


#Separate effects model without seasonality measures
source("MCMCglmm_MassDietMating.R")
MCMCglmm_MassDietMating_Separate(phylo=t100,data=x,savefile="MassDietMating_2023-fulldata.Rdata")

#Separate effects model with seasonality measures
source("MCMCglmm_MassDietMatingSeason.R")
MCMCglmm_MassDietMatingSeason_Separate(phylo=t100,data=x,savefile="MassDietMatingSeason_2023-fulldata.Rdata")
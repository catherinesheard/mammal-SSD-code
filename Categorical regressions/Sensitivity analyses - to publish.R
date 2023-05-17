library(MCMCglmm)
library(phangorn)

###################
#Preparing the data
###################

#Reading in the trees, and taking the first 100
trees<-read.tree("trees1k.nex")
t100<-trees[1:100]
tree<-t100[[1]]

#Define which families of variables will be used in this model
predictors<-c(
  "BodyMass",
  "Diet.Plants",
  "MatingSystem.updated",
  "PrecipVar",
  "TempVar"
)

#Read in and prepare the dataset
source("Prepare_Data.R")
x<-Prepare_Data(phylo=tree,filename="data_2023.csv",variables=predictors)
x<-subset(x,MatingSystem.updated!="cooperative_breeder")

#trim out everything from the tree that's not in the dataset
t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$UphamTreeName.full))
#t100 and x should now match

#If there's less than a 5% difference in mean body mass, override the SSD code to be monomorphic
x$zSSD[which(x$mass.sub.5==1)]<-"n"

####################
# Categorical Models
####################

#These models will loop through as many trees as are passed to them, so ensure
#the correct number of trees are input
#MAKE SURE TO CHECK YOUR SAVE FILE NAMES

source("MCMCglmm_MassDietMatingSeason_not_social.R")
MCMCglmm_MassDietMatingSeason_Separate(phylo=t100,data=x,savefile="MassDietMatingSeason_2023-5percent.Rdata")


###################
#Preparing the data
###################

#Reading in the trees, and taking the first 100
trees<-read.nexus("trees1k.nex")
t100<-trees[1:100]
tree<-t100[[1]]

#Define which families of variables will be used in this model
predictors<-c(
  "BodyMass",
  "Diet.Plants",
  "MatingSystem.updated",
  "PrecipVar",
  "TempVar"
)

#Read in and prepare the dataset
source("Prepare_Data.R")
x<-Prepare_Data(phylo=tree,filename="data_2023.csv",variables=predictors)
x<-subset(x,MatingSystem.updated!="cooperative_breeder")

#trim out everything from the tree that's not in the dataset
t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$UphamTreeName.full))
#t100 and x should now match

#If there's less than a 20% difference in mean body mass, override the SSD code to be monomorphic

x$zSSD[which(x$mass.sub.20==1)]<-"n"

####################
# Categorical Models
####################

#These models will loop through as many trees as are passed to them, so ensure
#the correct number of trees are input
#MAKE SURE TO CHECK YOUR SAVE FILE NAMES

source("MCMCglmm_MassDietMatingSeason_not_social.R")
MCMCglmm_MassDietMatingSeason_Separate(phylo=t100,data=x,savefile="MassDietMatingSeason_2023-20percent.Rdata")




###################
#Preparing the data
###################

#Reading in the trees, and taking the first 100
trees<-read.nexus("trees1k.nex")
t100<-trees[1:100]
tree<-t100[[1]]

#Define which families of variables will be used in this model
predictors<-c(
  "BodyMass",
  "Diet.Plants",
  "MatingSystem.updated",
  "PrecipVar",
  "TempVar"
)

#Read in and prepare the dataset
source("Prepare_Data.R")
x<-Prepare_Data(phylo=tree,filename="data_2023.csv",variables=predictors)
x<-subset(x,MatingSystem.updated!="cooperative_breeder")

#trim out everything from the tree that's not in the dataset
t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$UphamTreeName.full))
#t100 and x should now match

#If there's a disagreement between  the body mass and linear measurement values, set those species to monomorphic

x$zSSD[which(x$diagreement==1)]<-"n"

####################
# Categorical Models
####################

#These models will loop through as many trees as are passed to them, so ensure
#the correct number of trees are input
#MAKE SURE TO CHECK YOUR SAVE FILE NAMES

source("MCMCglmm_MassDietMatingSeason.R")
MCMCglmm_MassDietMatingSeason_Separate(phylo=t100,data=x,savefile="MassDietMatingSeason_2023-disagree.Rdata")
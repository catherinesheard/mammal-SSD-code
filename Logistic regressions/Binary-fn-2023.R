library(picante)
library(MCMCglmm)
library(phangorn)

###############
#Preparing data
###############


#Loading and subsetting trees
trees<-read.tree("trees1k.nex")
t100<-trees[1:100]
tree<-t100[[1]]

#Loading the data
x<-read.csv("data_2023.csv")
x$UphamTreeName.full<-paste(gsub(" ","_",x$PhyloName),toupper(x$Family),toupper(x$Order),sep="_")

#Dropping relevant rows: those with NA in key columns, and any with no SSD
if(sum(is.na(x$PhyloName)>0)){x<-x[-which(is.na(x$PhyloName)),]}
if(sum(is.na(x$dimorphic.size)>0)){x<-x[-which(is.na(x$dimorphic.size)),]}
if(sum(x$dimorphic.size=="m")>0){x<-x[-which(x$dimorphic.size=="m"),]}
if(sum(is.na(x$BodyMass)>0)){x<-x[-which(is.na(x$BodyMass)),]}
if(sum(is.na(x$Diet.Plants)>0)){x<-x[-which(is.na(x$Diet.Plants)),]}
if(sum(is.na(x$MatingSystem.updated)>0)){x<-x[-which(is.na(x$MatingSystem.updated)),]}
if(sum(is.na(x$PrecipVar)>0)){x<-x[-which(is.na(x$PrecipVar)),]}
if(sum(is.na(x$TempVar)>0)){x<-x[-which(is.na(x$TempVar)),]}
if(sum(x$MatingSystem.updated=="cooperative_breeder")>0){x<-x[-which(x$MatingSystem.updated=="cooperative_breeder"),]}


#Removing rows with duplicate species names
x<-x[!duplicated(x$UphamTreeName.full),]

#Removing any species in the sample which are not in the phylogeny
x[,1]<-x$UphamTreeName.full
row.names(x)<-x$UphamTreeName.full

trimmed<-match.phylo.data(tree,x)
x<-trimmed$data

#Removing any species in the phylogeny which are not in the sample
t100<-lapply(t100,drop.tip,tip=setdiff(tree$tip.label,x$UphamTreeName.full))

#t100 and x should now match

#Scaling relevant variables
x$zMass<-scale(log(as.numeric(x$BodyMass)))
x$zPlants<-scale(as.numeric(x$Diet.Plants))
x$zVert<-scale(as.numeric(x$Diet.Vert))
x$zInvert<-scale(as.numeric(x$Diet.Invert))
x$zTempVar<-scale(as.numeric(x$TempVar))
x$zPrecipVar<-scale(as.numeric(x$PrecipVar))

x$zSSD<-as.factor(x$dimorphic.size)

x$zSSD<-relevel(as.factor(x$zSSD),"n")
x$MatingSystem.updated<-relevel(as.factor(x$MatingSystem.updated),"monogamous")


###################
# Categorical Model
###################


i=1

tree<-t100[[i]] 
tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE)

animalA<-inverseA(tree)$Ainv 

#Setting up priors
gelmanprior<-list(B=list(mu=rep(0,10),
                         V=gelman.prior(~zMass+zPlants+zVert+zPrecipVar+zTempVar+
                                          as.factor(MatingSystem.updated), #if you edit this script, you'll have to adjust this equation to match your linear model, and mu will have to match the number of parameters being estimated
                                        data = x,  scale=1+pi^2/3)),
                  R=list(V=1,fix=1),G=list(G1=list(V=1E-10,nu=-1)))

#First model run
model.structure<-MCMCglmm(zSSD~zMass+zPlants+zVert+zPrecipVar+zTempVar+
                            as.factor(MatingSystem.updated), 
                          random=~UphamTreeName.full, 
                          ginverse=list(UphamTreeName.full=animalA), 
                          prior = gelmanprior, 
                          verbose=TRUE, 
                          family="categorical",
                          data = x,
                          nitt=55000*2,
                          thin=100,
                          burnin=5000*2,
                          pl=TRUE,
                          pr=TRUE,
                          slice=TRUE)

final.model<-model.structure
final.model$VCV[((i-1)*10+1):(i*10), ]<-model.structure$VCV[1:10,] 
final.model$Sol[((i-1)*10+1):(i*10), ]<-model.structure$Sol[1:10,] 
final.model$Liab[((i-1)*10+1):(i*10), ]<-model.structure$Liab[1:10,] 

nsamp.t<-nrow(model.structure$VCV)
start1.t=list(R=model.structure$VCV[nsamp.t,"units"], G=list(G1=model.structure$VCV[nsamp.t,"UphamTreeName.full"]))

setwd("C:/Users/tanys/Documents/Bristol MSci/Internship/MCMCglmm/Results")
save(final.model,file="MassDietMatingSeason_fn-2023.Rdata")



#Subsequent model runs
for(i in 1:100){
  tree<-t100[[i]]  
  tree<-nnls.tree(cophenetic(tree),tree,rooted=TRUE) #the tree is non-ultrametric, which is annoying, so we force it to be ultrametric
  
  
  animalA<-inverseA(tree)$Ainv 
  
  
  model.structure<-MCMCglmm(zSSD~zMass+zPlants+zVert+zPrecipVar+zTempVar+
                              as.factor(MatingSystem.updated), 
                            random=~UphamTreeName.full, 
                            ginverse=list(UphamTreeName.full=animalA), 
                            prior = gelmanprior, 
                            start= start1.t,
                            verbose=FALSE, 
                            family="categorical",
                            data = x,
                            nitt=55000,
                            thin=5000,
                            burnin=5000,
                            pl=TRUE,
                            pr=TRUE,
                            slice=TRUE)
  
  print(i)
  
  final.model$VCV[((i-1)*10+1):(i*10), ]<-model.structure$VCV[1:10,] 
  final.model$Sol[((i-1)*10+1):(i*10), ]<-model.structure$Sol[1:10,] 
  final.model$Liab[((i-1)*10+1):(i*10), ]<-model.structure$Liab[1:10,] 
  
  nsamp.t<-nrow(model.structure$VCV)
  start1.t=list(R=model.structure$VCV[nsamp.t,"units"], G=list(G1=model.structure$VCV[nsamp.t,"UphamTreeName.full"]))
  
  
  save(final.model,file="MassDietMatingSeason_fn-2023.Rdata")
  
  
}

save(final.model,file="MassDietMatingSeason_fn-2023.Rdata")
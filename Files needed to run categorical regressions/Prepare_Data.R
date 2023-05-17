Prepare_Data<-function(phylo,filename,variables,order=NULL){
  #phylo: the phylogeny the data will have to match
  #filename: file the data should be drawn from
  #variables: which variables will be used in the later model
  #order: if the model is being run for a single order, this will define the order
  
  library(picante)
  
  #Read in data
  data<-read.csv(filename)
  if(!is.null(order)){
    data<-subset(data,Order==order)
  }
  data$UphamTreeName.full<-paste(gsub(" ","_",data$PhyloName),toupper(data$Family),toupper(data$Order),sep="_")
  
  #Remove rows where the specified columns are NA
  if(sum(is.na(data$PhyloName)>0)){data<-data[-which(is.na(data$PhyloName)),]}
  if(sum(is.na(data$dimorphic.size)>0)){data<-data[-which(is.na(data$dimorphic.size)),]}
  for (i in 1:length(variables)){
    if(sum(is.na(data[,variables[i]]>0))){data<-data[-which(is.na(data[,variables[i]])),]}
  }
  
  data<-data[!duplicated(data$UphamTreeName.full),]
  data[,1]<-data$UphamTreeName.full
  row.names(data)<-data$UphamTreeName.full
  
  #Remove species which are not in the phylogeny
  trimmed<-match.phylo.data(phylo,data)
  data<-trimmed$data
  
  #Scale all of the potential explanatory variables
  data$zMass<-scale(log(as.numeric(data$BodyMass)))
  data$zPlants<-scale(as.numeric(data$Diet.Plants))
  data$zVert<-scale(as.numeric(data$Diet.Vert))
  data$zInvert<-scale(as.numeric(data$Diet.Invert))
  data$zTempVar<-scale(as.numeric(data$TempVar))
  data$zPrecipVar<-scale(as.numeric(data$PrecipVar))
  
  
  data$zSSD<-as.factor(data$dimorphic.size)
  
  return(data)
  
}
MCMCglmm_MassDietMating_Separate<-function(phylo,data,savefile){
  # phylo = a distribution of 100 trees, trimmed to match the data
  # data = a dataframe containing the species-level data
  # savefile = name of the file the results will be saved to
  
  i=1
  
  singletree<-phylo[[i]]
  
  #making the tree ultrametric
  singletree<-nnls.tree(cophenetic(singletree),singletree,rooted=TRUE)
  
  #Setting up variance/covariance matrix
  animalA<-inverseA(singletree)$Ainv 
  
  #Setting priors
  gelmanprior<-list(R=list(V=diag(2),fix=1),G=list(G1=list(V=1E-10,nu=-1)))
  
  #calculating number of iterations needing to be run
  no.it<-(length(phylo)*1000)+10000
  if (no.it<110000){
    no.it<-110000
  }
  
  #Releveling categorical variables
  data$zSSD<-relevel(as.factor(data$zSSD),"n")
  data$MatingSystem.updated<-relevel(as.factor(data$MatingSystem.updated),"monogamous")
  data$socialorg<-relevel(as.factor(data$socialorg),"group")
  
  #Preliminary model run
  model.structure<-MCMCglmm(zSSD~trait:zMass+trait:zPlants+trait:zVert+trait:as.factor(MatingSystem.updated)+trait:as.factor(socialorg)+trait, 
                            random=~UphamTreeName.full, 
                            ginverse=list(UphamTreeName.full=animalA), 
                            prior = gelmanprior, 
                            verbose=TRUE, 
                            family="categorical",
                            rcov = ~us(trait):units,
                            data = data,
                            nitt=no.it,
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
  start1.t=list(R=list(V=diag(2),fix=1), G=list(G1=model.structure$VCV[nsamp.t,"UphamTreeName.full"]))
 
  save(final.model,file=savefile)
  
  #Main model runs
  for(i in 1:length(phylo)){
    singletree<-phylo[[i]]  
    singletree<-nnls.tree(cophenetic(singletree),singletree,rooted=TRUE) #the tree is non-ultrametric, which is annoying, so we force it to be ultrametric
    
    animalA<-inverseA(singletree)$Ainv 
    
    model.structure<-MCMCglmm(zSSD~trait:zMass+trait:zPlants+trait:zVert+trait:as.factor(MatingSystem.updated)+trait:as.factor(socialorg)+trait, 
                              random=~UphamTreeName.full, 
                              ginverse=list(UphamTreeName.full=animalA), 
                              prior = gelmanprior, 
                              start= start1.t,
                              verbose=FALSE, 
                              family="categorical",
                              rcov = ~us(trait):units,
                              data = data,
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
    start1.t=list(R=diag(2), G=list(G1=model.structure$VCV[nsamp.t,"UphamTreeName.full"]))
    
    
    
    save(final.model,file=savefile)
    
    
  }
  
  save(final.model,file=savefile)
}
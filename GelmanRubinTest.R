
# First part of the code gets Gelman Rubin diagnostic for Mean behavior variables. Second part for the 144 
# covariance terms of the Phylo matrix corresponding to the 12 most common behaviors. Third part
# the Gelman Rubin test for the covariance terms of the individual matrix corresponding to the same 12 most
#common behaviors.

library(coda,lib.loc="Rpackages")

ListPsrf=c()
for(k in 1:134) {
FileName=paste("~/ownCloud/flytree/data/Anc134BhModel0PrTrueAllMinusYak00",1,sep="")
load(FileName)
Temp=model0$Sol
ChainMatrix=Temp[,k]

for(j in 2:20) {
  
  FileName=paste("~/ownCloud/flytree/data/Anc134BhModel0PrTrueAllMinusYak00",j,sep="")
  load(FileName)
  Temp=model0$Sol
  ChainT=Temp[,k]
  ChainMatrix=rbind(ChainT,ChainMatrix)
  rm(Temp)
}

ChainsMCMC=as.mcmc.list(lapply(as.data.frame(ChainMatrix), mcmc)) 
PSRFTemp=gelman.diag(ChainsMCMC, confidence = 0.95, transform=FALSE, autoburnin=FALSE,
            multivariate=TRUE)

ListPsrf=cbind(c(PSRFTemp$psrf[[1]], PSRFTemp$psrf[[2]]),ListPsrf)
print(ListPsrf)
rm(ChainT,ChainMatrix)
}

writeMat("~/ownCloud/flytree/data/MeanVariablesGelmanRubinTest.mat", ListPsrf=ListPsrf)

# Gelman Rubin for the 10% most common behaviors for Phylo

FileName="~/ownCloud/flytree/data/IndexPhyloHighBhForR.mat"
Index=readMat(FileName)

ListPsrfVCVPhylo=c()
for(k in 1:length(Index$InR)) {
  Filename=paste("~/ownCloud/flytree/data/Anc134BhModel0PrTrueAllMinusYak00",1,sep="")
  load(Filename)
  Temp=model0$VCV
  ChainMatrix=Temp[,Index$InR[[k]]]
  
  for(j in 2:20) {
    
    FileName=paste("~/ownCloud/flytree/data/Anc134BhModel0PrTrueAllMinusYak00",j,sep="")
    load(FileName)
    Temp=model0$VCV
    ChainT=Temp[,Index$InR[[k]]]
    ChainMatrix=rbind(ChainT,ChainMatrix)
    rm(Temp)
  }
  
  ChainsMCMC=as.mcmc.list(lapply(as.data.frame(ChainMatrix), mcmc)) 
  PSRFTempVCV=gelman.diag(ChainsMCMC, confidence = 0.95, transform=TRUE, autoburnin=FALSE,
                       multivariate=TRUE)
  
  ListPsrfVCVPhylo=cbind(c(PSRFTempVCV$psrf[[1]], PSRFTempVCV$psrf[[2]]),ListPsrfVCVPhylo)
  print(ListPsrfVCVPhylo)
  rm(ChainT,ChainMatrix)
}

writeMat("~/ownCloud/flytree/data/VCVPhyloVariablesGelmanRubinTest.mat", ListPsrfVCVPhylo=ListPsrfVCVPhylo)

# Gelman Rubin for the 10% most common behaviors for Individual matrix

FileName="~/ownCloud/flytree/data/IndexIndiHighBhForR.mat"
Index=readMat(FileName)

ListPsrfVCVIndi=c()
for(k in 1:length(Index$InR)) {
  Filename=paste("~/ownCloud/flytree/data/Anc134BhModel0PrTrueAllMinusYak00",1,sep="")
  load(Filename)
  Temp=model0$VCV
  ChainMatrix=Temp[,Index$InR[[k]]]
  
  for(j in 2:20) {
    
    FileName=paste("~/ownCloud/flytree/data/Anc134BhModel0PrTrueAllMinusYak00",j,sep="")
    load(FileName)
    Temp=model0$VCV
    ChainT=Temp[,Index$InR[[k]]]
    ChainMatrix=rbind(ChainT,ChainMatrix)
    rm(Temp)
  }
  
  ChainsMCMC=as.mcmc.list(lapply(as.data.frame(ChainMatrix), mcmc)) 
  PSRFTempVCV=gelman.diag(ChainsMCMC, confidence = 0.95, transform=TRUE, autoburnin=FALSE,
                          multivariate=TRUE)
  
  ListPsrfVCVIndi=cbind(c(PSRFTempVCV$psrf[[1]], PSRFTempVCV$psrf[[2]]),ListPsrfVCVIndi)
  print(ListPsrfVCVIndi)
  rm(ChainT,ChainMatrix)
}

writeMat("~/ownCloud/flytree/data/VCVIndiVariablesGelmanRubinTest.mat", ListPsrfVCVIndi=ListPsrfVCVIndi)



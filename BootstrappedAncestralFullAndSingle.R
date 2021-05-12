# ===Ancestral behaviors of melanogaster fly===
# =============================================
#MODEL WHERE COVARIANCE MATRICES ARE ONLY DIAGONAL. SO BEHAVIORAL TRAITS ARE NOT CO-EVOLVING IN THIS MODEL.

#library(MCMCglmm)
library(coda,lib.loc="Rpackages")
library(ape,lib.loc="Rpackages")
library(Matrix,lib.loc="Rpackages")
library(MCMCglmm) # in the r console:install.packages('MCMCglmm')
library(parallel)
##library(parallel,lib.loc="/Rpackages/")
#Install R.matlab package from the Pacakges in the right bottom window

#FUNCTIONS:
#SAMPLES USING MCMC
SimulatedSamplesSpSingleModel<-function(MCMCSol) {
  
  d=dim(MCMCSol)
  # Random Effects
  DataAnc=MCMCSol[,1]
  Rand1 = MCMCSol[,2] #Ind Random effects oldYasa
  Rand2 = MCMCSol[,3] #Ind Random effects oldMela
  Rand3 = MCMCSol[,4] #Ind Random effects oldSim
  Rand4 = MCMCSol[,5] #Ind Random effects Yakuba
  Rand5 = MCMCSol[,6] #Ind Random effects Santomea
  Rand6 = MCMCSol[,7] #Ind Random effects Sechellia
  Rand7 = MCMCSol[,8] #Ind Random effects Simulans
  Rand8 = MCMCSol[,9] #Ind Random effects Mauritiana
  Rand9 = MCMCSol[,10] #Ind Random effects Melanogaster
  
  SamplesYak = DataAnc+Rand1+Rand4  #Yak
  SamplesSant = DataAnc+Rand1+Rand5  #Sant
  SamplesSech = DataAnc+Rand2+Rand3+Rand6  #Sech
  SamplesSim = DataAnc+Rand2+Rand3+Rand7   #Sim
  SamplesMau = DataAnc+Rand2+Rand3+Rand8  #Mau
  SamplesMel = DataAnc+Rand2+Rand9  #Mel
  
  output<-list(SamplesYak,SamplesSant,SamplesSech,SamplesSim,SamplesMau,SamplesMel)
  return(output)
}

SimulatedSamplesSp<-function(MCMCSol) {
  
  d=dim(MCMCSol)
  # Random Effects
  DataAnc=MCMCSol[,1:134]
  In1 = seq(from = 136, to = d[2], by = 9) #Ind Random effects oldMela 
  Rand1 = MCMCSol[,In1]
  In2 = seq(from = 137, to = d[2], by = 9) #Ind Random effects oldSim
  Rand2 = MCMCSol[,In2]
  In3 = seq(from = 140, to = d[2], by = 9) #Ind Random effects Sech
  Rand3 = MCMCSol[,In3]
  In4 = seq(from = 141, to = d[2], by = 9) #Ind Random effects Sim
  Rand4 = MCMCSol[,In4]
  In5 = seq(from = 142, to = d[2], by = 9) #Ind Random effects Mau
  Rand5 = MCMCSol[,In5]
  In6 = seq(from = 135, to = d[2], by = 9) #Ind Random effects oldyasa
  Rand6 = MCMCSol[,In6]
  In7 = seq(from = 138, to = d[2], by = 9) #Ind Random effects Yak
  Rand7 = MCMCSol[,In7] 
  In8 = seq(from = 139, to = d[2], by = 9) #Ind Random effects Sant
  Rand8 = MCMCSol[,In8]
  In9 = seq(from = 143, to = d[2], by = 9) #Ind Random effects Mel
  Rand9 = MCMCSol[,In9]
  
  SamplesYak = DataAnc+Rand6+Rand7  #Yak
  SamplesSant = DataAnc+Rand6+Rand8  #Sant
  SamplesSech = DataAnc+Rand1+Rand2+Rand3  #Sech
  SamplesSim = DataAnc+Rand1+Rand2+Rand4   #Sim
  SamplesMau = DataAnc+Rand1+Rand2+Rand5  #Mau
  SamplesMel = DataAnc+Rand1+Rand9  #Mel
  
  output<-list(SamplesYak,SamplesSant,SamplesSech,SamplesSim,SamplesMau,SamplesMel)
  return(output)
}

JSDist<-function(Samples1, Samples2) {
  
  Dt = 0.01
  MAX = max(Samples1,Samples2)
  MIN = min(Samples1,Samples2)
  PS1 = hist(Samples1, breaks = seq(from = MIN, to = MAX+Dt, by =Dt), plot=FALSE)
  PS2 = hist(Samples2, breaks = seq(from = MIN, to = MAX+Dt, by =Dt),plot=FALSE)
  P1 = ((PS1$counts)+1)/(sum(PS1$counts) + length(PS1$breaks)-1)
  P2 = ((PS2$counts)+1)/(sum(PS2$counts) + length(PS2$breaks)-1)
  
  PM = (P1+P2)/2
  
  JSD = sum(P1*log(P1/PM))/2 + sum(P2*log(P2/PM))/2
  return(JSD)
}


# ---Defining phylogeny -> pedigree variable for MCMCglmm
tt="((yakuba:0.09,santomea:0.09)oldyasa:0.44,((sechelia:0.15,simulans:0.15,mauritiana:0.15)oldsim:0.14,melanogaster:0.29)oldmela:0.24);"


flytree<-read.tree(text=tt)

# ---Load the dataset: l_i= log p(behavior_i) - log p(behavior_0) 
# 134 behaviors + 0 behavior, 593 flies from 6 species (column-> animal)

LFlydat <- read.table("../../data/logDatFlySmoothAllMinusYak00NewNames.txt", header=TRUE, sep="\t", row.names="id")

# ---Defining prior (from mulTree, github, and R-sig-ME) nu==n
#phen.var<-var(LFlydat[,Beh+1])
#prior11<-list(G=list(G1=list(V=phen.var/4, n=1.002)), R=list(V=phen.var/4, n=1.002))


#---Run MCMCglmm
NN=10
DICSin = matrix(nrow = NN, ncol = 134)
DICCom = matrix(nrow = NN, ncol = 1)
JSDSinSum = matrix(nrow = NN, ncol = 134)
JSDComSum = matrix(nrow = NN, ncol = 134)

Nb= 134

IJ <- (1/(Nb+1))*(diag(Nb)+matrix(1,Nb,Nb))
priorCom<-list(G=list(G1=list(V=IJ/2,n=Nb)), R=list(V=IJ/2,n=Nb))


for(k in 1:NN){
  
  BootsDat=LFlydat[sample(nrow(LFlydat),size=dim(LFlydat)[1],replace=TRUE),]  
  
  IndYak = (BootsDat[,1] == 'yakuba')
  IndSant = (BootsDat[,1] == 'santomea')
  IndMel = (BootsDat[,1] == 'melanogaster')
  IndSech = (BootsDat[,1] == 'sechelia')
  IndSim = (BootsDat[,1] == 'simulans')
  IndMau = (BootsDat[,1] == 'mauritiana')
  
  SamplesExpYak = BootsDat[IndYak,]
  SamplesExpSant = BootsDat[IndSant,]
  SamplesExpSech = BootsDat[IndSech,]
  SamplesExpSim = BootsDat[IndSim,]
  SamplesExpMau = BootsDat[IndMau,]
  SamplesExpMel = BootsDat[IndMel,]
  
  SamplesExpSp <- c("SamplesExpYak","SamplesExpSant","SamplesExpSech","SamplesExpSim","SamplesExpMau","SamplesExpMel")
  
  # COMPLEX BEHAVIOR MODEL
  
  
  sfix="cbind("
  for (i in 1:(Nb-1)){sfix <-paste(sfix,colnames(BootsDat)[i+1],",",sep="")}
  sfix <-paste(sfix,colnames(BootsDat)[Nb+1],") ~ trait-1",sep="")
  fixed <- as.formula(sfix)
  #idh only fits the variance terms but not the off diagonal terms, unlike us that fits the entire matrix.
  # ---Run MCMCglmm
  modelc <- MCMCglmm(fixed=fixed,
                     random = ~ us(trait):animal,
                     rcov = ~ us(trait):units,
                     data = BootsDat,
                     family = rep("gaussian", Nb),
                     prior=priorCom,
                     pedigree=flytree,
                     pr=TRUE,
                     verbose = FALSE,
                     thin=20)
  
  DICCom[k,] = modelc$DIC
  
  MCMCSolCom=modelc$Sol
  SamplesSim = SimulatedSamplesSp(MCMCSolCom)
  
  
  #SINGLE BEHAVIOR MODEL
  
  for(b in 1:Nb){
   print(b)
   Beh = b
   # ---Defining prior (from mulTree, github, and R-sig-ME) nu==n
   phen.var<-var(BootsDat[,Beh+1])
   prior11<-list(G=list(G1=list(V=phen.var/4, n=1.002)), R=list(V=phen.var/4, n=1.002))    
  
   sfix="cbind("
   sfix <-paste(sfix,colnames(BootsDat)[Beh+1],") ~ 1",sep="")
   fixed <- as.formula(sfix)
  
   model0 <- MCMCglmm(fixed = fixed,
                     random = ~animal,
                     rcov = ~units,
                     data = BootsDat,
                     family = "gaussian",
                     prior=prior11,
                     pedigree=flytree,
                     pr=TRUE,
                     verbose = FALSE,
                     thin=20)
  
   #FileName=paste("../../data/SingleTraitResultsFeb2021/Anc134BhModel0PrTrueAllMinusYak00SingleBhFeb2021_",j,sep="")
   #save(model0, file = FileName)
   DICSin[k,b]<-model0$DIC
   
   MCMCSol=model0$Sol
   Samples = SimulatedSamplesSpSingleModel(MCMCSol)
   
   
   JSDSin <- matrix(0, 1, 6)
   JSDCom <- matrix(0, 1, 6)
   #CHECK
   for(l in 1:6) {
     #Samplest = eval(parse(text = SamplesSp[j]))
     Samplest = Samples[[l]] # Simple model samples for behavior b
     Samplestt = SamplesSim[[l]] # Complex model samples for behavior b
     SamplesExpt = eval(parse(text = SamplesExpSp[l]))
     
     #MDist[j] = abs(mean(Samplest)-mean(SamplesExpt))/sqrt(sd(Samplest)^2 + sd(SamplesExpt)^2)
     JSDSin[l] = JSDist(Samplest, SamplesExpt[,b+1])
     JSDCom[l] = JSDist(Samplestt[,b], SamplesExpt[,b+1])
     rm(Samplest,SamplesExpt,Samplestt)
   }
   
   JSDSinSum[k,b] = sum(JSDSin)
   JSDComSum[k,b] = sum(JSDCom)
   
   
  }
  print(k)
}
#}
#write.table(DICCh, file="../../data/SingleTraitResultsFeb2021/DICAncModelallMinusYak00AllSingleBhFeb2021.txt", row.names=FALSE, col.names=FALSE)
#write.table(model0$VCV, file="ShortScriptVCVMCMC.txt", row.names=FALSE, col.names=FALSE)
# attributes(model0$VCV) let us see the names of the elements in the covariance matrices





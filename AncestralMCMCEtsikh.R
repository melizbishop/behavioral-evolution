# ===Ancestral behaviors of melanogaster fly===
# =============================================

# ---Load the necessary libraries


library(coda,lib.loc="Rpackages")
library(ape,lib.loc="Rpackages")
library(Matrix,lib.loc="Rpackages")
library(MCMCglmm,lib.loc="Rpackages")
library(parallel)
#library(MCMCglmm)



# ---Global parameters
Nb = 134 #number of behaviors

# ---Set the output file
# sink("ancestral0.out")

# ---Defining phylogeny -> pedigree variable for MCMCglmm

tt="((yakuba:0.09,santomea:0.09)oldyasa:0.44,((sechelia:0.15,simulans:0.15,mauritiana:0.15)oldsim:0.14,melanogaster:0.29)oldmela:0.24);"

flytree<-read.tree(text=tt)


# ---Load the dataset: l_i= log p(behavior_i) - log p(behavior_0) 
# 134 behaviors + 0 behavior, 593 flies from 6 species (column-> animal)

LFlydat <- read.table("../../data/alphas.txt", header=TRUE, sep="\t", row.names="id")


# ---Defining prior
# non-informative prior
IJ <- (1/(Nb+1))*(diag(Nb)+matrix(1,Nb,Nb))
prior.1<-list(G=list(G1=list(V=IJ/2,n=Nb)), R=list(V=IJ/2,n=Nb))

# ---Creating fixed effects
sfix="cbind("
for (i in 1:(Nb-1)){sfix <-paste(sfix,colnames(LFlydat)[i+1],",",sep="")}
sfix <-paste(sfix,colnames(LFlydat)[Nb+1],") ~ trait-1",sep="")
fixed <- as.formula(sfix)

#---Run MCMCglmm
#model0 <- MCMCglmm(fixed=fixed,
#                   random = ~us(trait):animal,
#                   rcov = ~us(trait):units,
#                   data = LFlydat,
#                   family = rep("gaussian", Nb),
#                   prior=prior.1,
#                   verbose = FALSE)
#write.table(model0$Sol, file="ShortScriptMeanMCMC.txt", row.names=FALSE, col.names=FALSE)
#write.table(model0$VCV, file="ShortScriptVCVMCMC.txt", row.names=FALSE, col.names=FALSE)

#Running in parallel many chains
set.seed(1)
m6 <- mclapply(1:20, function(i) {
  MCMCglmm(fixed=fixed,
           random = ~us(trait):animal,
           rcov = ~us(trait):units,
           data = LFlydat,
           family = rep("gaussian", Nb),
           prior=prior.1,
           pedigree=flytree,
           verbose = FALSE,
           thin=50)
}, mc.cores=2)

MeansMCMC <- lapply(m6, function(m) m$Sol)
MeansMCMC <- do.call(mcmc.list, MeansMCMC)
MeanGelman=gelman.diag(MeansMCMC)

VarMCMC <- lapply(m6, function(m) m$VCV)
VarMCMC <- do.call(mcmc.list, VarMCMC)
#VarGelman=gelman.diag(VarMCMC)

VarGelman=gelman.diag(VarMCMC,multivariate=FALSE)


write.table(MeanGelman$psrf, file="ShortScriptMeanGelman.txt", row.names=FALSE, col.names=FALSE)
write.table(VarGelman$psrf, file="ShortScriptVarGelman.txt", row.names=FALSE, col.names=FALSE)



detach("package:MCMCglmm")

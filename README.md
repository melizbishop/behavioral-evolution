# behavioral-evolution

# Data 
LogBehaviorData.txt: Data set containing the logarithm of the probability of performing a particular behavior i (P_i) for each of the 561 flies.

watershedMap.mat: Data set containing which behavior corresponds to each particular region in the 2D behavioral representation.

Note: Dynamic data showing the sequence of behaviors performed by a fly along time cannot be uploaded due to its size, but can be obtained upon request at catarivera8@gmail.com

# Algorithms

AncestralGLMM.R: Code to fit a Generalized Linear Mixed Model (GLMM) to the behavioral data. 

AncestralMCMC.R and GelmanRubinTest.R: Generate different MCMC chains corresponding to different fitted GLMM and evaluate their convergence using Gelman Rubin test.

InformationBasedClustering: This folder contains the files necessary to apply Information Based Clustering given a distance matrix and obtained the corresponding 
Pareto front. First, use information_cluster_pairwise_distances.m and then use findParetoFrontInfoBasedClust.m.

BottleneckAlgorithmClustering: This folder contains the files necessary to apply the determinist Bottleneck Algorithm to cluster behaviors. First, use run_DIB.m and then, estimate the Pareto front using findParetoFront.m 

Code for Hernandez, Rivera, et al (2021)

# behavioral-evolution

# Data 
LogBehaviorData.txt: Data set containing the logarithm of the probability of performing a particular behavior i (P_i) for each of the 561 flies.

watershedMap.mat: Data set containing which behavior corresponds to each particular region in the 2D behavioral representation.

Note: Dynamic data showing the sequence of behaviors perfomed by a fly along time cannot be uploaded due to its size, but can be obtained upon request.

# Algorithms

AncestralGLMM.R: Code to fit a Generalized Linear Mixed Model (GLMM) to the behavioral data. 

AncestralMCMC.R and GelmanRubinTest.R: Generate different MCMC chains correspoding to the fitted GLMM and evaluate their convergence using Gelman Rubin test.



Code for Hernandez, Rivera, et al (2021)

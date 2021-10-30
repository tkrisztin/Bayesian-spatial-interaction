# Bayesian-spatial-interaction
Bayesian estimator for the Poisson or Negative Binomial spatial interaction model with optional spatial effects. Based on the mixture MCMC methods from the bayesf MATLAB package.

The file [mixturemcmc.R](mixturemcmc.R) is a direct implementation of the estimator in Frühwirth-Schnatter S. et. al (2009) Improved auxiliary mixture sampling for hierarchical models of non-Gaussian data. *Statistics and Computing*. **19**(4), p 479-492 https://link.springer.com/article/10.1007/s11222-008-9109-4

The file [mixturemcmc_fx_fdi.R](mixturemcmc_fx_fdi.R) estimates also origin- and destionation-specific spatial effects. For a detailed description of the algorithm, see the Krisztin T. and Piribauer P. (2021) Modelling European regional FDI flows using a Bayesian spatial Poisson interaction model. *The Annals of Regional Science*, **67**, p 593–616 http://link.springer.com/article/10.1007/s00168-021-01058-x

library(MASS)
library(cubature)

# global parameters: dimension
dim <- 3
source("src/references/cornerPeakFamily.R") # genz function to test

# Bayesian Quadrature method
# set number of new query points using sequential design
source("src/BARTBQ.R")
predictionBART <- mainBARTBQ()
meanValueBART <- predictionBART[1]
standardDeviationBART <- predictionBART

# Bayesian Quadrature with Monte Carlo integration method
source("src/monteCarloIntegration.R")
predictionMonteCarlo <- monteCarloIntegrationUniform(copeak)
meanValueMonteCarlo <- predictionMonteCarlo[1]
standardDeviationMonteCarlo <- predictionMonteCarlo[2]

# Bayesian Quadrature with Gaussian Process
source("src/GPBQ.R")
predictionGPBQ <- computeGPBQ(dim, epochs = 400, N=10, FUN = copeak)
meanValueGP <- predictionGPBQ[1]
standardDeviationGP <- predictionGPBQ[2]


# Exact integral of genz function in hypercube [0,1]^dim
real <- adaptIntegrate(copeak,lowerLimit = rep(0,dim), upperLimit = rep(1,dim))

# Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
percentageErrorBART <- meanValueBART - real[[1]]
percentageErrorMonteCarlo <- meanValueMonteCarlo - real[[1]]
percentageErrorGP <- meanValueGP - real[[1]]


#[5, 4, 3, 2, 1.5, 1.3, 1.2, 1.1]
# mean vs n
#log sd vs log n


# !/usr/bin/env R
# Load required packages
library(MASS)
library(cubature)
library(lhs)
library(data.tree)
library(dbarts)
library(matrixStats)
library(mvtnorm)
library(doParallel)
library(kernlab)

# define string formatting
`%--%` <- function(x, y) 
# from stack exchange:
# https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

# global parameters: dimension
dim <- 1
num_iterations <- 1
whichGenz <- 4
whichKernel <- "matern32"
# turn on/off sequential design
# 1 denotes TRUE to sequential
# 0 denotes FALSE to sequential
cat("\nBegin testing:\n")
sequential = TRUE
cat("Sequantial design set to", sequential, "\n")
# prior measure over the inputs
# uniform by default
measure = "gaussian"
cat("Prior measure:", measure, "\n")

print(c(dim, num_iterations, whichGenz))
source("src/genz/genz.R") # genz function to test
if (whichGenz == 4) { genz <- gaussian_weighted; genzFunctionName <-  deparse(substitute(gaussian)) }
print("Testing with: %s" %--% genzFunctionName)

# prepare training dataset
trainX <- replicate(dim, rtnorm(50, mean=0.5, lower=0, upper=1))
trainY <- genz(trainX)

for (num_cv in 1:5) {
  # Bayesian Quadrature method
  # set number of new query points using sequential design
  source("src/BARTBQ.R")
  t0 <- proc.time()
  predictionBART <- mainBARTBQ(dim, num_iterations, FUN = genz, trainX, trainY, sequential, measure)
  t1 <- proc.time()
  bartTime <- (t1 - t0)[[1]]
  
  # measure2 <- "gaussian"
  # genz2 <- gaussian_weighted
  # trainX <- replicate(dim, rtnorm(100, mean=0.5, lower=0, upper=1, sd=dim^2))
  # trainY <- genz2(trainX)
  # predictionBART2 <- mainBARTBQ(dim, num_iterations, FUN = genz2, trainX, trainY, sequential, measure2)
  
  # Bayesian Quadrature with Monte Carlo integration method
  print("Begin Monte Carlo Integration")
  source("src/monteCarloIntegration.R")
  
  t0 <- proc.time()
  predictionMonteCarlo <- monteCarloIntegrationUniform(FUN = genz, numSamples = num_iterations, dim, measure)
  t1 <- proc.time()
  
  MITime <- (t1 - t0)[[1]]
  
  # Bayesian Quadrature with Gaussian Process
  print("Begin Gaussian Process Integration")
  if (num_cv == 1) {
  library(reticulate)
  source("src/optimise_gp.R")
  lengthscale <- optimise_gp_r(trainX, trainY, kernel = whichKernel, epochs=500)
  }
  
  source("src/GPBQ.R")
  t0 <- proc.time()
  # need to add in function to optimise the hyperparameters
  predictionGPBQ <- computeGPBQ(
    trainX, 
    trainY, 
    dim, 
    epochs = num_iterations-1, 
    kernel = whichKernel, 
    FUN = genz, 
    lengthscale,
    sequential, 
    measure
  )  
  t1 <- proc.time()
  GPTime <- (t1 - t0)[[1]]
  
  # Read in analytical integrals
  dimensionsList <- c(1,2,3,5,10,20)
  whichDimension <- which(dim == dimensionsList)
  if (whichGenz <= 6){
    analyticalIntegrals <- read.csv("results/genz/integrals.csv", header = FALSE)
    real <- analyticalIntegrals[whichGenz, whichDimension]
  } else if (whichGenz == 7) {
    source("src/genz/analyticalIntegrals.R")
    real <- stepIntegral(dim, jumps)
  } else {
    if (whichGenz == 8 & dim ==1){ real <- 0.008327796}
    if (whichGenz == 8 & dim ==2){ real <- 0.008327796 * 2}
    if (whichGenz == 8 & dim ==3){ real <- 0.008327796 * 3}
  }
  
  
  # Read in analytical integrals
  dimensionsList <- c(1,2,3,5,10,20)
  whichDimension <- which(dim == dimensionsList)
  analyticalIntegrals <- read.csv("results/genz/integrals.csv", header = FALSE)
  real <- analyticalIntegrals[whichGenz, whichDimension]
  
  # Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
  print("Final Results:")
  print(c("Actual integral:", real))
  print(c("BART integral:", predictionBART$meanValueBART[num_iterations]))
  print(c("MI integral:", predictionMonteCarlo$meanValueMonteCarlo[num_iterations]))
  print(c("GP integral:", predictionGPBQ$meanValueGP[num_iterations]))
  
  print("Writing full results to results/genz%s" %--% c(whichGenz))
  results <- data.frame(
    "epochs" = c(1:num_iterations),
    "BARTMean" = predictionBART$meanValueBART, "BARTsd" = predictionBART$standardDeviationBART,
    "MIMean" = predictionMonteCarlo$meanValueMonteCarlo, "MIsd" = predictionMonteCarlo$standardDeviationMonteCarlo,
    "GPMean" = predictionGPBQ$meanValueGP, "GPsd" = sqrt(predictionGPBQ$varianceGP),
    "actual" = rep(real, num_iterations),
    "runtimeBART" = rep(bartTime, num_iterations),
    "runtimeMI" = rep(MITime, num_iterations),
    "runtimeGP" = rep(GPTime, num_iterations)
  )
  print("Begin Plots")
  if (!sequential){
      csvName <- "results/genz/%s/%sDim%sNoSequential%s_%s.csv" %--% c(
      whichGenz, 
      genzFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
  } else {
    csvName <- "results/genz/%s/%sDim%s%s_%s.csv" %--% c(
      whichGenz, 
      genzFunctionName,
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
  }
  write.csv(results, file = csvName, row.names=FALSE)
}

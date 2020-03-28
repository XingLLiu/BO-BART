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
library(msm)
library(MCMCglmm)
set.seed(0)

# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

# global parameters: dimension
args <- commandArgs(TRUE)
dim <- 1
num_iterations <- 1
n <- 20
whichGenz <- 7
whichKernel <- "matern32"
# turn on/off sequential design
# 1 denotes TRUE to sequential
# 0 denotes FALSE to sequential
sequential <- FALSE
cat("Sequantial design set to", sequential, "\n")
# prior measure over the inputs
# uniform by default
measure <- "uniform"
cat("Prior measure:", measure, "\n")
jumps <- 1
cat("Number of jumps for step function:", jumps, "\n")


print(c(dim, num_iterations, whichGenz))
source("src/genz/genz.R") # genz function to test

if (whichGenz == 7) {genz <- function(xx){return(step(xx, jumps=jumps))}; genzFunctionName <-  deparse(substitute(step)) }

print("Testing with: %s" %--% genzFunctionName)

# prepare training dataset
trainX <- replicate(dim, runif(n))
trainY <- genz(trainX)
source("src/BARTBQ.R")
t0 <- proc.time()
predictionBART <- mainBARTBQ(dim, num_iterations, FUN = genz, trainX, trainY, sequential, measure)
t1 <- proc.time()
bartTime <- (t1 - t0)[[1]]

# Bayesian Quadrature with Monte Carlo integration method
print("Begin Monte Carlo Integration")
source("src/monteCarloIntegration.R")

t0 <- proc.time()
predictionMonteCarlo <- monteCarloIntegrationUniform(FUN = genz, trainX, trainY, numSamples=num_iterations, dim, measure)
t1 <- proc.time()
  
MITime <- (t1 - t0)[[1]]
  
# Bayesian Quadrature with Gaussian Process
lengthscale <- 1
source("src/GPBQ.R")
t0 <- proc.time()
lengthscale=1
# need to add in function to optimise the hyperparameters
predictionGPBQ <- computeGPBQ(
  trainX, 
  trainY, 
  dim, 
  epochs = num_iterations, 
  kernel = whichKernel, 
  FUN = genz, 
  lengthscale,
  sequential, 
  measure
)  
t1 <- proc.time()
GPTime <- (t1 - t0)[[1]]

print("Writing full results to results/genz%s" %--% c(whichGenz))
results <- data.frame(
  "epochs" = c(1:num_iterations),
  "runtimeBART" = rep(bartTime, num_iterations),
  "runtimeMI" = rep(MITime, num_iterations),
  "runtimeGP" = rep(GPTime, num_iterations)
)
csvName <- "results/genz/%s/computational_complexity_%sDim%sNoSequential%s_num%s.csv" %--% c(
  whichGenz, 
  genzFunctionName,
  dim,
  tools::toTitleCase(measure),
  n
)
write.csv(results, file = csvName, row.names=FALSE)


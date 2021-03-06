# !/usr/bin/env R
# uncomment below and fix it according when in department cluster.
setwd(getwd())

# Load required packages
library(lhs)
library(mvtnorm)
library(matrixStats)
library(kernlab)
library(cubature)
set.seed(0)

# define string formatting
`%--%` <- function(x, y) 
# from stack exchange:
# https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

calc_lengthscale <- function(dim) {
  samples <- c()
  for (i in 1:dim) {
    samples <- cbind(samples, runif(1000))
  }
  lengthscale <- median(dist(samples))
}

# global parameters: dimension
args <- as.double(commandArgs(TRUE))
dim <- args[1]
num_iterations <- args[2]
whichGenz <- args[3]
num_training <- args[4]
lengthscale <- calc_lengthscale(dim)

if (num_iterations == 1) { stop("NEED MORE THAN 1 ITERATION") }

print(c(dim, num_iterations, whichGenz))
source("./genz/genz.R") # genz function to test

if (whichGenz < 1 | whichGenz > 6) { stop("undefined genz function. Change 3rd argument to 1-6") }
if (whichGenz == 3 & dim == 1) { stop("incorrect dimension. Discrete Genz function only defined for dimension >= 2") } 

if (whichGenz == 1) { genz <- cont; genzFunctionName <-  deparse(substitute(cont)) }
if (whichGenz == 2) { genz <- copeak; genzFunctionName <-  deparse(substitute(copeak)) }
if (whichGenz == 3) { genz <- disc; genzFunctionName <-  deparse(substitute(disc)) }
if (whichGenz == 4) { genz <- gaussian; genzFunctionName <-  deparse(substitute(gaussian)) }
if (whichGenz == 5) { genz <- oscil; genzFunctionName <-  deparse(substitute(oscil)) }
if (whichGenz == 6) { genz <- prpeak; genzFunctionName <-  deparse(substitute(prpeak)) }

print("Testing with: %s" %--% genzFunctionName)

# Bayesian Quadrature with Gaussian Process
print("Begin Gaussian Process Integration")
source("./GPBQ.R")

t0 <- proc.time()
predictionGPBQ <- computeGPBQ(dim, epochs = num_iterations-1, N=num_training, FUN = genz, lengthscale)  
t1 <- proc.time()
GPTime <- (t1 - t0)[[1]]

# read in analytical integrals
dimensionsList <- c(1,2,3,5,10,20)
whichDimension <- which(dim == dimensionsList)
analyticalIntegrals <- read.csv("./genz/integrals.csv", header = FALSE)
real <- analyticalIntegrals[whichGenz, whichDimension]

# Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
print("Final Results:")
print(c("Actual integral:", real))
print(c("GP integral:", predictionGPBQ$meanValueGP[num_iterations]))
print(c("lenthscale=", lengthscale))

print("Writing full results to ../results/genz%s_dim%s" %--% c(whichGenz, dim))
results <- data.frame(
		"epochs" = c(1:num_iterations),
		"GPMean" = predictionGPBQ$meanValueGP,
		"GPVariance" = predictionGPBQ$varianceGP,
        "runtimeGP" = GPTime
)

write.csv(results, file = "../results/genz/GPresults%sdim%s_training%s.csv" %--% c(whichGenz, dim, num_training),row.names=FALSE)

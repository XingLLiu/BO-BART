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
args <- commandArgs(TRUE)
dim <- as.double(args[1])
num_iterations <- as.double(args[2])
whichGenz <- as.double(args[3])
whichKernel <- as.character(args[4])
# turn on/off sequential design
# 1 denotes TRUE to sequential
# 0 denotes FALSE to sequential
cat("\nBegin testing:\n")
if (as.double(args[5]) == 1 | is.na(as.double(args[5]))) {
  sequential <- TRUE
} else {
  sequential <- FALSE
}
cat("Sequantial design set to", sequential, "\n")
# prior measure over the inputs
# uniform by default
if (as.character(args[6]) != "gaussian" | is.na(args[6])) {
  measure <- "uniform"
} else{
  measure <- as.character(args[6])
}
cat("Prior measure:", measure, "\n")
# extra parameter for step function
# 1 by default
if (whichGenz == 7 & is.na(args[7])) {
  jumps <- 1
  cat("Number of jumps for step function:", jumps, "\n")
} else if (whichGenz == 7){
  jumps <- as.double(args[7])
  cat("Number of jumps for step function:", jumps, "\n")
}

if (num_iterations == 1) { stop("NEED MORE THAN 1 ITERATION") }

print(c(dim, num_iterations, whichGenz))
source("src/genz/genz.R") # genz function to test

if (whichGenz < 1 | whichGenz > 8) { stop("undefined genz function. Change 3rd argument to 1-8") }
if (whichGenz == 3 & dim == 1) { stop("incorrect dimension. Discrete Genz function only defined for dimension >= 2") } 

if (whichGenz == 1) { genz <- cont; genzFunctionName <-  deparse(substitute(cont)) }
if (whichGenz == 2) { genz <- copeak; genzFunctionName <-  deparse(substitute(copeak)) }
if (whichGenz == 3) { genz <- disc; genzFunctionName <-  deparse(substitute(disc)) }
if (whichGenz == 4) { genz <- gaussian; genzFunctionName <-  deparse(substitute(gaussian)) }
if (whichGenz == 5) { genz <- oscil; genzFunctionName <-  deparse(substitute(oscil)) }
if (whichGenz == 6) { genz <- prpeak; genzFunctionName <-  deparse(substitute(prpeak)) }
if (whichGenz == 7) { genz <- function(xx){return(step(xx, jumps=jumps))}; genzFunctionName <-  deparse(substitute(step)) }
if (whichGenz == 8) { genz <- mix; genzFunctionName <-  deparse(substitute(mix)) }

print("Testing with: %s" %--% genzFunctionName)

# prepare training dataset
if (measure == "uniform") {
  trainX <- replicate(dim, runif(100))
} else if (measure == "gaussian") {
  trainX <- replicate(dim, rtnorm(100, mean = 0.5, lower=0, upper=1, dim^2/100))
}
trainY <- genz(trainX)

# Bayesian Quadrature method
# set number of new query points using sequential design
source("src/BARTBQ.R")
t0 <- proc.time()
predictionBART <- mainBARTBQ(dim, num_iterations, FUN = genz, trainX, trainY, sequential, measure)
t1 <- proc.time()
bartTime <- (t1 - t0)[[1]]

measure2 <- "gaussian"
trainX <- replicate(dim, rtnorm(100, mean = 0.5, lower=0, upper=1))
trainY <- genz(trainX)
predictionBART2 <- mainBARTBQ(dim, num_iterations, FUN = genz, trainX, trainY, sequential, measure2)

# Bayesian Quadrature with Monte Carlo integration method
print("Begin Monte Carlo Integration")
source("src/monteCarloIntegration.R")

t0 <- proc.time()
predictionMonteCarlo <- monteCarloIntegrationUniform(FUN = genz, numSamples=num_iterations, dim)
t1 <- proc.time()

MITime <- (t1 - t0)[[1]]

print("Begin Plots")
figName <- "Figures/GaussianUniformPriors%sDim%s%s.pdf" %--% c(
          genzFunctionName,
          dim,
          tools::toTitleCase(measure)
          )

# 1. Open jpeg file
pdf(figName, width = 10, height = 11)

# 2. Create the plot
par(mfrow = c(1,1), pty = "s")
plot(predictionMonteCarlo$meanValueMonteCarlo, ylim = c(0, 0.1), main = genzFunctionName)
points(predictionBART$meanValueBART, col = "red")
points(predictionBART2$meanValueBART, col = "blue")
legend("topright", legend = c("MC", "BART uniform", "BART Gaussian"), col = c("black", "red", "blue"), cex=0.8, lty = c(1,1,1))

# 3. Close the file
invisible(dev.off())

print("Please check {ROOT}/%s for plots" %--% figName)



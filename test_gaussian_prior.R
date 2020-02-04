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
# number of jumps if step function is used
jumps <- 1

print(c(dim, num_iterations, whichGenz))
source("src/genz/genz.R") # genz function to test

if (whichGenz == 4) { genz <- gaussian_weighted; genzFunctionName <-  deparse(substitute(gaussian)) }
if (whichGenz == 6) { genz <- prpeak_weighted; genzFunctionName <-  deparse(substitute(prpeak)) }
if (whichGenz == 7) { genz <- function(xx){return(step(xx, jumps=jumps))}; genzFunctionName <-  deparse(substitute(step)) }
if (whichGenz == 8) { genz <- function(xx){return(step_weighted(xx, jumps=jumps))}; genzFunctionName <-  deparse(substitute(step)) }
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
    figName <- "Figures/%s/%sDim%sNoSequential%s_%s.pdf" %--% c(
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
    figName <- "Figures/%s/%sDim%s%s_%s_.pdf" %--% c(
            whichGenz,
            genzFunctionName,
            dim,
            tools::toTitleCase(measure),
            num_cv
    )
  }

  write.csv(results, file = csvName, row.names=FALSE)

  print("Begin Plots")
  # 1. Open jpeg file
  pdf(figName, width = 10, height = 11)
  # 2. Create the plot
  par(mfrow = c(1,2), pty = "s")
  plot(x = c(1:num_iterations), y = predictionMonteCarlo$meanValueMonteCarlo,
       pch = 16, type = "l",
       xlab = "Number of epochs N", ylab = "Integral approximation", col = "blue",
       main = "Convergence of methods: mean vs N \nusing %s with %s epochs in %s dim" %--% c(genzFunctionName, num_iterations, dim),
       ylim = c(-real, 3 * real), 
       lty = 1,
       xaxs="i", yaxs="i"
       )
  lines(x = c(1:num_iterations), predictionBART$meanValueBART, type = 'l', col = "red", lty = 1)
  lines(x = c(1:num_iterations), predictionGPBQ$meanValueGP, type = 'l', col = "green", lty = 1)
  abline(a = real, b = 0, lty = 4)
  legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ", "Actual"),
         col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1,1))
  
  # 2. Create the plot
  plot(x = log(c(2:num_iterations)), y = log(predictionMonteCarlo$standardDeviationMonteCarlo[-1]),
       pch = 16, type = "l",
       xlab = "Number of epochs N", ylab = "Log standard deviation", col = "blue",
       main = "Convergence of methods: log(sd) vs log(N) \nusing %s with %s epochs in %s dim" %--% c(genzFunctionName, num_iterations, dim),
       lty = 1,
       xaxs="i", yaxs="i")
  lines(x = log(c(2:num_iterations)), log(predictionBART$standardDeviationBART[-1]), type = 'l', col = "red", lty = 1)
  lines(x = log(c(2:num_iterations)), log(sqrt(predictionGPBQ$varianceGP[-1])), type = 'l', col = "green", lty = 1)
  legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ"),
         col=c("blue", "red", "green"), cex=0.8, lty = c(1,1,1,1))
  # 3. Close the file
  invisible(dev.off())
  
  print("Please check {ROOT}/%s for plots" %--% figName)
}

# --------------------------------
# To run: input two command-line arguements
# Rscript test_step_func.R dim num_of_iters sequential
# --------------------------------

# !/usr/bin/env R
# setwd(getwd())
# uncomment the following when running the code for the first time to load real integral values
# source("./genz/saveComputeIntegrals.R")

# Load required packages
# source("src/packages/requiredPackages.R")
# requiredPackages()
library(yaml)
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
jump <- 2
whichGenz <- 7
whichKernel <- as.character(args[4])

# turn on/off sequential design
# 1 denotes TRUE to sequential
# 0 denotes FALSE to sequential
cat("\nBegin testing:\n")
if (as.double(args[3]) == 1 | is.na(as.double(args[3]))) {
  sequential <- TRUE
  print("Sequantial design set to TRUE.")
} else {
  sequential <- FALSE
  print("Sequantial design set to FALSE.")
}

if (num_iterations == 1) { stop("NEED MORE THAN 1 ITERATION") }

print(c(dim, num_iterations, whichGenz, sequential, whichKernel))
source("src/genz/genz.R") # genz function to test

if (whichGenz == 7) { genz <- step; genzFunctionName <-  deparse(substitute(step)) }

print("Testing with: %s" %--% genzFunctionName)

# prepare training dataset
trainX <- replicate(dim, runif(100))
trainY <- genz(trainX, jump = jump)


# Bayesian Quadrature method
# set number of new query points using sequential design

source("src/BARTBQ.R")
t0 <- proc.time()
predictionBART <- mainBARTBQ(dim, num_iterations, FUN = genz, trainX, trainY, sequential)
t1 <- proc.time()
bartTime <- (t1 - t0)[[1]]

# Bayesian Quadrature with Monte Carlo integration method
print("Begin Monte Carlo Integration")
source("src/monteCarloIntegration.R")
t0 <- proc.time()
predictionMonteCarlo <- monteCarloIntegrationUniform(FUN = genz, numSamples=num_iterations, dim)
t1 <- proc.time()

MITime <- (t1 - t0)[[1]]

# Bayesian Quadrature with Gaussian Process
print("Begin Gaussian Process Integration")
library(reticulate)
source("src/optimise_gp.R")
lengthscale <- optimise_gp_r(trainX, trainY, kernel = whichKernel, epochs=500)
source("src/GPBQ.R")
t0 <- proc.time()
# need to add in function to optimise the hyperparameters
predictionGPBQ <- computeGPBQ(trainX, trainY, dim, epochs = num_iterations-1, kernel = whichKernel, FUN = genz, lengthscale,sequential)  
t1 <- proc.time()
GPTime <- (t1 - t0)[[1]]

dimensionsList <- c(1,2,3,5,10,20)
whichDimension <- which(dim == dimensionsList)
real <- mean(predictionMonteCarlo$meanValueMonteCarlo)
############
plot(predictionGPBQ$meanValueGP, ylim = c(-0.01, 1.01))
############

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

if (!sequential){
  csvName <- "results/genz/%s/results%sdim%sNoSequential.csv" %--% c(
          whichGenz, 
          whichGenz, 
          dim
     )
  figName <- "Figures/%s/%s%sDimNoSequential.pdf" %--% c(whichGenz, genzFunctionName, dim)
} else {
  csvName <- "results/genz/%s/results%sdim%s.csv" %--% c(
          whichGenz, 
          whichGenz, 
          dim
     )
  figName <- "Figures/%s/%s%sDim.pdf" %--% c(whichGenz, genzFunctionName, dim)
}

write.csv(results, file = csvName, row.names=FALSE)

# 2. Create the plot
pdf(figName, width = 10, height = 11)
par(mfrow = c(1,2), pty = "s")
plot(x = c(1:num_iterations), y = predictionMonteCarlo$meanValueMonteCarlo,
     pch = 16, type = "l",
     xlab = "Number of epochs N", ylab = "Integral approximation", col = "blue",
     main = "Convergence of methods: mean vs N \nusing %s with %s epochs in %s dim" %--% c(genzFunctionName, num_iterations, dim),
     ylim = c(-real, 1 + real), 
     lty = 1,
     xaxs="i", yaxs="i"
     )
lines(x = c(1:num_iterations), predictionBART$meanValueBART, type = 'l', col = "red", lty = 1)
lines(x = c(1:num_iterations), predictionGPBQ$meanValueGP, type = 'l', col = "green", lty = 1)
abline(a = real, b = 0, lty = 4)
legend("bottomleft", legend=c("MC Integration", "BART BQ", "GP BQ", "Actual"),
       col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1,1))

# 2. Create the plot
plot(x = log(c(2:num_iterations)), y = log(predictionMonteCarlo$standardDeviationMonteCarlo[-1]),
     pch = 16, type = "l",
     xlab = "Number of epochs N", ylab = "Log standard deviation", col = "blue",
     main = "Convergence of methods: log(sd) vs log(N) \nusing %s with %s epochs in %s dim" %--% c(genzFunctionName, num_iterations, dim),
     ylim = c(-5, 2),
     lty = 1,
     xaxs="i", yaxs="i")
lines(x = log(c(2:num_iterations)), log(predictionBART$standardDeviationBART[-1]), type = 'l', col = "red", lty = 1)
lines(x = log(c(2:num_iterations)), log(sqrt(predictionGPBQ$varianceGP[-1])), type = 'l', col = "green", lty = 1)
legend("bottomleft", legend=c("MC Integration", "BART BQ", "GP BQ"),
       col=c("blue", "red", "green"), cex=0.8, lty = c(1,1,1,1))

# 3. Close the file
invisible(dev.off())

print("Please check {ROOT}/%s for plots" %--% figName)

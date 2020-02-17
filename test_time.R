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
dim <- as.double(args[1])
num_data <- 50*dim
num_iterations <- 1
whichKernel <- "matern32"
sequential <- FALSE
# prior measure over the inputs
# uniform by default
measure <- "uniform"
print(c(dim, num_iterations))

source("src/genz/fisher_integrands.R")
# C <- rep(0.45958094818045914, dim)
# R <- rep(0.13332051229486602, dim)
# H <- rep(1.5020584867022158, dim)
# F <- rep(4.624266954512351, dim)
# P <- rep(1, dim)

C <- replicate(1, runif(100, 0.1, 0.9))
R <- replicate(1, rbeta(100, 5, 2))
H <- replicate(1, runif(100, 0.5*exp(1), 1.5*exp(1)))
F <- replicate(1, runif(100, 0, 5))
P <- replicate(1, rbinom(100, 1, 0.5))
# C <- c(0.1706093, 0.5319923, 0.7117816)
# R <- c(0.6786221, 0.7207544, 0.5120249)
# H <- c(3.029867, 3.065427, 3.114357)
# F <- c(4.607562, 4.526243, 2.221768)
# P <- c(1, 0, 1)

cut_point <- 0.5
fisher_function_full <- create_fisher_function(C, R, H, F, P, 100)

fisher_function <- function(x) {
  x_in <- cbind(x, matrix(cut_point, nrow = nrow(x), ncol = 100-dim))
  return(fisher_function_full(x_in))
}


# prepare training dataset
if (measure == "uniform") {
  trainX<- replicate(dim, runif(num_data, 0, 1))
  trainY <- fisher_function(trainX)
} else if (measure == "gaussian") {
  trainX <- replicate(dim, rtnorm(num_data, mean=0.5, lower=0, upper=1))
  trainY <- fisher_function(trainX)
}

real <- 1
for (i in 1:dim) {
  fisher_1d <- create_fisher_function(C[i], R[i], H[i], F[i], P[i], 1)
  real <- real*estimate_real_integral(fisher_1d, 1, 1e7)
}
print(real)
if (dim != 100) {
  for (i in (dim+1):100) {
    fisher_1d <- create_fisher_function(C[i], R[i], H[i], F[i], P[i], 1)
    print(fisher_1d(matrix(cut_point)))
    real <- real * fisher_1d(matrix(cut_point))
  }
}
  

source("src/BARTBQ.R")
t0 <- proc.time()
predictionBART <- mainBARTBQ(dim, num_iterations, FUN = fisher_function, trainX, trainY, sequential, measure)
t1 <- proc.time()
bartTime <- (t1 - t0)[[1]]

# Bayesian Quadrature with Monte Carlo integration method
print("Begin Monte Carlo Integration")
source("src/monteCarloIntegration.R")

t0 <- proc.time()
predictionMonteCarlo <- monteCarloIntegrationUniform(FUN = fisher_function, trainX, trainY, numSamples=num_iterations, dim, measure)
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
predictionGPBQ <- computeGPBQ(
trainX, 
trainY, 
dim, 
epochs = num_iterations, 
kernel = whichKernel, 
FUN = fisher_function, 
lengthscale,
sequential, 
measure
)  
t1 <- proc.time()
GPTime <- (t1 - t0)[[1]]

# Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
print("Final Results:")
print(c("Actual integral:", real))
print(c("BART integral:", predictionBART$meanValueBART[num_iterations]))
print(c("MI integral:", predictionMonteCarlo$meanValueMonteCarlo[num_iterations]))
print(c("GP integral:", predictionGPBQ$meanValueGP[num_iterations]))

print("Writing full results to results/fisher_function")
cat(length(predictionBART$meanValueBART), length(predictionGPBQ$meanValueGP), length(predictionMonteCarlo$meanValueMonteCarlo))
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
csvName <- "results/fisher_function/TimeDim%sNoSequential%s.csv" %--% c(
    dim,
    tools::toTitleCase(measure)
)
} else {
csvName <- "results/fisher_function/TimeDim%s%s.csv" %--% c(
    dim,
    tools::toTitleCase(measure)
)
}

results_models <- list("BART"=predictionBART, "GP"=predictionGPBQ, "MC"=predictionMonteCarlo)
save(results_models, file = "results/fisher_function/TimeDim%s%s.RData" %--% c(
dim,
tools::toTitleCase(measure)
))

write.csv(results, file = csvName, row.names=FALSE)

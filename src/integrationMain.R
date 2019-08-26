# !/usr/bin/env R
setwd(getwd())
# uncomment the following when running the code for the first time to load real integral values
# source("./genz/saveComputeIntegrals.R")

# Load required packages
source("./packages/requiredPackages.R")
requiredPackages()

# define string formatting
`%--%` <- function(x, y) 
# from stack exchange:
# https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

# global parameters: dimension
args <- as.double(commandArgs(TRUE))
dim <- args[1]
num_iterations <- args[2]
whichGenz <- args[3]

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

# prepare training dataset
trainX <- randomLHS(100, dim)
trainY <- genz(trainX)


# Bayesian Quadrature method
# set number of new query points using sequential design

source("./BARTBQ.R")
t0 <- proc.time()
predictionBART <- mainBARTBQ(dim, num_iterations, FUN = genz, trainX, trainY)
t1 <- proc.time()
bartTime <- (t1 - t0)[[1]]
# Bayesian Quadrature with Monte Carlo integration method
print("Begin Monte Carlo Integration")
source("./monteCarloIntegration.R")
t0 <- proc.time()
predictionMonteCarlo <- monteCarloIntegrationUniform(FUN = genz, numSamples=num_iterations, dim)
t1 <- proc.time()

MITime <- (t1 - t0)[[1]]

# Bayesian Quadrature with Gaussian Process
print("Begin Gaussian Process Integration")
source("./GPBQ.R")

t0 <- proc.time()
predictionGPBQ <- computeGPBQ(dim, epochs = num_iterations-1, N=10, FUN = genz)  
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
print(c("BART integral:", predictionBART$meanValueBART[num_iterations]))
print(c("MI integral:", predictionMonteCarlo$meanValueMonteCarlo[num_iterations]))
print(c("GP integral:", predictionGPBQ$meanValueGP[num_iterations]))

print("Writing full results to ../results/genz%s" %--% c(whichGenz))
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
write.csv(results, file = "../results/genz/%s/results%sdim%s.csv" %--% c(whichGenz, whichGenz, dim),row.names=FALSE)

print("Begin Plots")

# 1. Open jpeg file
jpeg("../report/Figures/%s/convergenceMean%s%sDimensions.jpg" %--% c(whichGenz, genzFunctionName, dim), width = 700, height = 583)
# 2. Create the plot
par(mfrow = c(1,2), pty = "s")
plot(x = c(1:num_iterations), y = predictionMonteCarlo$meanValueMonteCarlo,
     pch = 16, type = "l",
     xlab = "Number of epochs N", ylab = "Integral approximation", col = "blue",
     main = "Convergence of methods: mean vs N \nusing %s with %s epochs in %s dim" %--% c(genzFunctionName, num_iterations, dim),
     ylim = c(-real, real + real), 
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
dev.off()

print("Please check {ROOT}/report/Figures for plots")



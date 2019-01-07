library(MASS)
library(cubature)

#define string formatting
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

source("./references/genz.R") # genz function to test
if (whichGenz < 1 | whichGenz > 6) stop("undefined genz function. Change 3rd argument to 1-6") 

if (whichGenz == 1) genz <- cont; genzFunctionName <-  deparse(substitute(cont))
if (whichGenz == 2) genz <- copeak; genzFunctionName <-  deparse(substitute(cont))
if (whichGenz == 3) genz <- prpeak; genzFunctionName <-  deparse(substitute(prpeak))
if (whichGenz == 4) genz <- gaussian; genzFunctionName <-  deparse(substitute(gaussian))
if (whichGenz == 5) genz <- oscil; genzFunctionName <-  deparse(substitute(oscil))
if (whichGenz == 6) genz <- disc; genzFunctionName <-  deparse(substitute(disc))

print("Testing with: %s" %--% genzFunctionName)

# Bayesian Quadrature method
# set number of new query points using sequential design
source("./BARTBQ.R")
predictionBART <- mainBARTBQ(dim, num_iterations, FUN = genz)

# Bayesian Quadrature with Monte Carlo integration method
source("./monteCarloIntegration.R")
predictionMonteCarlo <- monteCarloIntegrationUniform(FUN = genz, numSamples=num_iterations, dim)

# Bayesian Quadrature with Gaussian Process
source("./GPBQ.R")
predictionGPBQ <- computeGPBQ(dim, epochs = num_iterations-1, N=10, FUN = genz)


# Exact integral of genz function in hypercube [0,1]^dim
real <- adaptIntegrate(genz,lowerLimit = rep(0,dim), upperLimit = rep(1,dim))

# Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
percentageErrorBART <- predictionBART$meanValueBART - real[[1]]
percentageErrorMonteCarlo <- predictionMonteCarlo$meanValueMonteCarlo - real[[1]]
percentageErrorGP <- predictionGPBQ$meanValueGP - real[[1]]

print("Begin Plots")

# 1. Open jpeg file
jpeg("../report/Figures/convergence%sStandardDeviation%sDimensions.jpg" %--% c(genzFunctionName, dim), width = 700, height = 583)
# 2. Create the plot
plot(x = log(c(1:num_iterations)), y = log(predictionMonteCarlo$standardDeviationMonteCarlo),
     pch = 16, frame = FALSE, type = "l",
     xlab = "Number of epochs N", ylab = "Log standard deviation", col = "blue",
     main = "Convergence of methods: log(sd) vs log(N) \nusing %s with %s epochs in %s dim" %--% c(genzFunctionName, num_iterations, dim),
     lty = 1)
lines(x = log(c(1:num_iterations)), log(predictionBART$standardDeviationBART), type = 'l', col = "red", lty = 1)
lines(x = log(c(1:num_iterations)), log(predictionGPBQ$standardDeviationGP), type = 'l', col = "green", lty = 1)
legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ"),
       col=c("blue", "red", "green"), cex=0.8, lty = c(1,1,1,1))
# 3. Close the file
dev.off()

# 1. Open jpeg file
jpeg("../report/Figures/convergenceMean%s%sDimensions.jpg" %--% c(genzFunctionName, dim), width = 700, height = 583)
# 2. Create the plot
plot(x = c(1:num_iterations), y = predictionMonteCarlo$meanValueMonteCarlo,
     pch = 16, frame = FALSE, type = "l",
     xlab = "Number of epochs N", ylab = "Integral approximation", col = "blue",
     main = "Convergence of methods: mean vs N \nusing %s with %s epochs" %--% c(genzFunctionName, num_iterations, dim),
     ylim = c(0, real[[1]] + real[[1]]), 
     lty =1
     )
lines(x = c(1:num_iterations), predictionBART$meanValue, type = 'l', col = "red", lty = 1)
lines(x = c(1:num_iterations), predictionGPBQ$meanValueGP, type = 'l', col = "green", lty = 1)
abline(a = real[[1]], b = 0, lty = 4)
legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ", "Actual"),
       col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1,1))
# 3. Close the file
dev.off()

print("Please check {ROOT}/report/Figures for plots")
library(MASS)
library(cubature)

# global parameters: dimension
dim <- 3
source("src/references/cornerPeakFamily.R") # genz function to test

# Bayesian Quadrature method
# set number of new query points using sequential design
source("src/BARTBQ.R")
predictionBART <- mainBARTBQ(dim)

# Bayesian Quadrature with Monte Carlo integration method
source("src/monteCarloIntegration.R")
predictionMonteCarlo <- monteCarloIntegrationUniform(copeak, numSamples=400, dim)

# Bayesian Quadrature with Gaussian Process
source("src/GPBQ.R")
predictionGPBQ <- computeGPBQ(dim, epochs = 399, N=10, FUN = copeak)


# Exact integral of genz function in hypercube [0,1]^dim
real <- adaptIntegrate(copeak,lowerLimit = rep(0,dim), upperLimit = rep(1,dim))

# Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
percentageErrorBART <- meanValueBART - real[[1]]
percentageErrorMonteCarlo <- meanValueMonteCarlo - real[[1]]
percentageErrorGP <- meanValueGP - real[[1]]

# 1. Open jpeg file
jpeg("report/Figures/convergenceCopeakStandardDeviation.jpg", width = 700, height = 583)
# 2. Create the plot
plot(x = log(c(1:400)), y = log(predictionMonteCarlo$standardDeviation2),
     pch = 16, frame = FALSE, type = "l",
     xlab = "Number of epochs N", ylab = "Log standard deviation", col = "blue",
     main = "Convergence of methods: log(sd) vs log(N) using copeak with 400 epochs",
     lty = 1)
lines(x = log(c(1:400)), log(predictionBART$standardDeviation), type = 'l', col = "red", lty = 2)
lines(x = log(c(1:400)), log(predictionGPBQ$standardDeviationGP), type = 'l', col = "green", lty = 3)
legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ"),
       col=c("blue", "red", "green"), lty=1:3, cex=0.8)
# 3. Close the file
dev.off()

# 1. Open jpeg file
jpeg("report/Figures/convergenceCopeak.jpg", width = 700, height = 583)
# 2. Create the plot
plot(x = c(1:400), y = predictionMonteCarlo$meanValue2,
     pch = 16, frame = FALSE, type = "l",
     xlab = "Number of epochs N", ylab = "Integral approximation", col = "blue",
     main = "Convergence of methods: mean vs N using copeak with 400 epochs",
     ylim = c(0, real[[1]] + real[[1]]), lty = 1
     )
lines(x = c(1:400), predictionBART$meanValue, type = 'l', col = "red", lty = 2)
lines(x = c(1:400), predictionGPBQ$meanValueGP, type = 'l', col = "green", lty = 3)
abline(a = real[[1]], b = 0, lty = 4)
legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ", "Actual"),
       col=c("blue", "red", "green", "black"), lty=1:4, cex=0.8)
# 3. Close the file
dev.off()

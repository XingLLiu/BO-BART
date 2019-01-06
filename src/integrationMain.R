library(MASS)
library(cubature)

# global parameters: dimension
dim <- 3
source("src/references/cornerPeakFamily.R") # genz function to test

# Bayesian Quadrature method
# set number of new query points using sequential design
source("src/BARTBQ.R")
predictionBART <- mainBARTBQ()

# Bayesian Quadrature with Monte Carlo integration method
source("src/monteCarloIntegration.R")
predictionMonteCarlo <- monteCarloIntegrationUniform(copeak)

# Bayesian Quadrature with Gaussian Process
source("src/GPBQ.R")
predictionGPBQ <- computeGPBQ(dim, epochs = 400, N=10, FUN = copeak)


# Exact integral of genz function in hypercube [0,1]^dim
real <- adaptIntegrate(copeak,lowerLimit = rep(0,dim), upperLimit = rep(1,dim))

# Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
percentageErrorBART <- meanValueBART - real[[1]]
percentageErrorMonteCarlo <- meanValueMonteCarlo - real[[1]]
percentageErrorGP <- meanValueGP - real[[1]]

# 1. Open jpeg file
jpeg("report/Figures//convergenceCopeak.jpg", width = 700, height = 583)
# 2. Create the plot
plot(x = log(c(1:400)), y = log(standardDeviationMonteCarlo$standardDeviation2),
     pch = 16, frame = FALSE, type = "l",
     xlab = "number of samples N", ylab = "log standard deviation", col = "blue",
     main = "Convergence of methods: log(sd) vs log(N) using Copeak")
lines(x = log(c(1:400)), log(predictionBART$standardDeviation), type = 'l', col = "red")
lines(x = log(c(1:400)), log(predictionGPBQ$standardDeviationGP), type = 'l', col = "green")
legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ"),
       col=c("blue", "red", "green"), lty=1:2, cex=0.8)
# 3. Close the file
dev.off()


#[5, 4, 3, 2, 1.5, 1.3, 1.2, 1.1]
# mean vs n
#log sd vs log n


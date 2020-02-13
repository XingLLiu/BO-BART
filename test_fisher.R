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
num_iterations <- as.double(args[2])
whichKernel <- "matern32"
sequential <- FALSE
# prior measure over the inputs
# uniform by default
measure <- "uniform"
print(c(dim, num_iterations))

source("src/genz/fisher_integrands.R")
# C <- replicate(dim, runif(9, 0.1, 0.9))
C <- rep(0.45958094818045914, dim)
R <- rep(0.13332051229486602, dim)
H <- rep(1.5020584867022158, dim)
F <- rep(4.624266954512351, dim)
P <- rep(1, dim)


fisher_function <- create_fisher_function(C, R, H, F, P, dim)
# prepare training dataset
if (measure == "uniform") {
  trainX <- replicate(dim, runif(100, 0, 1))
  trainY <- fisher_function(trainX)
} else if (measure == "gaussian") {
  trainX <- replicate(dim, rtnorm(100, mean=0.5, lower=0, upper=1))
  trainY <- fisher_function(trainX)
}
fisher_1d <- create_fisher_function(C[1], R[1], H[1], F[1], P[1], 1)
real <- estimate_real_integral(fisher_1d, 1, 1e7)^dim

for (num_cv in 1:5) {
  cat("NUM_CV", num_cv, "\n")
  # Bayesian Quadrature method
  # set number of new query points using sequential design
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
    csvName <- "results/fisher_function/Dim%sNoSequential%s_%s.csv" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
    figName <- "Figures/fisher_function/Dim%sNoSequential%s_%s.pdf" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
  } else {
    csvName <- "results/fisher_function/Dim%s%s_%s.csv" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
    figName <- "Figures/fisher_function/Dim%s%s_%s_.pdf" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_cv
    )
  }
  
  results_models <- list("BART"=predictionBART, "GP"=predictionGPBQ, "MC"=predictionMonteCarlo)
  save(results_models, file = "results/fisher_function/Dim%s%s_%s.RData" %--% c(
    dim,
    tools::toTitleCase(measure),
    num_cv
  ))
  
  write.csv(results, file = csvName, row.names=FALSE)
  
  if (dim == 1) {
    
    # Plotting
    x_plot <- replicate(dim, runif(500))
    x_plot <- x_plot[order(x_plot),]
    y_pred <- predict(predictionBART$model, x_plot)
    y_pred_mean <- colMeans(y_pred)
    y_pred_sd <- sqrt(colVars(y_pred))
    
    # obtain posterior samples
    integrals <- sampleIntegrals(predictionBART$model, dim, measure)
    ymin <- min(predictionBART$trainData[, (dim + 1)]); ymax <- max(predictionBART$trainData[, (dim + 1)])
    integrals <- (integrals + 0.5) * (ymax - ymin) + ymin
    
    K <- predictionGPBQ$K
    X <- predictionGPBQ$X
    Y <- predictionGPBQ$Y
    Y <- Y[order(X)]
    maternKernel <- maternKernelWrapper(lengthscale = lengthscale)
    
    k_xstar_x <- kernelMatrix(maternKernel, matrix(x_plot, ncol=1), X)
    k_xstar_xstar <- kernelMatrix(maternKernel, 
                                  matrix(x_plot, ncol=1), 
                                  matrix(x_plot, ncol=1))
    jitter = 1e-6
    K_inv <- solve(K + diag(jitter, nrow(K)))
    
    gp_post_mean <- k_xstar_x %*% K_inv %*% Y
    gp_post_cov <- k_xstar_xstar - k_xstar_x %*% K_inv %*% t(k_xstar_x)
    gp_post_sd <- sqrt(diag(gp_post_cov))
    
    #plot of integrals
    GPdensity <- dnorm(
      seq(0, 1, 0.01), 
      mean = predictionGPBQ$meanValueGP[1], 
      sd = sqrt(predictionGPBQ$varianceGP[1])
    )
    plot(x_plot, y_pred_mean)
    hist(integrals)
    KDE_BART <- density(integrals)
    
    pdf(figName, width = 12, height = 8)
    par(mfrow = c(1,2), pty = "s", cex=1.5)
    plot(
      seq(0, 1, 0.01), 
      GPdensity, 
      ty="l", 
      col = "dodgerblue", 
      xlim = c(0,1), 
      ylim = c(0, 60),
      xlab = "x",
      ylab = "Posterior density",
      cex.lab = 1.5,
      cex.axis = 1.5,
      lwd=3
    )
    points(KDE_BART, ty="l", col = "orangered", lwd=3)
    abline(v=real)
    legend("topright", legend=c("BART-Int", "GP-BQ", "Actual"),
           col=c("orangered", "dodgerblue", "black"), cex=0.6, lty = c(1,1,1), lwd=3)
    
    a <-density(integrals)$y 
    plot(x_plot, 
         gp_post_mean, 
         col = "dodgerblue", 
         cex=0.5, 
         ty="l", 
         ylim=c(-2, 4),
         xlab = "x",
         cex.lab = 1.5,
         cex.axis = 1.5,
         ylab = "y",
         cex.lab = 1.5,
         cex.axis = 1.5,
         lwd=3
    )
    points(trainX[order(trainX),], trainY[order(trainX)], ty="l", lwd=3)
    points(trainX[order(trainX),], trainY[order(trainX)], col = "black", bg='black', pch=21, lwd=3, cex=0.5)
    polygon(c(x_plot, rev(x_plot)), 
            c(
              gp_post_mean + 2*gp_post_sd, 
              rev(gp_post_mean - 2*gp_post_sd)
            ), 
            col = adjustcolor("dodgerblue", alpha.f = 0.10), 
            border = "dodgerblue", lty = c("dashed", "solid"))
    # points(trainX, trainY, col = "blue")
    polygon(c(x_plot, rev(x_plot)), 
            c(
              y_pred_mean + 2*y_pred_sd, 
              rev(y_pred_mean - 2*y_pred_sd)
            ), 
            col = adjustcolor("orangered", alpha.f = 0.10),  
            border = "orangered", lty = c("dashed", "solid"))
    points(x_plot, y_pred_mean, col = "orangered", cex=0.5, ty="l", lwd=3)
  }
}




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
num_data <- as.double(args[2])
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

# C <- replicate(1, runif(3, 0.1, 0.9))
# R <- replicate(1, rbeta(3, 5, 2))
# H <- replicate(1, runif(3, 0.5*exp(1), 1.5*exp(1)))
# F <- replicate(1, runif(3, 0, 5))
# P <- replicate(1, rbinom(3, 1, 0.5))
C <- c(0.1706093, 0.5319923, 0.7117816)
R <- c(0.6786221, 0.7207544, 0.5120249)
H <- c(3.029867, 3.065427, 3.114357)
F <- c(4.607562, 4.526243, 2.221768)
P <- c(1, 0, 1)

cut_point <- 0.5
fisher_function_full <- create_fisher_function(C, R, H, F, P, 3)

fisher_function <- function(x) {
  x_in <- cbind(x, matrix(cut_point, nrow = nrow(x), ncol = 3-dim))
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
if (dim != 3) {
  for (i in (dim+1):3) {
    fisher_1d <- create_fisher_function(C[i], R[i], H[i], F[i], P[i], 1)
    print(fisher_1d(matrix(cut_point)))
    real <- real * fisher_1d(matrix(cut_point))
  }
}
  

for (num_cv in 1:20) {
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
    csvName <- "results/fisher_function/Dim%sNoSequential%s_%s_%s.csv" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_data,
      num_cv
    )
    figName <- "Figures/fisher_function/Dim%sNoSequential%s_%s_%s.pdf" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_data,
      num_cv
    )
    figName_convergence <- "Figures/fisher_function/convergence_Dim%sNoSequential%s_%s_%s.pdf" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_data,
      num_cv
    )
  } else {
    csvName <- "results/fisher_function/Dim%s%s_%s_%s.csv" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_data,
      num_cv
    )
    figName <- "Figures/fisher_function/Dim%s%s_%s_%s.pdf" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_data,
      num_cv
    )
    figName_convergence <- "Figures/fisher_function/convergence_Dim%s%s_%s_%s.pdf" %--% c(
      dim,
      tools::toTitleCase(measure),
      num_data,
      num_cv
    )
  }
  
  results_models <- list("BART"=predictionBART, "GP"=predictionGPBQ, "MC"=predictionMonteCarlo)
  save(results_models, file = "results/fisher_function/Dim%s%s_%s_%s.RData" %--% c(
    dim,
    tools::toTitleCase(measure),
    num_data,
    num_cv
  ))
  
  write.csv(results, file = csvName, row.names=FALSE)
  
  if (dim == 1) {
    
    # Plotting
    x_plot <- matrix(cut_point, nrow = 300, ncol = 3)
    x_plot[,1:dim] <- replicate(dim, runif(300))
    y_plot <- fisher_function(x_plot)
    x_plot <- x_plot[,1:dim]
    y_plot <- y_plot[order(x_plot)]
    x_plot <- x_plot[order(x_plot)]
    y_pred <- predict(predictionBART$model, x_plot)
    y_pred_mean <- colMeans(y_pred)
    y_pred_sd <- sqrt(colVars(y_pred))
    plot(x_plot, y_plot, ty="l")
    
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
    K_inv <- chol2inv(chol(K + diag(jitter,nrow(K))))
    
    gp_post_mean <- k_xstar_x %*% K_inv %*% Y
    gp_post_cov <- k_xstar_xstar - k_xstar_x %*% K_inv %*% t(k_xstar_x)
    gp_post_sd <- sqrt(diag(gp_post_cov))
    
    #plot of integrals
    GPdensity <- dnorm(
      seq(0, 1, 0.01), 
      mean = predictionGPBQ$meanValueGP[1], 
      sd = sqrt(predictionGPBQ$varianceGP[1])
    )
    plot(x_plot, y_pred_mean, ty="l")
    plot(x_plot, gp_post_mean, ty="l")
    
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
           col=c("orangered", "dodgerblue", "black"), cex=1, lty = c(1,1,1), lwd=3)
    
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
    points(x_plot[order(x_plot)], y_plot[order(x_plot)], ty="l", lwd=3)
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
    dev.off()
  }
}
num_data_lists <- list()
num_data_lists[[1]] <- c(20,50,100)
num_data_lists[[2]] <- c(50,100,200)
for (dim in 1:2) {
  num_data_list <- num_data_lists[[dim]]
  results <- data.frame(
  "num_data" = num_data_list,
  "BARTMape" = num_data_list,
  "BARTSE" = num_data_list,
  "GPMape" = num_data_list,
  "GPSE" = num_data_list,
  "MIMape" = num_data_list,
  "MISE" = num_data_list
  )
  for (num_data in num_data_list) {
    n <- 1
    mape <- 0
    sd <- 0
    bart_estimates <-c()
    gp_estimates <-c()
    mi_estimates <-c()
    for (num_cv in 1:20) {
      csv_name <- paste(
      "results/fisher_function/",
      "Dim",
      dim,
      "NoSequentialUniform_",
      num_data,
      "_",
      num_cv,
      ".csv",
      sep=""
      )
      result <- read.csv(csv_name)
      real <- result$actual
      mape <- mape + 1/20 * abs((c(result$BARTMean, result$GPMean, result$MIMean)-real)/real)
      bart_estimates[num_cv] <- result$BARTMean
      gp_estimates[num_cv] <- result$GPMean
      mi_estimates[num_cv] <- result$MIMean
    }
    se[1] <- 1/sqrt(num_cv) * sd(bart_estimates)
    se[2] <- 1/sqrt(num_cv) * sd(gp_estimates)
    se[3] <- 1/sqrt(num_cv) * sd(mi_estimates)
    results[n, c(2,4,6)] <- mape
    results[n, c(3,5,7)] <- se
    n <- n + 1
  }
  result_csv_name <- paste("results/fisher_function/results_fisher_function_dim_", dim, ".csv", sep="")
  write.csv(results, result_csv_name)
}


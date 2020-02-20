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

# define string formatting
`%--%` <- function(x, y) 
# from stack exchange:
# https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}


genz_pipeline_nosequential <- function(args)
{
  epochs <- args$epochs
  num_cv_total <- args$num_runs

  dim <- args$dim
  num_iterations <- args$num_iterations
  whichGenz <- args$whichGenz
  whichKernel <- args$whichKernel
  # turn on/off sequential design
  # 1 denotes TRUE to sequential
  # 0 denotes FALSE to sequential
  cat("\nBegin testing:\n")
  if (as.double(args$sequential) == 1 | is.null(as.double(args$sequential))) {
    sequential <- TRUE
  } else {
    sequential <- FALSE
  }
  cat("Sequantial design set to", sequential, "\n")
  # prior measure over the inputs
  # uniform by default
  if (args$measure != "gaussian" | is.null(args$measure)) {
    measure <- "uniform"
  } else{
    measure <- as.character(args$measure)
  }
  cat("Prior measure:", measure, "\n")

  # extra parameter for step function
  # 1 by default
  if (whichGenz == 7 & is.null(args$jumps)) {
    jumps <- 1
    cat("Number of jumps for step function:", jumps, "\n")
  } else if (whichGenz == 7){
    jumps <- as.double(args$jumps)
    cat("Number of jumps for step function:", jumps, "\n")
  }

  # extra parameter for additive Gaussian function
  if (whichGenz == 9){ add_gauss_a <- NA}

  if (num_iterations == 1 & sequential) { stop("NEED MORE THAN 1 ITERATION") }

  print(c(dim, num_iterations, whichGenz))
  source("src/genz/genz.R") # genz function to test

  if (whichGenz < 1 | whichGenz > 9) { stop("undefined genz function. Change 3rd argument to 1-9") }
  if (whichGenz == 3 & dim == 1) { stop("incorrect dimension. Discrete Genz function only defined for dimension >= 2") } 

  if (whichGenz == 1) { genz <- cont; genzFunctionName <-  deparse(substitute(cont)) }
  if (whichGenz == 2) { genz <- copeak; genzFunctionName <-  deparse(substitute(copeak)) }
  if (whichGenz == 3) { genz <- disc; genzFunctionName <-  deparse(substitute(disc)) }
  if (whichGenz == 4) { genz <- gaussian; genzFunctionName <-  deparse(substitute(gaussian)) }
  if (whichGenz == 5) { genz <- oscil; genzFunctionName <-  deparse(substitute(oscil)) }
  if (whichGenz == 6) { genz <- prpeak; genzFunctionName <-  deparse(substitute(prpeak)) }
  if (whichGenz == 7) { genz <- function(xx){return(step(xx, jumps=jumps))}; genzFunctionName <-  deparse(substitute(step)) }
  if (whichGenz == 8) { genz <- mix; genzFunctionName <-  deparse(substitute(mix)) }
  if (whichGenz == 9) { genz <- function(xx){return(additive_gaussian(xx, a=add_gauss_a))}; genzFunctionName <-  deparse(substitute(additive_gaussian)) }

  print("Testing with: %s" %--% genzFunctionName)

  # Initialize dataframe
  BARTMean <- BARTsd <- rep(NA, length(epochs))
  MIMean <- MIsd <- rep(NA, length(epochs))
  GPMean <- GPsd <- rep(NA, length(epochs))
  runtimeBART <- runtimeMI <- runtimeGP <- rep(NA, length(epochs))

  for (num_cv in 1:num_cv_total) {

    for (iter in 1:length(epochs)) {
      print("Epochs: [%i/%i]; CV: [%i/%i]" %--% c(iter, length(epochs), num_cv, num_cv_total))

      # prepare training dataset
      if (iter == 1){
        if (measure == "uniform") {
            trainX <- replicate(dim, runif(epochs[iter]))
            trainY <- genz(trainX)
        } else if (measure == "gaussian") {
            trainX <- replicate(dim, rtnorm(epochs[iter], mean=0.5, lower=0, upper=1))
            genz <- gaussian_weighted
            trainY <- genz(trainX)
        }
      } else {
        if (measure == "uniform") {
            trainX <- rbind(trainX, runif(epochs[iter] - epochs[iter - 1]))
            # trainX <- rbind(trainX, replicate(dim, runif(epochs[iter] - epochs[iter - 1])))
            trainY <- genz(trainX)
        } else if (measure == "gaussian") {
            trainX <- rbind(trainX, runif(epochs[iter] - epochs[iter - 1]))
            genz <- gaussian_weighted
            trainY <- genz(trainX)
        }
      }

      # Numerical integration methods
      # set number of new query points using sequential design
      source("src/BARTBQ.R")
      t0 <- proc.time()
      predictionBART <- mainBARTBQ(dim, num_iterations, FUN = genz, trainX, trainY, sequential, measure)
      t1 <- proc.time()
      runtimeBART[iter] <- (t1 - t0)[[1]]
      
      # Bayesian Quadrature with Monte Carlo integration method
      print("Begin Monte Carlo Integration")
      source("src/monteCarloIntegration.R")
      
      t0 <- proc.time()
      predictionMonteCarlo <- monteCarloIntegrationUniform(FUN = genz, trainX, trainY, numSamples=num_iterations, dim, measure)
      t1 <- proc.time()
      runtimeMI[iter] <- (t1 - t0)[[1]]
      
      # Bayesian Quadrature with Gaussian Process
      print("Begin Gaussian Process Integration")
      if (num_cv == 1){
        library(reticulate)
        source("src/optimise_gp.R")
        lengthscale <- optimise_gp_r(trainX, trainY, kernel = whichKernel, epochs=500)
        print("...Finished training for the lengthscale")
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
          FUN = genz, 
          lengthscale,
          sequential, 
          measure
      )  
      t1 <- proc.time()
      runtimeGP[iter] <- (t1 - t0)[[1]]

      # Store results
      BARTMean[iter] <- predictionBART$meanValueBART
      BARTsd[iter] <- predictionBART$standardDeviationBART
      MIMean[iter] <- predictionMonteCarlo$meanValueMonteCarlo
      MIsd[iter] <- predictionMonteCarlo$standardDeviationMonteCarlo
      GPMean[iter] <- predictionGPBQ$meanValueGP
      GPsd[iter] <- sqrt(predictionGPBQ$varianceGP)
    }

      # Read in analytical integrals
      source("src/genz/analyticalIntegrals.R")
      dimensionsList <- c(1,2,3,5,10,20)
      whichDimension <- which(dim == dimensionsList)
      if (whichGenz <= 6){
          analyticalIntegrals <- read.csv("results/genz/integrals.csv", header = FALSE)
          real <- analyticalIntegrals[whichGenz, whichDimension]
      } else if (whichGenz == 7) {
          real <- stepIntegral(dim, jumps)
      } else if (whichGenz == 8) {
          if (dim ==1){ real <- 0.008327796}
          if (dim ==2){ real <- 0.008327796 * 2}
          if (dim ==3){ real <- 0.008327796 * 3}
      } else if (whichGenz == 9) {
          real <- additiveGaussianIntegral(dim, a = add_gauss_a)
      }
      
      # # Bayesian Quadrature methods: with BART, Monte Carlo Integration and Gaussian Process respectively
      # print("Final Results:")
      # print(c("Actual integral:", real))
      # print(c("BART integral:", predictionBART$meanValueBART[num_iterations]))
      # print(c("MI integral:", predictionMonteCarlo$meanValueMonteCarlo[num_iterations]))
      # print(c("GP integral:", predictionGPBQ$meanValueGP[num_iterations]))
      
      print("Writing full results to results/genz%s" %--% c(whichGenz))
      results <- data.frame(
              "epochs" = epochs,
              "BARTMean" = BARTMean, "BARTsd" = BARTsd,
              "MIMean" = MIMean, "MIsd" = MIsd,
              "GPMean" = GPMean, "GPsd" = GPsd,
              "actual" = rep(real, num_iterations),
              "runtimeBART" = runtimeBART,
              "runtimeMI" = runtimeMI,
              "runtimeGP" = runtimeGP
      )
      # results_models <- list("BART"=predictionBART, "GP"=predictionGPBQ, "MC"=predictionMonteCarlo)
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
          # save(results_models, file = "results/genz/%s/%sDim%sNoSequential%s_%s.RData" %--% c(
          #         whichGenz,
          #         genzFunctionName,
          #         dim,
          #         tools::toTitleCase(measure),
          #         num_cv
          # ))
      } else {
          csvName <- "results/genz/%s/%sDim%s%s_%s.csv" %--% c(
                  whichGenz, 
                  genzFunctionName,
                  dim,
                  tools::toTitleCase(measure),
                  num_cv
          )
          figName <- "Figures/%s/%sDim%s%s_%s.pdf" %--% c(
                  whichGenz,
                  genzFunctionName,
                  dim,
                  tools::toTitleCase(measure),
                  num_cv
          )
      }

    write.csv(results, file = csvName, row.names=FALSE)
  }

  # plot results
  print("Start plotting results...")
  source("src/genz/drawGraphs_CV2.R")
  plot_args <- list(dims_list = c(dim), genz_list = c(9), sequential_list = c("NoSequential"),
                    measure = measure, num_cv = num_cv_total)
  plot_results(plot_args)

  print(paste("...done. Please find the plots in", "Figures/%s/" %--% c(whichGenz)))
}

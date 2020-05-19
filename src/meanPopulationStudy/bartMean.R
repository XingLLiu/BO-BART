library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
# library(docstring)

# define string formatting
`%--%` <- function(x, y) 
# from stack exchange:
# https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

computeBART <- function(trainX, trainY, candidateX, candidateY, num_iterations) 
#' BART-BQ for estimating average income
#' @description Compute mean for BART-BQ with
#' implementation of query sequential design 
#' of adding more training data to the original dataset.
#' @param trainX data frame. Covariates of training set.
#' @param trainY data frame. Regressors of training set.
#' number of rows must agree with number of rows of trianX.
#' @param candidateX data frame. Covariates of candidate set.
#' @param candidateY data frame. Regressor of candidate set.
#' @param num_iterations numeric. Number of new training data to be added.
#' @details For every iteration, we compute the test error.
#' @return A list of mean integral value, standard deviation of integral 
#' value (segmented and not segmented) and new training set. 
{
  dim <- ncol(trainX)
  print( c("Adding number of new training data:", num_iterations) )
  # outputs
  meanValue <- rep(0, num_iterations)
  standardDeviation <- rep(0, num_iterations)
  eduMeanValue <- matrix(0, 2*num_iterations, nrow=2, ncol=num_iterations)
  eduStandardDeviation <- matrix(0, 2*num_iterations, nrow=2, ncol=num_iterations)
  trainData <- cbind(trainX, trainY)
  fullData <- rbind(trainX, candidateX)
  
  colnames(trainData)[dim+1] <- "response"

  # index of segmentation
  index <- matrix(c(1,16,17,24), nrow = 2, ncol = 2)

  # generate extra training data using the scheme (see pdf)
  for (i in 1:num_iterations) {
    
    # set seed to enable reproduction of the results
    print(c("BART: Epoch=", i))
    # first build BART model
    sink("/dev/null")
    # model <- bart(trainData[, 1:dim], trainData[, dim+1], keeptrees=TRUE, keepevery=3L, 
    #               nskip=500, ndpost=2000, ntree=50, k=2, usequant=FALSE)
    model <- bart(trainData[,1:dim], trainData[, dim+1], keeptrees=TRUE, keepevery=3L, 
                  nskip=200, ndpost=5000, ntree=50, k=3, usequant=FALSE)              
    sink()

    # predict the values
    fValues <- predict(model, candidateX)

    # select the best candidate, find its response
    var <- colVars(fValues)
    index <- sample(which(var==max(var)), 1)
    response <- candidateY[index]

    # # data segmentation by education level
    # mat <- matrix(c(1,16,17,24), nrow = 2, ncol = 2)

    # for (cat in 1:2){

    #   eduCandidateX <- candidateX[(candidateX$Education >= mat[1, cat] & candidateX$Education <= mat[2, cat]), ]
      
    #   # make prediction
    #   fValues <- predict(model, eduCandidateX)

    #   # posterior mean and variance
    #   integrals <- rowMeans(fValues)
    #   eduMeanValue[cat, i] <- mean(integrals)
    #   eduStandardDeviation[cat, i] <- sd(integrals)

    # }
    
    # add new data to train set
    trainData <- rbind(trainData, cbind(candidateX[index, ], response))
    
    # Integral with respect to \Pi_n
    pred <- predict(model, fullData)
    meanValue[i] <- mean(pred)
    # standardDeviation[i] <- sd(trainData[, dim+1]) 
    standardDeviation[i] <- sum((rowMeans(pred) - meanValue[i])^2) / (nrow(pred) - 1)

    # remove newly added value from candidate set
    candidateX <- candidateX[-index,]
    candidateY <- candidateY[-index]

  }

  return(list("meanValueBART"=meanValue, "standardDeviationBART"=standardDeviation, 
              "eduMeanValueBART"=eduMeanValue, "eduStandardDeviationBART"=eduStandardDeviation, "trainData"=trainData))
}


stratified <- function(df, group, size, select=NULL, replace=FALSE, bothSets=FALSE) 
#' Stratified sampling
#' @description Method of stratification from CRAN package fifer.
{
  if (is.null(select)) {
    df <- df
  } else {
    if (is.null(names(select))) stop("'select' must be a named list")
    if (!all(names(select) %in% names(df)))
      stop("Please verify your 'select' argument")
    temp <- sapply(names(select),
                   function(x) df[[x]] %in% select[[x]])
    df <- df[rowSums(temp) == length(select), ]
  }
  df.interaction <- interaction(df[group], drop = TRUE)
  df.table <- table(df.interaction)
  df.split <- split(df, df.interaction)
  if (length(size) > 1) {
    if (length(size) != length(df.split))
      stop("Number of groups is ", length(df.split),
           " but number of sizes supplied is ", length(size))
    if (is.null(names(size))) {
      n <- setNames(size, names(df.split))
      message(sQuote("size"), " vector entered as:\n\nsize = structure(c(",
              paste(n, collapse = ", "), "),\n.Names = c(",
              paste(shQuote(names(n)), collapse = ", "), ")) \n\n")
    } else {
      ifelse(all(names(size) %in% names(df.split)),
             n <- size[names(df.split)],
             stop("Named vector supplied with names ",
                  paste(names(size), collapse = ", "),
                  "\n but the names for the group levels are ",
                  paste(names(df.split), collapse = ", ")))
    }
  } else if (size < 1) {
    n <- round(df.table * size, digits = 0)
  } else if (size >= 1) {
    if (all(df.table >= size) || isTRUE(replace)) {
      n <- setNames(rep(size, length.out = length(df.split)),
                    names(df.split))
    } else {
      message(
        "Some groups\n---",
        paste(names(df.table[df.table < size]), collapse = ", "),
        "---\ncontain fewer observations",
        " than desired number of samples.\n",
        "All observations have been returned from those groups.")
      n <- c(sapply(df.table[df.table >= size], function(x) x = size),
             df.table[df.table < size])
    }
  }
  temp <- lapply(
    names(df.split),
    function(x) df.split[[x]][sample(df.table[x],
                                     n[x], replace = replace), ])
  set1 <- do.call("rbind", temp)
  
  if (isTRUE(bothSets)) {
    set2 <- df[!rownames(df) %in% rownames(set1), ]
    list(SET1 = set1, SET2 = set2)
  } else {
    set1
  }
}


computeMI <- function(trainX, trainY, candidateX, candidateY, num_iterations, seed=NA)
#' Simple random sampling (Monte Carlo)
#' @description Compute sample mean of response by simple random sampling 
#' @param trainX data frame. Covariates of training set.
#' @param trainY data frame. Regressors of training set.
#' number of rows must agree with number of rows of trianX.
#' @param candidateX data frame. Covariates of candidate set.
#' @param candidateY data frame. Regressor of candidate set.
#' @param num_iterations numeric. Number of new training data to be added.
#' @param seed numeric. Seed for shuffling the data; no shuffling if not given.
#' @return A list of mean integral value, standard deviation of integral 
#' value (segmented and not segmented) and new training set. 
{

  if (!is.na(seed)){
    set.seed(seed)
    candidateY <- sample(candidateY)
  }
  MImean <- rep(NA, num_iterations)
  combinedY <- c(trainY, candidateY[1:num_iterations])
  combinedN <- length(combinedY)
  
  MImean <- cumsum(combinedY) / (1:combinedN)
  cumvar <- (cumsum(combinedY^2) - MImean^2 * (1:combinedN)) / (0:(combinedN - 1))
  MIstandardDeviation <- sqrt(cumvar) / sqrt(1:combinedN)
  
  MImean <- MImean[-(1:length(trainY))]
  MIstandardDeviation <- MIstandardDeviation[-c(1:length(trainY))]

  return(list("meanValueMI"=MImean, "standardDeviationMI"=MIstandardDeviation))

}


computeBRS <- function(trainX, trainY, candidateX, candidateY, group, num_iterations)
#' Block random sampling by group
#' @description Compute sample mean of response by block random sampling, 
#' blocks are buildt based on group.
#' @param trainX data frame. Covariates of training set.
#' @param trainY data frame. Regressors of training set.
#' number of rows must agree with number of rows of trianX.
#' @param candidateX data frame. Covariates of candidate set.
#' @param candidateY data frame. Regressor of candidate set.
#' @param group string. Name of the attribute to be stratified on.
#' @param num_iterations numeric. Number of new training data to be added.
#' @return A list of mean integral value, standard deviation of integral 
#' value (segmented and not segmented) and new training set. 
{
    dim <- ncol(trainX)
    data <- cbind(candidateX, candidateY)
    # sample query points
    samples <- stratified(data, group, num_iterations/length(candidateY))
    # randomise the order of the samples
    samples <- samples[sample(nrow(samples)), ]

    BRmean <- c()
    BRstandardDeviation <- c()

    # Monte Carlo in each block 
    # for (i in 1:num_iterations) {
      
    #     BRmean[i] <- mean(c(trainY, samples[1:i, dim+1]))

    #     BRstandardDeviation[i] <- sd(c(trainY, samples[1:i, dim+1])) /sqrt(length(c(trainY, samples[1:i, dim+1])))

    # }

    # Monte Carlo in each block 
    BRmean <- rep(NA, num_iterations)
    combinedY <- c(trainY, samples[1:num_iterations, dim + 1])
    combinedN <- length(combinedY)
    
    BRmean <- cumsum(combinedY) / (1:combinedN)
    cumvar <- (cumsum(combinedY^2) - BRmean^2 * (1:combinedN)) / (0:(combinedN - 1))
    BRstandardDeviation <- sqrt(cumvar) / sqrt(1:combinedN)
    
    BRmean <- BRmean[-(1:length(trainY))]
    BRstandardDeviation <- BRstandardDeviation[-c(1:length(trainY))]
    
    return(list("meanValueBRS"=BRmean, "standardDeviationBRS"=BRstandardDeviation))
}
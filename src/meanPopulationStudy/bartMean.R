library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)

# define string formatting
`%--%` <- function(x, y) 
# from stack exchange:
# https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

computePopulationMean <- function(dim, trainX, trainY, condidateX, candidateY, num_iterations) 
# compute mean for BART-BQ with
# implementation of query sequential design to add
# more training data to the original dataset
# For every iteration, we compute the test error
# input:
#   dim: dimension
#   trainX: covariates of training data
#   trainY: response of training data
#   numNewTraining: number of new training points to be added
#
# output:
#   list of mean integral value, standard deviation of integral value and new traiing set
{
  print( c("Adding number of new training data:", num_iterations) )
  # outputs
  meanValue <- rep(0, num_iterations)
  standardDeviation <- rep(0, num_iterations)
  trainData <- cbind(trainX, trainY)
  eduMeanValue <- matrix(0, 2*num_iterations, nrow=2, ncol=num_iterations)
  eduStandardDeviation <- matrix(0, 2*num_iterations, nrow=2, ncol=num_iterations)
  
  colnames(trainData)[dim+1] <- "response"

  # index of segmentation
  index <- matrix(c(1,16,17,24), nrow = 2, ncol = 2)

  # generate extra training data using the scheme (see pdf)
  for (i in 1:num_iterations) {

    set.seed(i)

    print(c("BART: Epoch=", i))
    # first build BART model
    sink("/dev/null")
    model <- bart(trainData[, 1:dim], trainData[, dim+1], keeptrees=TRUE, keepevery=5L, nskip=100, ndpost=200, ntree=50, k=10, usequant=FALSE) #50 10; 60 8; 100 10 
    sink()

    # predict the values
    fValues <- predict(model, candidateX)

    # posterior mean and variance
    integrals <- rowMeans(fValues)
    meanValue[i] <- mean(integrals)
    standardDeviation[i] <- sd(integrals) 

    # select the best candidate, find its response
    var <- colVars(fValues)
    index <- sample(which(var==max(var)), 1)
    response <- candidateY[index]

    # data segmentation by education level
    mat <- matrix(c(1,16,17,24), nrow = 2, ncol = 2)

    for (cat in 1:2){

      eduCandidateX <- candidateX[(candidateX$Education >= mat[1, cat] & candidateX$Education <= mat[2, cat]), ]
      
      # make prediction
      fValues <- predict(model, eduCandidateX)

      # posterior mean and variance
      integrals <- rowMeans(fValues)
      eduMeanValue[cat, i] <- mean(integrals)
      eduStandardDeviation[cat, i] <- sd(integrals)

    }
    
    # add new data to train set
    trainData <- rbind(trainData, cbind(candidateX[index, ], response))

    # remove newly added value from candidate set
    candidateX <- candidateX[-index,]
    candidateY <- candidateY[-index]

  }

  return (list("meanValueBART"=meanValue, "standardDeviationBART"=standardDeviation, 
               "eduMeanValueBART"=eduMeanValue, "eduStandardDeviationBART"=eduStandardDeviation, "trainData"=trainData))
}

computeBART <- function(trainX, trainY, candidateX, candidateY, num_iterations) 
# main method
# input:
#   dim
#   num_iterations:
#   trainX: covariates of training set
#   trainY: regressors of training set
#   testX: covariates of test set
#   testY: regressor of test set
#
# returns prediction as a list
{
  # prepare training data and parameters
    numNewTraining <- num_iterations
    dim <- ncol(trainX)
    # compute population mean response
    BARTResults <- computePopulationMean(dim, trainX, trainY, candidateX, candidateY, num_iterations=numNewTraining) 

    return (BARTResults)
}

stratified <- function(df, group, size, select=NULL, replace=FALSE, bothSets=FALSE) {

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

computeMI <- function(trainX, trainY, candidateX, candidateY, num_iterations)
{

  MImean <- c()
  MIstandardDeviation <- c()

  for (i in 1:num_iterations) {
    
    MImean[i] <- mean(c(trainY, candidateY[1:i]))

    n = length(c(trainY, candidateY[1:i]))
    MIstandardDeviation[i] <- sqrt( var(c(trainY, candidateY[1:i]))/(n-1) )

  }

  return(list("meanValueMI"=MImean, "standardDeviationMI"=MIstandardDeviation))

}


computeBRS <- function(trainX, trainY, candidateX, candidateY, group, num_iterations){

    dim <- ncol(trainX)
    data <- cbind(candidateX, candidateY)
    samples <- stratified(data, group, num_iterations/length(candidateY))
    samples <- samples[sample(nrow(samples)), ]

    BRmean <- c()
    BRstandardDeviation <- c()

    # Monte Carlo in each block 
    for (i in 1:num_iterations) {
      
        BRmean[i] <- mean(c(trainY, samples[1:i, dim+1]))

        n = length(c(trainY, samples[1:i, dim+1]))
        BRstandardDeviation[i] <- sqrt( var(c(trainY, samples[1:i, dim+1])) /(n-1) )

    }
    
    return(list("meanValueBRS"=BRmean, "standardDeviationBRS"=BRstandardDeviation))

}
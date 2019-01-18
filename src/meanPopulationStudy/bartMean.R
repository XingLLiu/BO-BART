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

# Tree code

terminalProbability <- function(currentNode) 
  # probabiltity ending up in terminal node
{
  prob <- currentNode$probability
  
  while ( !isRoot(currentNode$parent) ) {
    currentNode <- currentNode$parent
    prob <- prob*currentNode$probability
  }
  
  return (prob)
}

fillProbabilityForNode <- function(oneTree, cutPoints, cut) 
  # Drop data set into the tree and assign them to different nodes 
{
  if ( !is.null(oneTree$leftChild) ) {
    
    decisionRule <- cutPoints[[oneTree$splitVar]][oneTree$splitIndex]
    
    #oneTree$leftChild$probability <- (decisionRule - cut[1, oneTree$splitVar]) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
    oneTree$leftChild$probability <- pnorm(decisionRule) - pnorm(cut[1, oneTree$splitVar])
  
    #oneTree$rightChild$probability <- (cut[2, oneTree$splitVar] - decisionRule) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
    oneTree$rightChild$probability <- pnorm(cut[2, oneTree$splitVar]) - pnorm(decisionRule)

    range <- cut[, oneTree$splitVar]
    
    cut[, oneTree$splitVar] = c(range[1], decisionRule)
    
    fillProbabilityForNode(oneTree$leftChild, cutPoints, cut)
    
    cut[, oneTree$splitVar] = c(decisionRule, range[2])
    
    fillProbabilityForNode(oneTree$rightChild, cutPoints, cut)
    
  } else if( is.null(oneTree$probability) ) {
    oneTree$probability <- 1
  }
  return (oneTree)
}

terminalProbabilityStore <- function(Tree)
  # store probability of getting to terminal node 
{
  terminalNodes = Traverse(Tree, filterFun = isLeaf)
  
  for (i in 1:length(terminalNodes)) {
    probability2 <- terminalProbability(terminalNodes[[i]])
    terminalNodes[[i]]$terminal_probability <- probability2
  }
  
  return (Tree)
}

getTree <- function(sampler, chainNum, sampleNum, treeNum)
  # create tree
{
  cutPoints <- dbarts:::createCutPoints(sampler)
  
  if (sampler$control@keepTrees) {
    treeString <- sampler$state[[chainNum]]@savedTrees[treeNum, sampleNum]
    treeFits <- sampler$state[[chainNum]]@savedTreeFits[, treeNum, sampleNum]
  }
  else {
    treeString <- sampler$state[[chainNum]]@trees[treeNum]
    treeFits <- sampler$state[[chainNum]]@treeFits[,treeNum]
  }                           
  
  tree <- dbarts:::buildTree(strsplit(gsub("\\.", "\\. ", treeString),
                                      " ", fixed = TRUE)[[1]])
  tree$remainder <- NULL
  
  tree$indices <- seq_len(length(sampler$data@y))
  tree <- dbarts:::fillObservationsForNode(tree, sampler, cutPoints)
  
  tree <- dbarts:::fillPlotInfoForNode(tree, sampler, treeFits)
  maxDepth <- dbarts:::getMaxDepth(tree)
  
  tree <- dbarts:::fillPlotCoordinatesForNode(tree, maxDepth, 1L, 1L)
  numEndNodes <- tree$index - 1L
  tree$index <- NULL
  
  return (tree)
}

singleTreeSum <- function(treeNum, model, drawNum, dim, trainX) 
  # Sum over a single tree's terminal nodes
{
  cutPoints<-dbarts:::createCutPoints(model$fit)
  
  trainX_mins <- apply(trainX,2,min)
  trainX_maxes <- apply(trainX,2,max)
  cut <- rbind(trainX_mins, trainX_maxes)
  
  treeList <- getTree(model$fit, 1, drawNum, treeNum)
  
  selectedTree <- FromListSimple(treeList) 
  
  #Modify tree by the functions written above
  selectedTree <- fillProbabilityForNode(selectedTree, cutPoints, cut)
  selectedTree <- terminalProbabilityStore(selectedTree)
  
  
  terminalNodeList <- Traverse(selectedTree, filterFun = isLeaf)
  
  #Calculate approximation of integreal in the single tree 
  integral <- 0
  for (node in terminalNodeList) {
    #We use the mean of prediction value Y's in the terminal node as u
    #rescale this step
    integral <- integral + node$terminal_probability * node$mu
  }
  return (integral)
}

posteriorSum <- function(drawNum, model, dim, trainX)
  # Sum over all the trees in one posterior draws
  # input:
  #   drawNum: which draw of p trees
  #   model:  set of tree
{
  nTree <- ncol(model$fit$state[[1]]@treeFits)
  treeNum <- seq(1, nTree, length.out=nTree)
  
  #Extra variables
  var <- list(model, drawNum, dim, trainX)
  
  #Calculate integration over all trees in the draw by mapply
  integral <- sum( unlist( mapply(singleTreeSum, treeNum, MoreArgs=var) ) )
  
  return (integral)
}


sampleIntegrals <- function(model, dim, trainX) 
  # sum over all posterior draws 
  # input: 
  #     model: BART model
  #
  # output:
  #     integrals: mean integral values for each tree as a vector
{
  nDraw <- dim(model$fit$state[[1]]@savedTreeFits)[3]
  drawNum <- seq(1, nDraw, length.out=nDraw)
  
  #Extra Variables
  var <- list(model, dim, trainX)
  integrals <- mapply(posteriorSum, drawNum, MoreArgs=var)
  return (integrals)
}

computeBART <- function(dim, trainX, trainY, condidateX, candidateY, numNewTraining) 
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
  
  print( c("Adding number of new training data:", numNewTraining ) )
  # outputs
  meanValue <- rep(0, numNewTraining)
  standardDeviation <- rep(0, numNewTraining)
  trainData <- cbind(trainX, trainY)
  colnames(trainData)[dim+1] <- "INCOME"
  
  set.seed(123)

  # generate extra training data using the scheme (see pdf)
  for (i in 1:numNewTraining) {

    print(c("BART: Epoch=", i))
    # find the min and max range of y
    # ymin <- min(trainData[, dim+1]); ymax <- max(trainData[, dim+1])
    # first build BART and scale mean and standard deviation
    sink("/dev/null")
    model <- bart(trainData[,1:dim], trainData[,dim+1], keeptrees=TRUE, keepevery=5L, nskip=100, ndpost=50, ntree = 10, k = 5, usequant = TRUE)
    sink()
    # # obtain posterior samples
    # integrals <- sampleIntegrals(model, dim, trainData[, 1:dim])
    # integrals <- (integrals + 0.5) * (ymax - ymin) + ymin
    
    # meanValue[i] <- mean(integrals)
    # standardDeviation[i] <- sqrt( sum((integrals - meanValue[i])^2) / (length(integrals) - 1) )

    # predict the values
    fValues <- predict(model, candidateX)
    
    probability = 1 #uniform probability
    
    var <- colVars(fValues)
    index <- sample(which(var==min(var)), 1)
    INCOME <- candidateY[index]
    
    # remove newly added value from candidate set
    trainData <- rbind(trainData, cbind(candidateX[index, ], INCOME))

    meanValue[i] <- mean(trainData[, dim+1])
    standardDeviation[i] <- sqrt( var(trainData[, dim+1]) )

    candidateX <- candidateX[-index,]
    candidateY <- candidateY[-index]
    
  }

  return (list("meanValueBART"=meanValue, "standardDeviationBART"=standardDeviation, 
               "trainData" = trainData))
}

computePopulationMean <- function(trainX, trainY, candidateX, candidateY, num_iterations) 
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
    # compute population mean income
    BARTResults <- computeBART(dim, trainX, trainY, candidateX, candidateY, numNewTraining) 

    return (BARTResults)
}

BRcomputeMean <- function(trainX, trainY, candidateX, candidateY, num_iterations = num_new_surveys){

    # sample the population by sex
    maleCandidateY <- candidateY[candidateX$Sex == 1]
    femaleCandidateY <- candidateY[candidateX$Sex == 2]

    maleRatio <- sum(trainX$Sex == 1) / nrow(trainX)

    numMaleCandidate <- floor(num_iterations * maleRatio)
    numFemaleCandidate <- num_iterations - numMaleCandidate

    BRmean <- c()
    BRstandardDeviation <- c()

    # Monte Carlo in each block 
    for (i in 1:numMaleCandidate) {
      
        BRmean[i] <- mean(c(trainY, maleCandidateY[1:i]))
        BRstandardDeviation[i] <- sqrt( var(c(trainY, maleCandidateY[1:i])) )
    }

    for (i in 1:numFemaleCandidate) {
        BRmean[i+numMaleCandidate] <- mean(c(trainY, maleCandidateY[1:numMaleCandidate], femaleCandidateY[1:i]))
        BRstandardDeviation[i+numMaleCandidate] <- sqrt( var(c(trainY, maleCandidateY[1:numMaleCandidate], femaleCandidateY[1:i])) )
    }
    
    return(list("BRmean" = BRmean, "BRstandardDeviation" = BRstandardDeviation))

}
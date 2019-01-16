###
# Functions for BART
###
library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)

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

    oneTree$leftChild$probability <- (decisionRule - cut[1, oneTree$splitVar]) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
    oneTree$rightChild$probability <- (cut[2, oneTree$splitVar] - decisionRule) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
    
    cut[, oneTree$splitVar] = c(0, decisionRule)
    
    fillProbabilityForNode(oneTree$leftChild, cutPoints, cut)
    
    cut[, oneTree$splitVar] = c(decisionRule, 1)
    
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

singleTreeSum <- function(treeNum, model, drawNum, dim) 
# Sum over a single tree's terminal nodes
{
  cutPoints<-dbarts:::createCutPoints(model$fit)
  cut <- array(c(0, 1), c(2,dim))
  
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

posteriorSum <- function(drawNum, model, dim)
# Sum over all the trees in one posterior draws
# input:
#   drawNum: which draw of p trees
#   model:  set of tree
{
  nTree <- ncol(model$fit$state[[1]]@treeFits)
  treeNum <- seq(1, nTree, length.out=nTree)
  
  #Extra variables
  var <- list(model, drawNum, dim)
  
  #Calculate integration over all trees in the draw by mapply
  integral <- sum( unlist( mapply(singleTreeSum, treeNum, MoreArgs=var) ) )

  return (integral)
}


sampleIntegrals <- function(model, dim) 
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
  var <- list(model, dim)
  integrals <- mapply(posteriorSum, drawNum, MoreArgs=var)
  return (integrals)
}

BARTBQSequential <- function(dim, trainX, trainY, numNewTraining, FUN) 
# compute integral for BART-BQ with
# implementation of query sequential design to add
# more training data to the original dataset
# input:
#   dim: dimension
#   trainX: covariates of training data
#   trainY: response of training data
#   numNewTraining: number of new training points to be added
#   FUN: function that we are integrating over
#
# output:
#   list of mean integral value, standard deviation of integral value and new traiing set
{
  
  print( c("Adding number of new training data:", numNewTraining ) )
  # outputs
  meanValue <- rep(0, numNewTraining)
  standardDeviation <- rep(0, numNewTraining)
  trainData <- data.frame(trainX, trainY)
  
  # generate extra training data using the scheme (see pdf)
  for (i in 1:numNewTraining) {
    
    print(c("BART: Epoch=", i))
    # find the min and max range of y
    ymin <- min(trainY); ymax <- max(trainY)
    # first build BART and scale mean and standard deviation
    sink("/dev/null")
    model <- bart(trainData[,1:dim], trainData[,dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=1000, ntree=50, k = 5)
    sink()
    # obtain posterior samples
    integrals <- sampleIntegrals(model, dim)
    integrals <- (integrals + 0.5) * (ymax - ymin) + ymin
    meanValue[i] <- mean(integrals)
    standardDeviation[i] <- sqrt(sum((integrals - meanValue[i])^2) / (length(integrals) - 1))

    # sequential design section, where we build the new training data
    candidateSet <- randomLHS(1000, dim)
    
    # predict the values
    fValues <- predict(model, candidateSet)
    
    probability = 1 #uniform probability
    #expectedValue <- colMeans(fValues*probability)
    
    var <- colVars(fValues)
    index <- sample(which(var==max(var)), 1)
    value <- FUN(t(candidateSet[index,]))
    trainData <- rbind(trainData, c(candidateSet[index,], value))
}

  return (list("meanValueBART"=meanValue, "standardDeviationBART"=standardDeviation, 
               "trainData" = trainData))
}

mainBARTBQ <- function(dim, num_iterations, FUN, trainX, trainY) 
# main method
# input:
#   dim
#   num_iterations:
#   FUN:
#   trainX: covariates of training set
#   trainY: regressors of training set
#
# returns prediction as a list
{
  # prepare training data and parameters
  genz <- FUN #select genz function
  numNewTraining <- num_iterations
  prediction <- BARTBQSequential(dim, trainX, trainY, numNewTraining, FUN = genz) 

  return (prediction)
}


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
  }
  else {
    treeString <- sampler$state[[chainNum]]@trees[treeNum]
  }                           
  if (sampler$control@keepTrees) {
    treeFits <- sampler$state[[chainNum]]@savedTreeFits[, treeNum, sampleNum]
  }
  else {
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

singleTreeSum <- function(treeNum, model, drawNum) 
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
    integral <- integral + node$terminal_probability*(node$mu) 
  }
  return (integral)
}

posteriorSum <- function(drawNum, model)
# Sum over all the trees in one posterior draws
{
  integral <- 0
  nTree <- ncol(model$fit$state[[1]]@treeFits)
  treeNum <- seq(1, nTree, length.out=nTree)
  
  #Extra variables
  var <- list(model, drawNum)
  
  #Calculate integration over all trees in the draw by mapply
  integral <- sum( unlist( mapply(singleTreeSum, treeNum, MoreArgs=var, SIMPLIFY = TRUE) ) )
  
  return (integral)
}


sampleIntegrals <- function(model) 
# sum over all posterior draws 
{
  nDraw <- dim(model$fit$state[[1]]@savedTreeFits)[3]
  drawNum <- seq(1, nDraw, length.out=nDraw)
  
  #Extra Variables
  var <- list(model)
  
  integrals <- mapply(posteriorSum, drawNum, MoreArgs=var, SIMPLIFY = TRUE)
  
  return (integrals)
}

BARTBQSequential <- function(dim, trainX, trainY, numNewTraining=1) 
# compute integral for BART-BQ with
# implementation of query sequential design to add
# more training data to the original dataset
{
  
  print( c("Adding number of new training data:", numNewTraining ) )
  # outputs
  meanValue <- rep(0, numNewTraining)
  standardDeviation <- rep(0, numNewTraining)
  trainData <- data.frame(trainX, trainY)
  
  # generate extra training data using the scheme (see pdf)
  for (i in 1:numNewTraining) {
  
  print(c("Epoch=", i))
  
  # first build BART and scale mean and standard deviation
  sink("/dev/null")
  model <- bart(trainData[1:dim], trainData[,dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=1000, ntree=50, k = 5)
  sink()
  
  # obtain posterior samples
  invisible(integrals <- sampleIntegrals(model))
  
  # find the min and max range of y
  ymin <- min(trainY); ymax <- max(trainY)
  
  scaledMean <- (mean(integrals) + 0.5)*(ymax-ymin) + ymin 
  scaledStandardDeviation <- sqrt(var(integrals))*(ymax-ymin)
  meanValue[i] <- scaledMean
  standardDeviation[i] <- scaledStandardDeviation

  # sequential design section, where we build the new training data
  fits <- model$fit$state[[1]]@savedTreeFits
  candidateSet <- randomLHS(1000, dim)

  # predict the values
  fValues <- predict(model, candidateSet)
  
  probability = 1 #uniform probability
  
  expectedValue <- colMeans(fValues*probability)
  
  var <- colVars(fValues)
  index <- sample(which(var==max(var)), 1)
  value <- copeak(candidateSet[index,])
  trainData <- rbind(trainData, c(candidateSet[index,], value))
  
}

  return (list("meanValue"=meanValue, "standardDeviation"=standardDeviation, 
               "trainData" = trainData))

}

mainBARTBQ <- function() 
# main method
# returns prediction as a list
{
  # prepare training data and parameters
  genz <- copeak #select genz function
  trainX <- randomLHS(10, dim) # pick X values from a hypercube (uniform) [a,b]^10
  trainY <- genz(trainX) # test values of y obtained by genz functino
  numNewTraining <- 400
  dim <- 3
  prediction <- BARTBQSequential(dim, trainX, trainY, numNewTraining) 

  return (prediction)

}



#findModel <- function(df, n)
## iterate over nodes
#{
#  iter = 0
#  while (iter <= n) {
#    df <- Findx(df)
#    iter <- iter + 1
#    print (iter)
#  }
#  return (df)
#}


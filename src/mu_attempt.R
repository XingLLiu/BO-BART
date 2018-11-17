library("lhs")
library("dbarts")
library("data.tree")
library("matrixStats")
library("mvtnorm")
library("cubature")
library("truncnorm")
library("mlegp")
library("MASS")
library("base")
namedList<-treatSens:::namedList


#dimension of each observation
dim = 1
#number of observations
observation = 8


#testing continuous Gen z function 
contGenz <- function(x, u = rep(0.5, 1, length(x)), a = rep(5, 1, length(x))){
  
  sum <- sum(a * abs(x-u))
  
  y <- exp(-sum)
  
  return(y)
  
}

#return y for multiple x
contGenzReturn <- function(x){
  
  contGenzY <- NULL
  
  iter <- 1
  
  while (iter <= nrow(x)) {
    
    contGenzY <- rbind(contGenzY, contGenz(x[iter, ]))
    
    iter <- iter + 1
    
  }
  
  return(contGenzY)
  
}


#previous function weighted by probability distribution standard normal
f <- function(x){
  
  return (contGenz(x) * dmvnorm(x))
  
}


#build Tree from characters of tree
buildTree <- function(treeChars){
  
  if (treeChars[1] == ".") {
    return(list(remainder = treeChars[-1]))
    }
  
  splitVar <- as.integer(treeChars[1]) + 1L
  splitIndex <- as.integer(treeChars[2]) + 1L
  
  leftChild <- buildTree(treeChars[-c(1, 2)])
  rightChild <- buildTree(leftChild$remainder)
  leftChild$remainder <- NULL
  remainder <- rightChild$remainder
  rightChild$remainder <- NULL
  
  result <- namedList(splitVar, splitIndex, leftChild, rightChild, remainder)
  leftChild$parent <- result
  rightChild$parent <- result
  
  result

}


#assign probability to each terminal nodes and assign unique name to them
fillProbabilityForNode <- function(oneTree, cutPoints, cut){
  
  if (!is.null(oneTree$leftChild)){
    
    decisionRule <- cutPoints[[oneTree$splitVar]][oneTree$splitIndex]
    
    oneTree$leftChild$probability <- ptruncnorm(decisionRule, a = cut[1, oneTree$splitVar], b = cut[2, oneTree$splitVar], mean = 0, sd = 1)
    
    oneTree$rightChild$probability <- (1 - ptruncnorm(decisionRule, a = cut[1, oneTree$splitVar], b = cut[2, oneTree$splitVar], mean = 0, sd = 1))
    
    cutLeft <- cut
    cutLeft[2, oneTree$splitVar] <- decisionRule
    
    #assign probability to each terminal nodes and assign unique name to them
    fillProbabilityForNode(oneTree$leftChild, cutPoints, cutLeft)
    
    cutRight <- cut
    cutRight[1, oneTree$splitVar] <- decisionRule
    
    #assign probability to each terminal nodes and assign unique name to them
    fillProbabilityForNode(oneTree$rightChild, cutPoints, cutRight)
    
  }else if(is.null(oneTree$probability)){
    
    oneTree$probability <- 1
    
  }
  
  return(oneTree)
}


#drop data set into the tree and assign them to different nodes 
fillObservationsForNode <- function(node, sampler, cutPoints){
  
  node$indices <- seq_len(nrow(model$fit$data@x))
  
  if (!is.null(node$leftChild)) {
    
    goesLeft <- sampler$data@x[node$indices, node$splitVar] <= cutPoints[[node$splitVar]][node$splitIndex]
    node$leftChild$indices  <- node$indices[goesLeft]
    node$rightChild$indices <- node$indices[!goesLeft]
    
    node$leftChild  <- fillObservationsForNode(node$leftChild,  sampler, cutPoints)
    node$rightChild <- fillObservationsForNode(node$rightChild, sampler, cutPoints)
    
  }
  
  return(node)
  
}


#associate terminal value mu
fillPlotInfoForNode <- function(node, sampler, treeFits){
  
  if (!is.null(node$leftChild)) {
    
    node$leftChild  <- fillPlotInfoForNode(node$leftChild,  sampler, treeFits)
    node$rightChild <- fillPlotInfoForNode(node$rightChild, sampler, treeFits)
    
  } else {
    
    node$mu <- treeFits[node$indices[1]]
    
  }
  
  return(node)
  
}


#sum over a single tree's terminal nodes
singleTree <- function(model){
  
  trees <- model$fit$state[[1]]@savedTrees
  fits <- model$fit$state[[1]]@savedTreeFits
  
  cutPoints <- dbarts:::createCutPoints(model$fit)
  cut <- array(c(-Inf, Inf), c(2, dim))
  
  treeList <- buildTree(strsplit(gsub("\\.", "\\. ", trees), " ", fixed = TRUE)[[1]])
  
  #modify tree by the functions written above
  treeList <- fillProbabilityForNode(treeList, cutPoints, cut)
  treeList <- fillObservationsForNode(treeList, model$fit, cutPoints)
  treeList <- fillPlotInfoForNode(treeList, model$fit, model$yhat.train.mean)

  return (treeList)
  
}


#sum terminal value of a single tree
singleTreeSum <- function(treeNum, model, oneTree){
  
  if (!is.null(oneTree$leftChild)){
    
    sum <- singleTreeSum(treeNum, model, oneTree$leftChild) + singleTreeSum(treeNum, model, oneTree$rightChild)
    
  } else {
    
    sum <- oneTree$mu * oneTree$probability
    
  }
  
  return(sum)
  
}


#sum over all the trees in one posterior draws
posteriorSum <- function(drawNum,model){
  
  oneTree <- singleTree(model)
  integral <- 0
  nTree <- ncol(model$fit$state[[1]]@treeFits)
  treeNum <- seq(1, nTree, length.out = nTree)
  
  #extra variables
  var <- list(model, oneTree)
  
  #calculate integration over all trees in the draw by mappy
  integral <- sum(unlist((mapply(singleTreeSum, treeNum, MoreArgs = var, SIMPLIFY = TRUE))))
  
  return (integral)
  
}


#sum over all posterior draws 
sampleIntegrals <- function(model){
  
  nDraw <- dim(model$fit$state[[1]]@savedTreeFits)[3]
  drawNum <- seq(1, nDraw, length.out = nDraw)
  
  #Extra Variables
  var <- list(model)
  
  integrals <- mapply(posteriorSum, drawNum, MoreArgs = var, SIMPLIFY = TRUE)
  
  return (integrals)
  
}

x = seq(-2,2,.02)
y = c(-1,rep(-1,50),rep(0,50),rep(2,50),rep(0,50))
length(x)
length(y)
plot(x,y)

model<-bart(x,y,keeptrees =TRUE,keepevery=20L,nskip=1000,ndpost=2000,ntree=100,k=5)
x.test = seq(-2,2,.01)
y.test = predict(model,x.test)
points(x.test,colMeans(y.test),col="red")

sampleIntegrals(model)

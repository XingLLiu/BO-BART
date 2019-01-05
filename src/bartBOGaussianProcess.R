###
# functions for BART, BO and Gaussian process regression
#
###
library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)

# probabiltity ending up in leaf node

prob2<-function(currentNode){
  prob<-currentNode$probability;
  
  while (!isRoot(currentNode$parent)){
    currentNode<-currentNode$parent
    prob<-prob*currentNode$probability;
  }
  
  return (prob)
}


# Drop data set into the tree and assign them to different nodes 
fillProbabilityForNode <- function(oneTree,cutPoints, cut){
  
  if (!is.null(oneTree$leftChild)){
    
    decisionRule <- cutPoints[[oneTree$splitVar]][oneTree$splitIndex]
    
    oneTree$leftChild$probability  <-(decisionRule - cut[1, oneTree$splitVar]) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
    
    oneTree$rightChild$probability <- (cut[2, oneTree$splitVar] - decisionRule) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
    
    
    cut[,oneTree$splitVar]=c(0,decisionRule)
    
    fillProbabilityForNode(oneTree$leftChild, cutPoints, cut)
    
    cut[, oneTree$splitVar]=c(decisionRule,1)
    
    fillProbabilityForNode(oneTree$rightChild,  cutPoints, cut)
    
  }else if(is.null(oneTree$probability)){
    
    oneTree$probability <- 1
    
  }
  
  return(oneTree)
}

# store probability of getting to terminal node 
# !!?? What does this do @August Shen
terminalProbability<-function(Tree)
{
  terminalNodes = Traverse(Tree,filterFun = isLeaf)
  
  for (i in 1:length(terminalNodes)){
    probability2<-prob2(terminalNodes[[i]])
    terminalNodes[[i]]$terminal_probability<-probability2
  }
  
  return (Tree)
}

# create tree
# @August Shen
getTree <- function(sampler, chainNum, sampleNum, treeNum)
{
  cutPoints <- dbarts:::createCutPoints(sampler)
  
  treeString <-
    if (sampler$control@keepTrees)
      sampler$state[[chainNum]]@savedTrees[treeNum, sampleNum]
  else                           sampler$state[[chainNum]]@trees[treeNum]
  treeFits <-
    if (sampler$control@keepTrees)
      sampler$state[[chainNum]]@savedTreeFits[,treeNum, sampleNum]
  else                           sampler$state[[chainNum]]@treeFits[,treeNum]
  
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
# Sum over a single tree's terminal nodes
SingleTreeSum<-function(treeNum,model,drawNum) {
  
  cutPoints<-dbarts:::createCutPoints(model$fit)
  cut<-array(c(0, 1),c(2,dim))
  
  treeList<-getTree(model$fit,1,drawNum,treeNum)
  
  selectedTree<-FromListSimple(treeList) 
  
  #Modify tree by the functions written above
  selectedTree<-fillProbabilityForNode(selectedTree,cutPoints,cut)
  selectedTree <- terminalProbability(selectedTree)
  
  
  terminalNodeList<-Traverse(selectedTree,filterFun = isLeaf)
  
  #Calculate approximation of integreal in the single tree 
  integral<-0;
  for (node in terminalNodeList){
    
    #We use the mean of prediction value Y's in the terminal node as u
    integral<-integral+node$terminal_probability*(node$mu) 
    
  }
  return (integral)
}


# Sum over all the trees in one posterior draws
PosteriorSum<-function(drawNum,model){
  
  integral<-0;
  nTree<-ncol(model$fit$state[[1]]@treeFits)
  treeNum<-seq(1,nTree,length.out=nTree)
  
  #Extra variables
  var<-list(model,drawNum)
  
  #Calculate integration over all trees in the draw by mapply
  integral<-sum( unlist( mapply(SingleTreeSum,treeNum,MoreArgs=var,SIMPLIFY = TRUE) ) )
  
  return (integral)
}


#Sum over all posterior draws 
sampleIntegrals <- function(model) {
  
  nDraw<- dim(model$fit$state[[1]]@savedTreeFits)[3]
  drawNum<-seq(1,nDraw,length.out=nDraw)
  
  #Extra Variables
  var<-list(model)
  
  integrals<-mapply(PosteriorSum,drawNum,MoreArgs=var,SIMPLIFY = TRUE)
  
  return (integrals)
}




# iterate
Findmodel<-function(df,n){
  iter=0;
  while (iter<=n){
    df<-Findx(df)
    iter<-iter+1
    print (iter)
  }
  return (df)
}

# probability of multivariate normal distribution
p_t_mvnorm<-function(mean, dim){
  return (pmvnorm(lowerx = rep(0,dim),upperx = rep(1,dim),mean, sigma = diag(dim)))
}

# the kernel function of the prior f(x)
# The exp squared covariance function, hence calculating covariance matrix cov
covFunction<-function(x,y) {
  cov<-c()
  for (i in 1:nrow(y)){
    cov[i] <- sum((x-y[i,])^2)
  }
  return (cov)
}

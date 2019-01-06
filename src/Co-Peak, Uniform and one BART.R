library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
namedList<-treatSens:::namedList
dim = 3


f1 <- function(xx, u = rep(0.5, 1, ncol(xx)), a = rep(5, 1, ncol(xx))){
  y<-c();
  for (i in 1:nrow(xx)){
    sum<-0
    for (j in 1:ncol(xx)){
      sum<-sum+a[j]*xx[i, j]
    }
    y[i]<-(1+sum)^(-(ncol(xx)+1))
  }
  return(10000*y)
}

realf1 <- function(xx, u = rep(0.5, 1, length(xx)), a = rep(5, 1, length(xx))){
  sum<-0
  for (i in 1:length(xx)){
    sum<-sum+a[i]*xx[i]
  }
  y<-(1+sum)^(-(length(xx)+1))
  return(10000*y)
}




prob2<-function(currentNode){
  prob<-currentNode$probability;
  
  while (!isRoot(currentNode$parent)){
    currentNode<-currentNode$parent
    prob<-prob*currentNode$probability;
  }
  
  return (prob)
}


#Drop data set into the tree and assign them to different nodes 
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

terminalProbability<-function(Tree)
{
  terminalNodes=Traverse(Tree,filterFun = isLeaf)
  
  for (i in 1:length(terminalNodes)){
    probability2<-prob2(terminalNodes[[i]])
    terminalNodes[[i]]$terminal_probability<-probability2
  }
  
  return (Tree)
}



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
  
  tree
}
#Sum over a single tree's terminal nodes
SingleTreeSum<-function(treeNum,model,drawNum){
  
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


#Sum over all the trees in one posterior draws
PosteriorSum<-function(drawNum,model){
  
  integral<-0;
  nTree<-ncol(model$fit$state[[1]]@treeFits)
  treeNum<-seq(1,nTree,length.out=nTree)
  
  #Extra variables
  var<-list(model,drawNum)
  
  #Calculate integration over all trees in the draw by mappy
  integral<-sum(unlist((mapply(SingleTreeSum,treeNum,MoreArgs=var,SIMPLIFY = TRUE))))
  
  return (integral)
}


#Sum over all posterior draws 
sampleIntegrals<-function(model) {
  
  nDraw<- dim(model$fit$state[[1]]@savedTreeFits)[3]
  drawNum<-seq(1,nDraw,length.out=nDraw)
  
  #Extra Variables
  var<-list(model)
  
  integrals<-mapply(PosteriorSum,drawNum,MoreArgs=var,SIMPLIFY = TRUE)
  
  return (integrals)
}




#iterate
Findmodel<-function(df,n){
  iter=0;
  while (iter<=n){
    df<-Findx(df)
    iter<-iter+1
    print (iter)
  }
  return (df)
}

meanValue<-numeric();
sd<-numeric()
trainX<-randomLHS(10,dim)
trainY<-f1(trainX)
df<-data.frame(trainX,trainY)


for (i in 1:400){
  
  trainY <- df$trainY
  
  
  model<-bart(df[1:dim],trainY,keeptrees = TRUE,keepevery=20L,nskip=1000,ndpost=1000,ntree=50, k = 5)
  
  integrals<-sampleIntegrals(model)
  
  ymin<-min(trainY); ymax<-max(trainY)
  
  scaledMean<-(mean(integrals)+0.5)*(ymax-ymin)+ymin 
  
  sDeviation<-sqrt(var(integrals))*(ymax-ymin)
  
  meanValue<-append(meanValue,scaledMean)
  
  sd<-append(sd,sDeviation)
  
  fits<-model$fit$state[[1]]@savedTreeFits
  candidateSet<-randomLHS(1000,dim)
  fValues<-predict(model,candidateSet);
  #probability=as.vector(dnorm(candidateSet,mean=0,sd=1));
  probability=1
  #expectedValue<-colMeans(fValues%*%diag(probability));
  expectedValue<-colMeans(fValues*probability)
  #var<-colMeans((fValues*probability-expectedValue)^2);
  var<-colVars(fValues);
  index<-sample(which(var==max(var)),1);
  value<-realf1(candidateSet[index,])
  df<-rbind(df,c(candidateSet[index,],value))
  
}



#Monte Carlo of 300 points
meanValue2<-numeric();
sd2<-numeric();

for (i in 1:400){
  
  CandidateSet<- randomLHS(i, dim)
  integration<-mean(f1(CandidateSet))
  meanValue2<-append(meanValue2,integration);
  sd2<-append(sd2,sqrt(var(meanValue2)))
}




p_t_mvnorm<-function(mean){
  
  return (ptmvnorm(lowerx = rep(0,dim),upperx = rep(1,dim),mean, sigma = diag(dim)))
}




#GP
meanValue3<-c()
sd3<-c()

N=10

X<-randomLHS(N, dim);
Y<-f1(X)

K<-matrix(0,nrow=N,ncol=N)

#The exp squared covariance function, hence calculating covariance matrix cov
covFunction<-function(x,y){
  cov<-c()
  for (i in 1:nrow(y)){
<<<<<<< HEAD
    # Can use pairwise difference function to avoid loops
    # Same for the calculation of K below
=======
>>>>>>> parent of e50fa22... [Harrison Zhu 05/01/2019] 1st phase corrections and edits to BO-BART + BO-GP
    cov[i] <- sum((x-y[i,])^2)
  }
  return (cov)
}

for (i in 1:N){
  K[i,]<-covFunction(X[i,],X)
}

z<-c();

for(i in 1:N){
  z[i]<-ptmvnorm(rep(0,dim), rep(1,dim) , mean = X[i,],sigma = diag(dim))[[1]]
}

#loop

for (p in 1:400){
  
  
  meanValue3[p]<-t(z)%*%ginv(K)%*%Y
  
  sd3[p]<-t(z)%*%ginv(K)%*%z
  
  candidateSet<-randomLHS(100,dim)
  
  candidate_Var <- c()
  
  candidate_p <- mapply(p_t_mvnorm,alply(candidateSet,1))
  
  K_prime <- diag(N+p)
  
  K_prime[1:(N+p-1),1:(N+p-1)] <- K
  
  for (i in 1:100){
    
    K_prime[1:(N+p-1),N+p] <- covFunction(candidateSet[i,], X)
    
    K_prime[N+p,1:(N+p-1)] <- covFunction(candidateSet[i,], X)
    
    z[N+p]<-candidate_p[[i]]
    
    candidate_Var[i] <- t(z)%*%ginv(K_prime)%*%z
    
  }
  
  index <- which(candidate_Var == max(candidate_Var))[1]
  
  K_prime[N+p,1:(N+p-1)] <- covFunction(candidateSet[index,], X)
  
  K_prime[1:(N+p-1),N+p] <- covFunction(candidateSet[index,], X)
  
  X<- rbind(X,candidateSet[index,])
  
  Y <- c(Y, realf1(candidateSet[index,]))
  
  K<-K_prime
  
  z[N+p] <- ptmvnorm(rep(0,dim), rep(1,dim) , mean = X[N+p,],sigma = diag(dim))[[1]] 
  
}

real<-adaptIntegrate(realf1,lowerLimit = rep(0,dim),upperLimit = rep(1,dim))
percentageError<-abs((meanValue-real[[1]])/real[[1]])*100
percentageError2<-abs((meanValue2-real[[1]])/real[[1]])*100
percentageError3<-abs((meanValue3-real[[1]])/real[[1]])*100






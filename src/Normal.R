library("lhs")
library("dbarts")
library("data.tree")
library("matrixStats")
library("mvtnorm")
library("cubature")
library("truncnorm")
library("mlegp")
library("MASS")
namedList<-treatSens:::namedList


#Set paramters and dimensions 
dim=1



#Build Tree from characters of tree
buildTree <- function(treeChars){
  if(treeChars[1] == "."){
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

#Assign probability to each terminal nodes and assign unique name to them
namedTree<-function(Tree,base,power){
  nodeList<-Traverse(Tree)
  
  terminalNodes=Traverse(Tree,filterFun = isLeaf)
  
  for(i in 1:length(terminalNodes)){
    probability2<-prob2(terminalNodes[[i]])
    terminalNodes[[i]]$prob2<-probability2
    terminalNodes[[i]]$name<-paste("Terminal",toString(i),sep="")
  }
  
  return (Tree)
}

#Calculate terminal node probability by multiplying along the path
prob2<-function(currentNode){
  prob<-currentNode$probability;
  
  while(!isRoot(currentNode$parent)){
    currentNode<-currentNode$parent
    prob<-prob*currentNode$probability;
  }
  
  return (prob)
}


#Drop data set into the tree and assign them to different nodes 
passData <- function(oneTree, cutPoints,dataPoints, cut){
  
  if(!is.null(oneTree$leftChild)){
    
    decisionRule <- cutPoints[[oneTree$splitVar]][oneTree$splitIndex]
    
    oneTree$leftChild$probability  <-ptruncnorm(decisionRule,a=cut[1,oneTree$splitVar],b=cut[2,oneTree$splitVar],mean=0,sd=1)
    oneTree$leftChild$xData <-subset(dataPoints,dataPoints[oneTree$splitVar]<=decisionRule)
    
    oneTree$rightChild$probability <-(1-ptruncnorm(decisionRule,a=cut[1,oneTree$splitVar],b=cut[2,oneTree$splitVar],mean=0,sd=1))
    oneTree$rightChild$xData <- subset(dataPoints,dataPoints[oneTree$splitVar]>decisionRule)
    
    cutLeft<-cut
    cutLeft[2,oneTree$splitVar]<-decisionRule
    
    passData(oneTree$leftChild, cutPoints, oneTree$leftChild$xData,cutLeft)
    
    cutRight<-cut
    cutRight[1,oneTree$splitVar]<-decisionRule
    
    passData(oneTree$rightChild,  cutPoints, oneTree$rightChild$xData,cutRight)
    
  } else if(is.null(oneTree$probability)){
    
    oneTree$probability <- 1
    oneTree$xData <- dataPoints
    
  }
}

#Predicted Y's for the data X's dropped in each terminal nodes
Yprediction<-function(tree,model,treeNum,drawNum,trainX){
  
  fits<-model$fit$state[[1]]@savedTreeFits
  terminal<-Traverse(tree,filterFun = isLeaf);
  
  for(node in terminal){
    xTrain<-node$xData;
    predictedY<-fits[as.integer(rownames(xTrain)),treeNum,drawNum]
    data<-data.frame(xTrain,predictedY);
    node$data<-data;
  }
}


#Sum over a single tree's terminal nodes
SingleTreeSum<-function(treeNum,model,trainX,drawNum,base,power){
  
  trees<-model$fit$state[[1]]@savedTrees
  fits<-model$fit$state[[1]]@savedTreeFits
  
  cutPoints<-dbarts:::createCutPoints(model$fit)
  cut<-array(c(-Inf,Inf),c(2,dim))
  
  treeList<-buildTree(strsplit(gsub("\\.", "\\. ", trees[treeNum,drawNum]), " ", fixed = TRUE)[[1]])
  selectedTree<-FromListSimple(treeList) 
  
  #Modify tree by the functions written above
  passData(selectedTree,cutPoints,trainX,cut)
  Yprediction(selectedTree,model,treeNum,drawNum,trainX)
  namedTree(selectedTree,base,power)
  
  terminalNodeList<-Traverse(selectedTree,filterFun = isLeaf)
  
  #Calculate approximation of integreal in the single tree 
  integral<-0;
  for(node in terminalNodeList){
    
    #We use the mean of prediction value Y's in the terminal node as u
    integral<-integral+node$prob2*(mean(node$data$predictedY)) 
    
  }
  return (integral)
}


#Sum over all the trees in one posterior draws
PosteriorSum<-function(drawNum,model,trainX,base,power){
  
  integral<-0;
  nTree<-ncol(model$fit$state[[1]]@treeFits)
  treeNum<-seq(1,nTree,length.out=nTree)
  
  #Extra variables
  var<-list(model,trainX,drawNum,base,power)
  
  #Calculate integration over all trees in the draw by mappy
  integral<-sum(unlist((mapply(SingleTreeSum,treeNum,MoreArgs=var,SIMPLIFY = TRUE))))
  
  return (integral)
}


#Sum over all posterior draws 
sampleIntegrals<-function(model,trainX,base,power) {
  
  nDraw<- dim(model$fit$state[[1]]@savedTreeFits)[3]
  drawNum<-seq(1,nDraw,length.out=nDraw)
  
  #Extra Variables
  var<-list(model,trainX,base,power)
  
  integrals<-mapply(PosteriorSum,drawNum,MoreArgs=var,SIMPLIFY = TRUE)
  
  return (integrals)
}

#Find the next inqury point using maximum variance
Findx<-function(data)
{
  trainX<-data[1:dim];
  trainY<-data$trainY
  model<-bart(trainX,trainY,keeptrees = TRUE,keepevery=20L,nskip=1000,ndpost=2000,ntree=50);
  fits<-model$fit$state[[1]]@savedTreeFits
  
  #Choose a candidate set, find the point of maximum posterior variance
  #Add that point to the dataframe
  candidateSet<-rmvnorm(1000, mean = rep(0, dim))
  fValues <- predict(model, candidateSet)
  probability <- dmvnorm(candidateSet, mean = rep(0,dim))
  var <- colVars(fValues %*% diag(probability))
  index <- sample( which(var == max(var)), 1 )
  value <- f1(candidateSet[index,])
  data <- rbind(data, c(candidateSet[index,], value))
  
  return (data)
}

#How many points we want to add in the dataframe, controlled by n

Findmodel <- function(data,n)
{
  iter=0;
  while (iter < n)
  {
    data <- Findx(data)
    iter <- iter + 1
    print (iter)
  }
  return (data)
}

error <- list()
variance <- list()
prediction <- list()
realValue <- numeric()

#The function to be integrated, the case here is a discontinuous one
f1<-function(x)
{
  y <- c();
  for( i in 1:length(x) )
  {
    if(x[i] < 0) {
      y[i] <- x[i]^2;
    } else {
      y[i] <- cos(x[i]) + x
    }
  }
  return (y)
}

#Previous function weighted by probability distribution standard normal
f <- function(x)
{
  return ( f1(x) * dmvnorm(x) )
}


#Monte Carlo of 300 points
meanValue2 <- numeric()
sd2 <- numeric()

for (i in 1:300)
{
  CandidateSet <- rmvnorm(1*i, mean=rep(0,dim))
  integration <- mean(f1(CandidateSet))
  meanValue2 <- append(meanValue2, integration)
  sd2 <- append(sd2, sqrt(var(meanValue2)))
}




#BART Quadrature

#numerics to store standard deviation and mean value(approxiamtion)
meanValue <- numeric();
sd <- numeric()

#Initialize dataset, size of dataset depends on dimensions
trainX <- rmvnorm(5, mean=rep(0,dim))
trainY <- f1(trainX)
data <- data.frame(trainX, trainY)

#Run BART, keep finding approximation as new data point added in data
#store approxiamtion in mean value, posterior sd in standard deviation

for (i in 1:400){
  
  data <- Findmodel(data, 1);
  model <- bart( data[1:dim], data$trainY,keeptrees =TRUE,keepevery=20L,nskip=1000,ndpost=2000,ntree=50,k=5 )
  integrals <- sampleIntegrals(model, data[1:dim], 0.95, 2)
  
  #Scale back the approxiamtion and standard deviation
  ymin <- min(data$trainY); ymax<-max(data$trainY)
  sDeviation <- sqrt(var(integrals)) * (ymax - ymin)
  scaledMean <- (mean(integrals) + 0.5) * (ymax - ymin) + ymin 
  
  meanValue <- append(meanValue, scaledMean)
  sd <- append(sd, sDeviation)
  
}


#Bayesian Quadrature
#Calculate the approxiamtion and posterior variance in GP
#using analytical solutions provided in the paper
meanValue3<-numeric();
sd3<-numeric();

#Find 300 approxiamtion in GP, first 10 points are for initialization
for (p in 10:310){
  
  N=1*p;
  X<-rmvnorm(N,mean=rep(0,dim));
  Y<-f1(X)
  
  #Adjust the parameters using mlegp function  (doubt!?)
  model<-mlegp(X,Y)
  Beta<-model$beta
  sigma<-model$sig2
  
  K<-matrix(0,nrow=N,ncol=N)
  
#The exp squared covariance function, hence calculating covariance matrix cov
covFunction<-function(x,y){
  cov<-1
  for (i in 1:dim){
    cov<-cov*exp(-Beta[i]*sum((x[i]-y[i])^2))
  }
  return (cov)
}

for (i in 1:N){
  for (j in 1:N){
    cov<-sigma*covFunction(X[i,],X[j,])
    K[i,j]<-cov
  }
}
  
#calculate approximation and posterior variance analyticall using formula
#in the paper
b<-0;
B<-diag(dim);
a<-X;
A<-diag(1/(2*Beta),dim)
z<-c();
for(i in 1:N){
  z[i]<-sigma*(det(solve(A)%*%B+diag(dim))^-0.5)*exp(-0.5*(a[i,]+b)%*%ginv(A+B)%*%(a[i,]-b))
}

meanValue3[p]<-t(z)%*%ginv(K)%*%Y
sd3[p]<-sigma*det(2*solve(A)%*%B+diag(dim))^(-0.5)-t(z)%*%ginv(K)%*%z



#This is to find real value, only work in low dimensions (up to 3);

real<-adaptIntegrate(f,lowerLimit = rep(-5,dim),upperLimit = rep(5,dim))
percentageError<-abs((meanValue-real[[1]])/real[[1]])*100

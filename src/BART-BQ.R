
#Set paramters and dimensions 
dim=3


#Build Tree from characters of tree
buildTree <- function(treeChars)
{
  if (treeChars[1] == ".") return(list(remainder = treeChars[-1]))
  
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
namedTree<-function(Tree,base,power)
{
  nodeList<-Traverse(Tree)
  
  terminalNodes=Traverse(Tree,filterFun = isLeaf)
  
  for (i in 1:length(terminalNodes)){
    probability2<-prob2(terminalNodes[[i]])
    terminalNodes[[i]]$prob2<-probability2
    terminalNodes[[i]]$name<-paste("Terminal",toString(i),sep="")
  }
  
  return (Tree)
}

#Calculate terminal node probability by multiplying along the path
prob2<-function(currentNode){
  prob<-currentNode$probability;
  
  while (!isRoot(currentNode$parent)){
    currentNode<-currentNode$parent
    prob<-prob*currentNode$probability;
  }
  
  return (prob)
}


#Drop data set into the tree and assign them to different nodes 
passData <- function(oneTree, cutPoints,dataPoints, cut){
  
  if (!is.null(oneTree$leftChild)){
    
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
    
  }else if(is.null(oneTree$probability)){
    
    oneTree$probability <- 1
    oneTree$xData <- dataPoints
    
  }
}

#Predicted Y's for the data X's dropped in each terminal nodes
Yprediction<-function(tree,model,treeNum,drawNum,trainX){
  
  fits<-model$fit$state[[1]]@savedTreeFits
  terminal<-Traverse(tree,filterFun = isLeaf);
  
  for (node in terminal){
    xTrain<-node$xData;
    predictedY<-fits[as.integer(rownames(xTrain)),treeNum,drawNum]
    df<-data.frame(xTrain,predictedY);
    node$data<-df;
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
  for (node in terminalNodeList){
    
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
  drawNum<-seq(1, ,length.out=nDraw)
  
  #Extra Variables
  var<-list(model,trainX,base,power)
  
  integrals<-mapply(PosteriorSum,drawNum,MoreArgs=var,SIMPLIFY = TRUE)
  
  return (integrals)
}



rescale = function(x) (x-min(x))/(max(x)-min(x)) - .5

#Find the next inqury point using maximum variance
Findx<-function(df){
  
  trainX<-df[1:dim];
  trainY<-df$trainY
  model<-bart(trainX,trainY,keeptrees = TRUE,keepevery=20L,nskip=1000,ndpost=2000,ntree=50,verbose=F);
  fits<-model$fit$state[[1]]@savedTreeFits
  
  #Choos a candidate set, find the point of maximum posterior variance
  #Add that point to the dataframe
  candidateSet<-rmvnorm(1000,mean=rep(0,dim))
  fValues<-predict(model,candidateSet);
  probability=dmvnorm(candidateSet,mean=rep(0,dim));
  var<-colVars(fValues%*%diag(probability));
  index<-sample(which(var==max(var)),1);
  value<-realf1(candidateSet[index,])
  df<-rbind(df,c(candidateSet[index,],value))
  
  return (df)
  
}

#How many points we want to add in the dataframe, controlled by n

Findmodel<-function(df,n){
  iter=0;
  while (iter<n){
    df<-Findx(df)
    iter<-iter+1
    print (iter)
  }
  return (df)
}

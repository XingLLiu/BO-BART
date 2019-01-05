# prepare training data
meanValue <- numeric() 
standardDeviation <- numeric() # standard deviation
trainX <- randomLHS(10, dim) # pick X values from a hypercube (uniform) [a,b]^10
trainY <- copeak(trainX) # test values of y obtained by genz functino
trainData <- data.frame(trainX, trainY)

# generate extra training data using the scheme (see pdf)
for (i in 1:epochs) {
  
  trainY <- trainData$trainY
  print(c("Epoch=", i))
  
  model<-bart(trainData[1:dim],trainY,keeptrees = TRUE,keepevery=20L,nskip=1000,ndpost=1000,ntree=50, k = 5)
  
  # obtain posterior samples
  integrals<-sampleIntegrals(model)
  
  # find the min and max range of y
  ymin<-min(trainY); ymax<-max(trainY)
  
  scaledMean <- (mean(integrals)+0.5)*(ymax-ymin)+ymin 
  sDeviation <- sqrt(var(integrals))*(ymax-ymin)
  meanValue<-append(meanValue,scaledMean)
  standardDeviation<-append(standardDeviation,sDeviation)

  fits<-model$fit$state[[1]]@savedTreeFits
  candidateSet<-randomLHS(1000,dim)

  # predict the values
  fValues<-predict(model,candidateSet)
  
  #probability=as.vector(dnorm(candidateSet,mean=0,standardDeviation=1));
  probability = 1
  
  #expectedValue<-colMeans(fValues%*%diag(probability));
  expectedValue<-colMeans(fValues*probability)
  
  #var<-colMeans((fValues*probability-expectedValue)^2)
  var<-colVars(fValues)
  index<-sample(which(var==max(var)),1)
  value<-copeak(candidateSet[index,])
  trainData<-rbind(trainData,c(candidateSet[index,],value))
  
}

# model is the BART we use for predictions
# trainData is the training set we will use to train our model



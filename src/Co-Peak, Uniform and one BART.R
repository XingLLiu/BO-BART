source("src/bartBOGaussianProcess.R")
source("src/references/cornerPeakFamily.R")
library(treatSens)
library(mvtnorm)
library(MASS)
library(cubature)
# main method
# training parameters
namedList<-treatSens:::namedList
dim <- 3
epochs <- 10

# prepare training data
meanValue <- numeric() 
standardDeviation <- numeric() # standard deviation
trainX <- randomLHS(10, dim) # pick X values from a hypercube (uniform) [a,b]^10
trainY <- copeak(trainX) # test values of y obtained by genz functino
trainData <- data.frame(trainX, trainY)

# train
for (i in 1:epochs) {
  
  trainY <- trainData$trainY
  print(c("Epoch=", i))
  
  model<-bart(trainData[1:dim],trainY,keeptrees = TRUE,keepevery=20L,nskip=1000,ndpost=1000,ntree=50, k = 5)
  integrals<-sampleIntegrals(model)
  
  ymin<-min(trainY); ymax<-max(trainY)
  
  scaledMean <- (mean(integrals)+0.5)*(ymax-ymin)+ymin 
  sDeviation <- sqrt(var(integrals))*(ymax-ymin)
  meanValue<-append(meanValue,scaledMean)
  standardDeviation<-append(standardDeviation,sDeviation)
  
  fits<-model$fit$state[[1]]@savedTreeFits
  candidateSet<-randomLHS(1000,dim)
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



# Monte Carlo sampling of 300 points
# @August
meanValue2<-numeric();
standardDeviation2<-numeric();

for (i in 1:epochs){
  
  CandidateSet<- randomLHS(i, dim)
  integration<-mean(copeak(CandidateSet))
  meanValue2<-append(meanValue2,integration);
  standardDeviation2<-append(standardDeviation2,sqrt(var(meanValue2)))
}


# GP
meanValue3 <- c()
standardDeviation3 <- c()

N=10

X<-randomLHS(N, dim)
Y<-copeak(X)

K<-matrix(0,nrow=N,ncol=N)



for (i in 1:N){
  K[i,]<-covFunction(X[i,],X)
}

# we use bayesian quadrature
# take samples from Gaussian process integral z[i]
# see pdf on Bayesian quadrature
#!!can be vectorised using mapply
z<-c();
for(i in 1:N) {
  z[i]<-pmvnorm(rep(0,dim), rep(1,dim) , mean = X[i,], sigma = diag(dim))[[1]] * (2*pi)^(dim/2) # add in extra term obtained by integration
}
meanValue3[1]<-t(z) %*% ginv(K) %*% Y

standardDeviation3[1]<-t(z)%*%ginv(K)%*%z #not quite right, missed out first term

# train
for (p in 1:epochs){
  
  
  candidateSet<-randomLHS(100,dim)
  
  candidate_Var <- c()
  
  candidate_p <- c()
  
  for(i in 1:100) {
    candidate_p[i]<-pmvnorm(rep(0,dim), rep(1,dim) , mean = candidateSet[i,], sigma = diag(dim))[[1]] * (2*pi)^(dim/2) 
    # add in extra term obtained by integration
  }
  
  K_prime <- diag(N+p)
  
  K_prime[1:(N+p-1), 1:(N+p-1)] <- K
  
  for (i in 1:100) {
    
    K_prime[1:(N+p-1),(N+p)] <- covFunction(candidateSet[i,], X)
    
    K_prime[(N+p),1:(N+p-1)] <- covFunction(candidateSet[i,], X)
    
    z[N+p] <- candidate_p[i]
    
    candidate_Var[i] <- t(z)%*%ginv(K_prime)%*%z # where is the integral here
    
  }
  
  index <- which(candidate_Var == max(candidate_Var))[1]
  
  K_prime[N+p,1:(N+p-1)] <- covFunction(candidateSet[index,], X)
  
  K_prime[1:(N+p-1),N+p] <- covFunction(candidateSet[index,], X)
  
  X<- rbind(X,candidateSet[index,])
  
  Y <- c(Y, copeak(candidateSet[index,]))
  
  K<-K_prime
  
  # add in extra term obtained by integration
  z[N+p] <- pmvnorm(rep(0,dim), rep(1,dim) , mean = X[N+p,],sigma = diag(dim))[[1]]* (2*pi)^(dim/2)
  
  meanValue3[p+1]<-t(z) %*% ginv(K) %*% Y
  
  standardDeviation3[p+1]<-t(z)%*%ginv(K)%*%z #not quite right, missed out first term
  
}

# cubature contains adaptIntegrate function
real <- adaptIntegrate(copeak,lowerLimit = rep(0,dim),upperLimit = rep(1,dim))
percentageError<-abs((meanValue-real[[1]])/real[[1]])*100
percentageError2<-abs((meanValue2-real[[1]])/real[[1]])*100
percentageError3<-abs((meanValue3-real[[1]])/real[[1]])*100

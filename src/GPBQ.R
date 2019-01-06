############
#Bayesian Quadrature with Gaussian Process
############
library(mvtnorm)

gaussianKernel <- function(X, xPrime) 
  # The exp squared covariance function, hence calculating covariance matrix cov
  # input:
  #     X: Initial design matrix
  #     xPrime: a covariate. 
  # output:
  #     K: covariance matrix i.e. cov(X, xPrime)
{
  K <- c()
  
  for (i in 1:nrow(X)) {
    K[i] <- exp( -0.5 * norm(xPrime - X[i,], type = "2") ^ 2 )
  }

  return (K)
}

#covariance matrix with exponential Gaussian kernel
Gaussian_kernel <- function(x){
  
  K <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  
  for (i in 1:nrow(x)){
    
    K[i, ] <- kernel(x[i, ], x)
    
  }
  
  return(K)
  
}



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
for (p in 1:epochs) {
  
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
  z[N+p] <- pmvnorm(rep(0,dim), rep(1,dim), mean = X[N+p,], sigma = diag(dim))[[1]] * (2*pi)^(dim/2)
  
  meanValue3[p+1]<-t(z) %*% ginv(K) %*% Y
  
  standardDeviation3[p+1]<-t(z)%*%ginv(K)%*%z #not quite right, missed out first term
  
}



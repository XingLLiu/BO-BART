############
#Bayesian Quadrature with Gaussian Process
############
library(mvtnorm)

gaussianKernel <- function(xPrime, X) 
  # The exp squared covariance function, hence calculating covariance matrix cov
  # input:
  #     X: Initial design matrix
  #     xPrime: a covariate. 
  # output:
  #     K: covariance i.e. cov(X, xPrime)
{
  K <- c()
  
  for (i in 1:nrow(X)) {
    K[i] <- exp( -0.5 * norm(xPrime - X[i,], type = "2") ^ 2 )
  }

  return (K)
}


computeGPBQ <- function(dim, epochs, N=10, FUN) 
# method for computation of the integration
# includes query sequential design
# input:
#     dim:
#     epochs: number of training instances
#     N:
{
  #define genz function
  genz <- FUN
  meanValueGP <- c()
  standardDeviationGP <- c()

  X <- randomLHS(N, dim)
  Y <- genz(X)

  K <- matrix(0,nrow=N,ncol=N)

  # compute kernel matrix for starting dataset
  for (i in 1:N) {
    K[i,] <- gaussianKernel(X[i,], X)
  }

  
  z<-c();
  for(i in 1:N) {
    z[i] <- pmvnorm(rep(0,dim), rep(1,dim) , mean = X[i,], sigma = diag(dim))[[1]] * (2*pi)^(dim/2) # add in extra term obtained by integration
  }
  meanValueGP[1] <- t(z) %*% ginv(K) %*% Y

  standardDeviationGP[1] <- t(z)%*%ginv(K)%*%z #not quite right, missed out first term


  # train
  for (p in 1:epochs) {
  print(c("Epochs no.:", p))
  candidateSet <- randomLHS(100,dim)
  
  candidate_Var <- c()
  
  candidate_p <- c()
  
  for(i in 1:100) {
    candidate_p[i] <- pmvnorm(rep(0, dim), rep(1, dim) , mean = candidateSet[i,], sigma = diag(dim))[[1]] * (2*pi)^(dim/2) 
    # add in extra term obtained by integration
  }
  
  K_prime <- diag(N+p)
  
  K_prime[1:(N+p-1), 1:(N+p-1)] <- K
  
  for (i in 1:100) {
    
    
    K_prime[1:(N+p-1),(N+p)] <- gaussianKernel(candidateSet[i,], X)
    
    K_prime[(N+p),1:(N+p-1)] <- gaussianKernel(candidateSet[i,], X)
    
    z[N+p] <- candidate_p[i]
    
    candidate_Var[i] <- t(z)%*%ginv(K_prime)%*%z # where is the integral here
    
  }
  
  index <- which(candidate_Var == max(candidate_Var))[1]
  
  K_prime[N+p,1:(N+p-1)] <- gaussianKernel(candidateSet[index,], X)
  
  K_prime[1:(N+p-1),N+p] <- gaussianKernel(candidateSet[index,], X)
  
  X <- rbind(X,candidateSet[index,])
  
  Y <- c(Y, genz( matrix( candidateSet[index,], ncol = length(candidateSet[index,])) ) )
  
  K <- K_prime
  
  # add in extra term obtained by integration
  z[N+p] <- pmvnorm(rep(0,dim), rep(1,dim), mean = X[N+p,], sigma = diag(dim))[[1]] * (2*pi)^(dim/2)
  
  meanValueGP[p+1] <- t(z) %*% ginv(K) %*% as.matrix(Y)
  
  standardDeviationGP[p+1] <- t(z)%*%ginv(K)%*%z #not quite right, missed out first term
  
  }

  return (list("meanValueGP" = meanValueGP, "standardDeviationGP" = standardDeviationGP))

}



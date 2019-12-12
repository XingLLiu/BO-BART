#Bayesian Quadrature with Gaussian Process
############
library(mvtnorm)
library(MASS)
library(kernlab)
gaussianKernel <- function(xPrime, X, lengthscale = 1) 
  
  #'Standard Kernel Function in GP
  #' 
  #'@description This function calculates the covariance of xPrime with respect
  #'to each row of X
  #' 
  #'@param xprime Numeric Array; 
  #'@param X Matrix
  #'@param lengthscale Real; parameter in standard kernel
  #'
  #'@return K; Covariance matrix between xprime and X
  
{
  K <- c()
  
  for (i in 1:nrow(X)) {
    K[i] <- exp( -0.5 * sum((xPrime - X[i,]))^2 / lengthscale^2)
  }

  return (K)
}

rescale <- function(x) {x * attr(x, 'scaled:scale') + attr(x, 'scaled:center')}

computeGPBQ <- function(X, Y_unscaled, dim, epochs, FUN, lengthscale=1, sequential) 
  
  #'Gaussian Process with Bayesian Quadrature
  #' 
  #'@description This function calculates the approxiamtion of integration using
  #'Gaussian Process, Bayesian Quadrature and Sequential Design
  #' 
  #'@param dim Integer; Dimension of input X 
  #'@param epochs Integer; Number of new data points
  #'@param FUN Function; The function to be integrated
  #'@param lengthscale Integer; The parameter in standard kernel
  #'
  #'@return List; A list containing meanValue (apprimation) and variance of GP method

{
  #define genz function
  genz <- FUN
  meanValueGP <- c()
  varianceGP <- c()
   
  N <- dim(X)[1]
  Y <- scale(Y_unscaled)

  K <- matrix(0,nrow=N,ncol=N)
  jitter = 1e-6

  # compute kernel matrix for starting dataset
  #for (i in 1:N) {
  #  K[i,] <- gaussianKernel(X[i,], X)
  #}
  K = kernelMatrix(rbfdot(.5/lengthscale^2), X)
  
  # compute the variance
  k = function(arg) {
    return(exp(-.5 * sum((arg[1]-arg[2])^2)/lengthscale^2))
  }
  int.points <- randomLHS(500000, 2)
  vv = vapply(1:nrow(int.points), function(i) k(int.points[i,]),0)
  var.firstterm = mean(vv)^dim
  
  z<-c()
  for(i in 1:N) {
    z[i] <- pmvnorm(rep(0,dim), rep(1,dim) , mean = X[i,], sigma = diag(lengthscale^2, dim))[[1]] * (2*pi*lengthscale^2)^(dim/2) # add in extra term obtained by integration
  }
  meanValueGP[1] <- t(z) %*% solve(K + diag(jitter, N) , rescale(Y))
  varianceGP[1] <- var.firstterm - t(z)%*% solve(K + diag(jitter, N) , z) #not quite right, missed out first term

  # train
  for (p in 1:epochs) {
   
    print(paste("GPBQ: Epoch =", p))
    candidateSetNum <- 100
    candidateSet <- randomLHS(candidateSetNum,dim)
    
    candidate_Var <- c()
    
    candidate_p <- c()
    for(i in 1:candidateSetNum) {
      candidate_p[i] <- pmvnorm(rep(0, dim), rep(1, dim) , mean = candidateSet[i,], sigma = diag(lengthscale^2, dim))[[1]] * (2*pi*lengthscale^2)^(dim/2) 
      # add in extra term obtained by integration
    }
    
    K_prime <- diag(N+p)
    
    K_prime[1:(N+p-1), 1:(N+p-1)] <- K
    for (i in 1:candidateSetNum) {

      kernel_new_entries <- kernelMatrix(rbfdot(.5/lengthscale^2), matrix(candidateSet[i,], nrow = 1), X)


      K_prime[1:(N+p-1),(N+p)] <- kernel_new_entries
      K_prime[(N+p),1:(N+p-1)] <- kernel_new_entries
      
      z[N+p] <- candidate_p[i]
      
      if (sequential){
        candidate_Var[i] <- t(z)%*%solve(K_prime + diag(jitter,nrow(K_prime)), z) # where is the integral here
      }
    }
    
    if (sequential){
      index <- which(candidate_Var == max(candidate_Var))[1]
    }
    else {
      index <- sample(1:candidateSetNum, 1)
    }
    
    kernel_new_entry <- kernelMatrix(rbfdot(.5/lengthscale^2), matrix(candidateSet[index,], nrow=1), X)
    
    K_prime[N+p,1:(N+p-1)] <- kernel_new_entry
    K_prime[1:(N+p-1),N+p] <- kernel_new_entry
    X <- rbind(X,candidateSet[index,])
    additionalResponse <- as.matrix( t(candidateSet[index,]), ncol = length(candidateSet[index,]) )

    Y_unscaled <- c(Y_unscaled, genz(additionalResponse))
    Y <- scale(Y_unscaled)
    K <- K_prime
    
    # add in extra term obtained by integration
    z[N+p] <- pmvnorm(rep(0,dim), rep(1,dim), mean = X[N+p,], sigma = diag(lengthscale^2, dim))[[1]] * (2*pi*lengthscale^2)^(dim/2)
    meanValueGP[p+1] <- t(z) %*% solve(K + diag(jitter,nrow(K)), as.matrix(rescale(Y)))
    varianceGP[p+1] <- var.firstterm - t(z)%*%solve(K + diag(jitter,nrow(K)), z) #not quite right, missed out first term
    
  }

  return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K))
}
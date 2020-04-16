#Bayesian Quadrature with Gaussian Process
############
library(mvtnorm)
library(MASS)
library(kernlab)
library(MCMCglmm)
maternKernelWrapper <- function(lengthscale = 1, sigma = 1) {
  maternKernel <- function(x, y) 
    #'Matern 3/2 Kernel Function in GP
    #' 
    #'@description This function calculates the covariance k(x,y)
    #' 
    #'@param x
    #'@param y
    #'@param lengthscale Real; parameter in standard kernel
    #'@param sigma Real; parameter in standard kernel
    #'
    #'@return K; Covariance matrix between xprime and X
  {
    d <- sqrt(sum((x - y)^2))
    k <- sigma^2 * (1 + sqrt(3) * d / lengthscale) * exp(- sqrt(3) * d / lengthscale)
    return (k)
  }
  return (maternKernel)
}

rescale <- function(x) {x * attr(x, 'scaled:scale') + attr(x, 'scaled:center')}

computeGPBQEmpirical <- function(X, Y, dim, candidateSet, candidateY, epochs, kernel="rbf", lengthscale, sequential=TRUE) 
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
  meanValueGP <- c()
  varianceGP <- c()
   
  N <- dim(X)[1]

  K <- matrix(0,nrow=N,ncol=N)
  jitter = 1e-7

  if (kernel == "rbf") {
    kernel <- rbfdot(.5/lengthscale^2)
  }
  else if (kernel == "matern32") {
     kernel <- maternKernelWrapper(lengthscale)
  }  
  
  K = kernelMatrix(kernel, X)
  # compute the variance
  var.firstterm <- sum(colSums(K))/(nrow(X)^2)
  z <- matrix(colSums(K)/nrow(X))
  covInverse <- chol2inv(chol(K + diag(jitter, nrow(K))))
  meanValueGP[1] <- t(z) %*% covInverse %*% Y
  tmp <- t(z)%*% covInverse %*% z 
  varianceGP[1] <- var.firstterm - tmp
  cat(var.firstterm, tmp,"\n")

  # train
  if (epochs == 1){
    return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K))
  }
  for (p in 2:epochs) {
   
    print(paste("GPBQ: Epoch =", p))
    K_prime <- diag(N+p-1)
    K_prime[1:(N+p-2), 1:(N+p-2)] <- K
    
    K_star_star <- kernelMatrix(kernel, candidateSet)
    K_star <- kernelMatrix(kernel, candidateSet, X)

    candidate_Var <- diag(K_star_star - K_star %*% solve(K + diag(jitter, nrow(K)), t(K_star)))
    
    index <- which(candidate_Var == max(candidate_Var))[1]
    
    kernel_new_entry <- kernelMatrix(kernel, matrix(candidateSet[index,], nrow=1), X)
    
    K_prime[N+p-1,1:(N+p-2)] <- kernel_new_entry
    K_prime[1:(N+p-2),N+p-1] <- kernel_new_entry
    X <- rbind(X,candidateSet[index,])

    Y <- c(Y, candidateY[index])
    K <- K_prime
    
    # add in extra term obtained by integration
    z <- matrix(colSums(K))/nrow(X)
    covInverse <- chol2inv(chol(K + diag(jitter, N+p-1)))
    meanValueGP[p] <- t(z) %*% covInverse %*% Y
    var.firstterm <- sum(colSums(K))/(nrow(X)^2)
    varianceGP[p] <- var.firstterm - t(z) %*% covInverse %*% z
  }

  return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K))
}
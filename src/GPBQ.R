#Bayesian Quadrature with Gaussian Process
############
library(mvtnorm)
library(MASS)
library(kernlab)
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

computeGPBQ <- function(X, Y, dim, epochs, kernel="rbf", FUN, lengthscale=1, sequential=TRUE, measure) 
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

  K <- matrix(0,nrow=N,ncol=N)
  jitter = 1e-8

  if (kernel == "rbf") {
    kernel <- rbfdot(.5/lengthscale^2)
  }
  else if (kernel == "matern32") {
     kernel <- maternKernelWrapper(lengthscale)
  }  
  
  K = kernelMatrix(kernel, X)
  # compute the variance
  if (measure == "uniform"){
    int.points.1 <- replicate(dim, runif(5000))
    int.points.2 <- replicate(dim, runif(5000))
  } else if (measure == "gaussian") {
    int.points.1 <- replicate(dim, rtnorm(5000, mean = 0.5, lower=0, upper=1))
    int.points.2 <- replicate(dim, rtnorm(5000, mean = 0.5, lower=0, upper=1))
  }
  cov <- kernelMatrix(kernel, int.points.1, int.points.2)
  var.firstterm <- mean(cov[upper.tri(cov)])
  cov <- kernelMatrix(kernel, int.points.1, X)
  z <- colMeans(cov)
  meanValueGP[1] <- t(z) %*% solve(K + diag(jitter, N) , Y)
  tmp <- t(z)%*% solve(K + diag(jitter, N) , z) 
  varianceGP[1] <- var.firstterm - tmp
  cat(var.firstterm, tmp,"\n")

  # train
  if (epochs == 0){
    return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K))
  }
  for (p in 1:epochs) {
   
    print(paste("GPBQ: Epoch =", p))
    candidateSetNum <- 100
	  candidateSet <- replicate(dim, runif(candidateSetNum))
    K_prime <- diag(N+p)
    K_prime[1:(N+p-1), 1:(N+p-1)] <- K
    


    
    if (sequential){
      K_star_star <- kernelMatrix(kernel, candidateSet)
      K_star <- kernelMatrix(kernel, candidateSet, X)
      candidate_Var <- diag(K_star_star - K_star %*% solve(K + diag(jitter, nrow(K)), t(K_star)))
      index <- which(candidate_Var == max(candidate_Var))[1]
    }
    else {
      index <- sample(1:candidateSetNum, 1)
    }
    
    kernel_new_entry <- kernelMatrix(kernel, matrix(candidateSet[index,], nrow=1), X)
    
    K_prime[N+p,1:(N+p-1)] <- kernel_new_entry
    K_prime[1:(N+p-1),N+p] <- kernel_new_entry
    X <- rbind(X,candidateSet[index,])
    additionalResponse <- as.matrix( t(candidateSet[index,]), ncol = length(candidateSet[index,]) )

    Y <- c(Y, genz(additionalResponse))
    K <- K_prime
    
    # add in extra term obtained by integration
    cov <- kernelMatrix(kernel, int.points.1, X)
    z <- colMeans(cov)
    meanValueGP[p+1] <- t(z) %*% solve(K + diag(jitter, N+p) , Y)
    varianceGP[p+1] <- var.firstterm - t(z)%*% solve(K + diag(jitter, N+p) , z) #not quite right, missed out first term
  }

  return (list("meanValueGP" = meanValueGP, "varianceGP" = varianceGP, "X" = X, "Y" = Y, "K" = K))
}

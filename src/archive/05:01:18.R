#Bayesian Quadrature with Gaussian Process

library(lhs)
library(data.tree)
library(matrixStats)
library(cubature)
library(tmvtnorm)
library(MASS)
library(plyr)
namedList<-treatSens:::namedList

#genz function for set of vectors (matrix)
f1 <- function(xx, u = rep(0.5, 1, ncol(xx)), a = rep(5, 1, ncol(xx))){
  
  y <- c()
  
  for (i in 1:nrow(xx)){
    
    y[i] <- realf1(xx[i, ])
    
  }
  
  return(y)
  
}

#genz function for a single vector
realf1 <- function(xx, u = rep(0.5, 1, length(xx)), a = rep(5, 1, length(xx))){
  
  sum <- 0
  
  for (i in 1:length(xx)){
    
    sum <- sum + a[i] * xx[i]
    
  }
  
  y <- (1 + sum) ^ (-(length(xx) + 1))
  
  return(10000*y)
  
}

#set parameters
N=10
dim = 3

#train data for X and Y
X <- randomLHS(N, dim);
Y <- f1(X)


p_t_mvnorm <- function(mean){
  
  return (ptmvnorm(lowerx = rep(0, dim), upperx = rep(1, dim), mean, sigma = diag(dim)))
  
}

#initialise outcomes
meanValue3 <- c()
sd3 <- c()

#Gaussian kernel
kernel <- function(x, y){
  
  #initialise kernel
  ker <- c()
  
  for (i in 1:nrow(y)){
    
    ker[i] <- exp(-norm(x - y[i, ], type = "2") ^ 2)
    
  }
  
  return(ker)
  
}

#covariance matrix with exponential Gaussian kernel
Gaussian_kernel <- function(x){
  
  K <- matrix(0, nrow = nrow(x), ncol = nrow(x))
  
  for (i in 1:nrow(x)){
    
    K[i, ] <- kernel(x[i, ], x)
    
  }
  
  return(K)
  
}
  
my_z <- function(x){
  
  z <- c()
  
  for(i in 1:nrow(x)){
    
    z[i] <- p_t_mvnorm(mean = x[i, ])[1]
    
  }
  
  return(z)
  
}

#initial covariance and z with training data
K <- Gaussian_kernel(X)
z <- my_z(X)

#loop
for (p in 1:100){
  
  meanValue3[p] <- t(z) %*% ginv(K) %*% Y
  
  sd3[p] <- t(z) %*% ginv(K) %*% z
  
  #randomly generate candidates
  candidate_set <- randomLHS(100, dim)
  
  #corresponding z of each candidate
  candidate_p <- mapply(p_t_mvnorm, alply(candidate_set, 1))
  
  #initialise variance of candidates
  candidate_var <- c()
  
  for (i in 1:100){
    
    #append a candidate to the training data
    candidate_data <- rbind(X, candidate_set[i, ])
    
    candidate_z <- c(z, candidate_p[[i]])
    
    #new covariance
    candidate_K <- Gaussian_kernel(candidate_data)
    
    #variance of the candidate
    candidate_var[i] <- t(candidate_z) %*% ginv(candidate_K) %*% candidate_z
    
  }
  
  #find the best candidate
  index <- which(candidate_var == min(candidate_var))[1]
  
  #append a new training data
  X<- rbind(X, candidate_set[index, ])
  Y <- c(Y, realf1(candidate_set[index, ]))
  
  #new selected kernel
  K <- Gaussian_kernel(X)
  
  #new z
  z <- my_z(X)

}

real <- adaptIntegrate(realf1, lowerLimit = rep(0, dim), upperLimit = rep(1, dim))

percentageError3 <- abs((meanValue3 - real[[1]]) / real[[1]]) * 100

plot(1:81, percentageError3)

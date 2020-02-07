h_func <- function(x) 
{
  return (1/(1 + exp(-80*x)))
}

gf_func <- function(x, f) 
{
  result <- rep(0, nrow(x))
  result[abs(x) >= 1] <- exp(-1/(1 - abs(x[abs(x)>=1])^2) + cos(f*pi*abs(x[abs(x)>=1])))
  return (result)
}

fisher_function <- function(x, C, R, H, F, P) {
  integrand <- 1
  for (i in 1:dim) {
    integrand <- integrand*(
      H[i]*gf_func((matrix(x[,i]) - C[i])/R[i], F[i]) + (-1)^P[i]*(
        0.5 - h_func(matrix(x[,i]) - C[i])
      )
    )
  }
  return (integrand)
}

dim <- 1
n <- 9
C <- replicate(dim, runif(9, 0.1, 0.9))
R <- replicate(dim, rbeta(9, 5, 2))
H <- replicate(dim, runif(9, 0.5*exp(1), 1.5*exp(1)))
F <- replicate(dim, runif(9, 0, 5))
P <- replicate(dim, rbinom(9, 1, 0.5))
trainX <- replicate(dim, runif(500, 0, 1))
par(mfrow = c(3, 3))
for (i in 1:9) {
  trainY <- fisher_function(trainX, C[i,1], R[i,1], H[i,1], F[i,1], P[i,1])
  plot(trainX[order(trainX),], trainY[order(trainX)], ty= "l")
} 




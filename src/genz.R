
#https://www.sfu.ca/~ssurjano/cont.html
#The gen z function in this script can be found in the above website.

#Continuous Genz function, take xx as a matrix of input (e.g. if we have 4 5-dim X),
#The input will be a 4*5 matrix, each row represents an input X
#It then returns the corresponding vector output of Y's
f1 <- function(xx, u = rep(0.5, 1, ncol(xx)), a = rep(5, 1, ncol(xx)))
{
  y <- c()
  for (i in 1:nrow(xx)){
    sumNum <- sum(5 * abs(xx[i, ] - 0.5))
    y[i] <- exp(-sumNum)
  }

  return(y)
}

#This function realf1 takes xx as a single input, or a vector. It returns a single Y
realf1 <- function(xx, u = rep(0.5, 1, length(xx)), a = rep(5, 1, length(xx)))
{
  sumNum <- sum(5 * abs(xx - 0.5))
  y <- exp(-sumNum)

  return(y)
}

#Previous function weighted by probability distribution standard normal
f<-function(x){
  
  return (f1(x)*dmvnorm(x));
}
  
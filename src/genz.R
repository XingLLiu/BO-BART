
#https://www.sfu.ca/~ssurjano/cont.html
#The gen z function in this script can be found in the above website.

#Continuous Genz function, take xx as a matrix of input (e.g. if we have 4 5-dim X),
#The input will be a 4*5 matrix, each row represents an input X
#It then returns the corresponding vector output of Y's
f1 <- function(xx, u = rep(0.5, 1, ncol(xx)), a = rep(5, 1, ncol(xx))){
  y<-c();
  for (i in 1:nrow(xx)){
    sum<-0
    for (j in 1:ncol(xx)){
      sum<-sum+a[j]*abs(xx[i,j]-u[j])
    }
    y[i]<-exp(-sum)
  }
  return(y)
}

#This function realf1 takes xx as a single input, or a vector.It returns a single Y
realf1 <- function(xx, u = rep(0.5, 1, length(xx)), a = rep(5, 1, length(xx))){
  sum<-0
  for (i in 1:length(xx)){
    sum<-sum+a[i]*abs(xx[i]-u[i])
  }
  y<-exp(-sum)
  return(y)
}

#Previous function weighted by probability distribution standard normal
f<-function(x){
  
  return (f1(x)*dmvnorm(x));
}
  
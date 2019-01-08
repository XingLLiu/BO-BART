
error<-list()
variance<-list()
prediction<-list()
realValue<-numeric()

#The function to be integrated, the case here is a discontinuous one
f1 <- function(xx, u = rep(0.5, 1, ncol(xx)), a = rep(5, 1, ncol(xx))){
  y<-c();
  for (i in 1:nrow(xx)){
    sum<-0
    for (j in 1:ncol(xx)){
      sum<-sum+a[j]*xx[j]
    }
    y[i]<-(1+sum)^(-(ncol(xx)+1))
  }
  return(y)
}

realf1 <- function(xx, u = rep(0.5, 1, length(xx)), a = rep(5, 1, length(xx))){
  sum<-0
  for (i in 1:length(xx)){
    sum<-sum+a[i]*xx[i]
  }
  y<-(1+sum)^(-(length(xx)+1))
  return(y)
}

#Previous function weighted by probability distribution standard normal
f<-function(x){
  
  return (f1(x)*dmvnorm(x));
  
}

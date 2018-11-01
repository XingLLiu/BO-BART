library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(mvtnorm)
library(cubature)
library(truncnorm)

dim=3
exp<-c()
sd3<-c()


distanceSquared<-function (x,c){
  dis<-(rowSums(t((t(x)-c)^2)));
  return (dis)
}

k<-function(x,c){
  y<-exp(-0.5*(1/1.28)*distanceSquared(x,c))
  return (y);
}


f1<-function(x){
  sum<-rep(0,dim(x)[1]);
  for (i in 1:10){
    sum<-sum+alpha[i]*k(x,c[i,])
  }
  return (sum)
}

f<-function(x){
  y<-f1(x)*dmvnorm(as.matrix(x));
  return (y);
}

c<-rmvnorm(n=10,mean=rep(0,dim));
alpha<-rnorm(10);

for (p in 1:40){
  

N=10*p;






  
X<-rmvnorm(N,mean=rep(0,dim));
Y<-f1(X)
model<-mlegp(X,Y)
Beta<-model$beta
sigma<-model$sig2
K<-matrix(0,nrow=N,ncol=N)



covFunction<-function(x,y){
  cov<-1
  for (i in 1:dim){
    cov<-cov*exp(-Beta[i]*sum((x[i]-y[i])^2))
  }
  return (cov)
}

for (i in 1:N){
  for (j in 1:N){
    cov<-sigma*covFunction(X[i,],X[j,])
    K[i,j]<-cov
  }
}





#real<-integrate(f,-Inf,Inf)[[1]]
b<-0;
B<-diag(dim);
a<-X;
A<-diag(1/(2*Beta),dim)
z<-c();
for(i in 1:N){
  z[i]<-sigma*(det(solve(A)%*%B+diag(dim))^-0.5)*exp(-0.5*(a[i,]+b)%*%solve(A+B)%*%(a[i,]-b))
}

exp[p]<-t(z)%*%ginv(K)%*%Y
sd3[p]<-sigma*det(2*solve(A)%*%B+diag(dim))^(-0.5)-t(z)%*%ginv(K)%*%z

}

realf1<-function(x){
  sum1<-0;
  for (i in 1:10){
    distance=sum((x-c[i,])^2);
    sum1<-sum1+alpha[i]*exp(-0.5*(1/1.28)*distance)
  }
  return (sum1)
}

realf<-function(x){
  return (realf1(x)*dmvnorm(x))
}
real<-adaptIntegrate(realf,lowerLimit = rep(-5,dim),upperLimit = rep(5,dim))[[1]]



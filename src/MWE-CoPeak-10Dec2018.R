library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(mvtnorm)
library(cubature)
library(truncnorm)
library(mlegp)
library(MASS)
namedList<-treatSens:::namedList

source("BART-BQ.R")
source("CoPeak-functions.R")
#BART Quadrature

#numerics to store standard deviation and mean value(approxiamtion)
meanValue<-numeric();
sd<-numeric()

#Initialize dataset, size of dataset depends on dimensions
set.seed(1)
trainX<-rmvnorm(99,mean=rep(0,dim))
trainY<-f1(trainX)
df<-data.frame(trainX,trainY)

#Run BART, keep finding approximation as new data point added in df
#store approxiamtion in mean value, posterior sd in standard deviation

for (i in 1:4){
  
  df<-Findmodel(df,1);
  print(sprintf("************: %d MAX trainY: %.04f", i, max(df$trainY)))
  model<-bart(df[1:dim],df$trainY,keeptrees =TRUE,keepevery=20L,nskip=1000,ndpost=2000,ntree=50,k=4,verbose=F)
  integrals<-sampleIntegrals(model,df[1:dim],0.95,2)
  
  #Scale back the approxiamtion and standard deviation
  ymin<-min(df$trainY); ymax<-max(df$trainY)
  print(sprintf("************* %d variance: %.04f", i, var(integrals)))
  sDeviation<-sqrt(var(integrals))*(ymax-ymin)
  scaledMean<-(mean(integrals)+0.5)*(ymax-ymin)+ymin 
  
  meanValue<-append(meanValue,scaledMean)
  sd<-append(sd,sDeviation)
  
}

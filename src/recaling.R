f1<-function(xprime){
  x=2*xprime+0.5
  y<-(sin(10*pi*x)/(2*x))+(x-1)^(4)
} 

trainX<-c(0,maximinLHS(8,1),1)
trainY<-f1(trainX)
rescale = function(x) (x-min(x))/(max(x)-min(x)) - .5
trainY = rescale(trainY)
range(trainY)
df<-data.frame(trainX,trainY)

model<-bart(df$trainX,df$trainY,keeptrees = TRUE,ntree=50,ndpost=500)
state<-model$fit$state
fits<-state[[1]]@savedTreeFits
posteriorDraw<-fits[1,1:50,1:500]
#this gives 500 posterior draws, 
#each consisting of 50 trees at the first point x=0
sumOverAllTrees<-colSums(posteriorDraw)
#sum the trees in each posterior draws, this returns E(y|x=0) in each posterior draw
m1<-mean(sumOverAllTrees)
#take mean in all posterior draw should give E(y|x=0) hopefully
prediction<-predict(model,0);
#this gives the prediction at x=0 in each posterior draw
m2<-mean(prediction)
#m1 and m2 should be the same...??
m1
m2

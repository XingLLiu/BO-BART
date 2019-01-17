library(dbarts)
library(data.tree)
library(hydroGOF)
source("/Users/christine798/Documents/GitHub/BO-BART/src/BART-BQ_reg.R")

netTrain <- read.csv("~/Documents/2017-2018/UROP/report/network/net_train.csv")
netTrain <- netTrain[-which(netTrain$V27 == 0), 2:29]
netTest <- read.csv("~/Documents/2017-2018/UROP/report/network/net_test.csv")
netTest <- netTest[-which(netTest$V27 == 0), 2:29]

n <- 50

highDimMSE <- c("benchmark", "BARTBQ")

for (dim in c(2)) {
   
   set.seed(1223)
   col <- c(sample(c(1:25, 27, 28), dim), 26)
   train <- netTrain[sample(nrow(netTrain), n, replace = FALSE), col]
   test <- netTest[sample(nrow(netTest), n/2, replace = FALSE), col]

   ymax <- max(train[, dim+1])
   ymin <- min(train[, dim+1])

   # trainX <- cbind(sample(1000, 100)/1000, sample(2000, 100)/2000)
   # trainY <- sin(trainX[, 1] + trainX[, 2])
   # trainFull <- cbind(trainX, trainY)

   # testX <- cbind(sample(1000, 10)/1000, sample(2000, 10)/2000)
   print(c("dim = ", dim))

   bart <- bart(train[, 1:dim], train[, dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=1000, ntree=50, k = 5)
   pred <- bartBQPredict(bart, train, test[, 1:dim], ymax, ymin)

   #bart <- bart(trainFull[, 1:dim], trainFull[, dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=1000, ntree=50, k = 2)

   bartBench <- bart(train[, 1:dim], train[, dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=1000, ntree=50, k = 5)
   benchPred <- colMeans(predict(bartBench, test[, 1:dim]))
   #pred <- bartBQPredict(bart, trainFull, testX[, 1:dim], max(trainFull[, dim+1]), min(trainFull[, dim+1]))
   
   benchMSE <- mse(test[, dim+1], benchPred)
   MSE <- mse(test[, dim+1], pred)

   highDimMSE <- rbind(highDimMSE, c(benchMSE, MSE))

}




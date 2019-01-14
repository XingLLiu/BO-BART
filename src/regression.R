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

for (dim in c(2, 3, 5, 10, 20)) {
   
   set.seed(1223)
   col <- c(sample(c(1:25, 27, 28), dim), 26)
   train <- netTrain[sample(nrow(netTrain), n, replace = FALSE), col]
   test <- netTest[sample(nrow(netTest), n/2, replace = FALSE), col]

   bartBench <- bart(train[, 1:dim], train[, dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=1000, ntree=50, k = 5)
   bart <- bart(rescale(train[, 1:dim]), train[, dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=1000, ntree=50, k = 5)
   
   print(c("dim = ", dim))

   benchPred <- colMeans(predict(bartBench, test[, 1:dim]))
   pred <- bartBQPredict(bart, rescale(test[, 1:dim]))

   benchMSE <- mse(test[, dim+1], benchPred)
   MSE <- mse(test[, dim+1], pred)

   highDimMSE <- rbind(highDimMSE, c(benchMSE, MSE))

}




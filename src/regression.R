library(dbarts)
library(data.tree)
library(hydroGOF)
source("/Users/christine798/Documents/GitHub/BO-BART/src/BART-BQ_reg.R")

netTrain <- read.csv("~/Documents/2017-2018/UROP/report/network/net_train.csv")
netTrain <- netTrain[, 2:29]
netTest <- read.csv("~/Documents/2017-2018/UROP/report/network/net_test.csv")
netTest <- netTest[, 2:29]

n <- 50
dim <- 3

set.seed(1223)
col <- c(sample(28, dim), 26)
train <- netTrain[sample(nrow(netTrain), n, replace = FALSE), col]
test <- netTest[sample(nrow(netTest), n/2, replace = FALSE),  col]

bartBench <- bart(train[, 1:3], train[,dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=1000, ntree=50, k = 5)
bart <- bart(rescale(train[, 1:3]), train[,dim+1], keeptrees=TRUE, keepevery=20L, nskip=1000, ndpost=1000, ntree=50, k = 5)


benchPred <- colMeans(predict(bartBench, test))
pred <- bartBQPredict(bart, test)

benchMSE <- mse(test[, 4], benchPred)
MSE <- mse(test[, 4], pred)

print(c(benchMSE, MSE))
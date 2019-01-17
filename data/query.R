<<<<<<< HEAD
trainData <- read.csv("../data/train.csv")
candidateData <- read.csv("../data/candidates.csv")
populationData <- read.csv("../data/fullData.csv")
=======
trainData <- read.csv("../../data/train.csv")
candidateData <- read.csv("../../data/candidates.csv")
populationData <- read.csv("../../data/fullData.csv")
>>>>>>> 1b9934dba4506a121b703a4208cf1edd9eaa7cc6

# extract covariates and response
cols <- dim(trainData)[2] - 1
trainX <- trainData[,-c(1, (cols+1))]
trainY <- trainData[, (cols+1)]

candidateX <- candidateData[,-c(1, (cols+1))]
candidateY <- candidateData[, cols+1]

# Compute population mean
<<<<<<< HEAD
source("/Users/christine798/Documents/GitHub/Bo-BART/src/regression/bartMean.R")
=======
source("./bartMean.R")
>>>>>>> 1b9934dba4506a121b703a4208cf1edd9eaa7cc6
populationMean <- mean(populationData[,ncol(populationData)]) # 38115.99
BARTResults_population <- computePopulationMean(trainX, trainY, candidateX, candidateY, num_iterations = 5)
#random sampling
MImean <- c()
for (i in 1:1000) {
  MImean[i] <- mean(c(trainY, candidateY[1:i]))
}

plot(MImean, type = "l")
plot(BARTResults$meanValueBART, type = "l")

cat("BART Population Mean", BARTResults$meanValueBART)
cat("True Population Mean", populationMean)

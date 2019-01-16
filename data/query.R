trainData <- read.csv("../data/train.csv")
candidateData <- read.csv("../data/candidates.csv")
populationData <- read.csv("../data/fullData.csv")

# extract covariates and response
cols <- dim(trainData)[2] - 1
trainX <- trainData[,-c(1, (cols+1))]
trainY <- trainData[, (cols+1)]

candidateX <- candidateData[,-c(1, (cols+1))]
candidateY <- candidateData[, cols+1]

dim <- dim(trainX)
# Compute population mean
source("/Users/christine798/Documents/GitHub/Bo-BART/src/regression/bartMean.R")
populationMean <- mean(populationData[,ncol(populationData)]) # 38115.99
BARTResults <- computePopulationMean(dim, trainX, trainY, candidateX, candidateY, num_iterations = 100)
#random sampling
MImean <- c()
for (i in 1:1000) {
  MImean[i] <- mean(c(trainY, candidateY[1:i]))
}

plot(MImean, type = "l")
plot(BARTResults$meanValueBART, type = "l")

cat("BART Population Mean", BARTResults$meanValueBART)
cat("True Population Mean", populationMean)

# male mean
trainX <- trainX[trainX$SEX == 1,]
trainX <- trainY[trainX$SEX == 1]


# MALE MEAN
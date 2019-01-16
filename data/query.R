library(sqldf)
library(readr)

trainData <- read.csv("../../data/train.csv")
candidateData <- read.csv("../../data/candidates.csv")
populationData <- read.csv("../../data/fullData.csv")


# extract covariates and response
cols <- dim(trainData)[2] - 1
trainX <- trainData[,-c(1, (cols+1))]
trainY <- trainData[, (cols+1)]

candidateX <- candidateData[,-c(1, (cols+1))]
candidateY <- candidateData[, cols+1]

# Compute population mean
source("bartMean.R")
populationMean <- mean(populationData[,dim+1]) # 38115.99
BARTResults <- computePopulationMean(dim, trainX, trainY, candidateX, candidateY, num_iterations = 5)

BARTResults$meanValueBART



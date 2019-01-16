library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)

# read in data
trainData <- read.csv("../../data/train.csv")
testData <- read.csv("../../data/validation.csv")
candidateData <- read.csv("../../data/candidates.csv")
populationData <- rbind(trainData, testData, candidateData)

dim <- dim(trainData)[2] - 1

# extract covariates and response
trainX <- trainData[,-(dim+1)]
trainY <- trainData[, dim+1]

candidateX <- candidateData[,-(dim+1)]
candidateY <- candidateData[, dim+1]


# Compute population mean
populationMean <- mean(populationData[,dim+1]) # 38115.99
source("bartMean.R")


# Racial mean income


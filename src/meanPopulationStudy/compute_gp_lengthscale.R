library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(caret)
set.seed(0)
source("src/meanPopulationStudy/bartMean2.R")
source("src/meanPopulationStudy/gpMean_2.R")
source("src/optimise_gp.R")

# paths to save results and plots
resultPath <- "results/populationStudy/"
plotPath <- "Figures/populationStudy/"

args <- as.double(commandArgs(TRUE))
num_cv_start <- args[1]
num_cv_end <- args[2]
num_data <- args[3]  # set to 2000 for this lengthscale
num_design <- args[4]

# read in data
trainData <- read.csv("data/train2.csv")
candidateData <- read.csv("data/candidate2.csv")

# convert num to factor, log income
convert <- function(data) {
  
  log_Total_person_income <- log(data[, ncol(data)])
  data <- sapply(data[, 2:(ncol(data)-1)], as.factor)
  data <- data.frame(data)
  data <- cbind(data, log_Total_person_income)
}

trainData <- convert(trainData)
candidateData <- convert(candidateData)

trainData <- trainData[!is.infinite(trainData$log_Total_person_income),]
candidateData <- candidateData[!is.infinite(candidateData$log_Total_person_income),]
trainData <- trainData[complete.cases(trainData),]
candidateData <- candidateData[complete.cases(candidateData),]
# compute the real population mean log income
# poptMean <- mean(c(trainData$log_Total_person_income, candidateData$log_Total_person_income))
lengthscales <- c()
for (num_cv in num_cv_start:num_cv_end) {
  # set new seed
  set.seed(num_cv)
  print(num_cv)
  trainData <- trainData[sample(c(1:dim(trainData)[1]), num_design),]
  candidateData <- candidateData[1:num_data, ]
  # poptMean <- mean(c(trainData$log_Total_person_income, candidateData$log_Total_person_income))

  # extract covariates and response
  cols <- ncol(trainData)
  trainX <- trainData[, -cols]
  trainY <- trainData[, cols]
  candidateX <- candidateData[, -cols]
  candidateY <- candidateData[, cols]

  # one-hot encoding
  trainX.num <- trainX
  candidateX.num <- candidateX
  dummyFullData <- dummyVars("~.", data = rbind(trainX, candidateX))
  trainX <- data.frame(predict(dummyFullData, newdata = trainX))
  candidateX <- data.frame(predict(dummyFullData, newdata = candidateX))
  lengthscale <- optimise_gp_r(as.matrix(trainX), trainY, kernel = "rbf", epochs = 500)
  lengthscales[num_cv] <- lengthscale
}
lengthscales <- data.frame(lengthscales)
lengthscales$num_cv <- c(num_cv_start:num_cv_end)

write.csv(lengthscales, file = paste0(plotPath, "lengthscales_", num_design, "_",  num_data, ".csv"), row.names=FALSE)
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
num_new_surveys <- args[1]
num_cv_start <- args[2]
num_cv_end <- args[3]
num_data <- args[4]   # size of candidate set
num_design <- args[5]   # size of initial training set
jitter <- 1e-6
cutoff <- 10


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
# encode as continuous
trainData$Education <- as.double(as.character(trainData$Education))
candidateData$Education <- as.double(as.character(candidateData$Education))
# remove missing values
trainData_full <- trainData[complete.cases(trainData),]
candidateData <- candidateData[complete.cases(candidateData),]
# keep the first num_data observations as candidate set
candidateData <- candidateData[1:num_data,]
# change response to binary
binary_response <- ifelse(trainData_full$log_Total_person_income >= cutoff, 1 - jitter, jitter)
trainData_full$log_Total_person_income <- binary_response
binary_response <- ifelse(candidateData$log_Total_person_income >= cutoff, 1 - jitter, jitter)
candidateData$log_Total_person_income <- binary_response

# compute the real population mean log income
# lengthscales <- read.csv("Figures/populationStudy/lengthscales_hicov_10000.csv")
# ground_truths <- read.csv(paste("results/populationStudy/popt_hicov", num_design,"_", num_data, ".csv", sep=""))

for (num_cv in num_cv_start:num_cv_end) {
  # set new seed
  set.seed(num_cv)
  print(num_cv)
  trainData <- trainData_full[sample(c(1:dim(trainData_full)[1]), num_design),]
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
  # save(trainX, trainX.num, trainY, candidateX, candidateX.num, candidateY, file = "data/survey_data.RData")
  # load(file = "data/survey_data.RData")
  # compute population average income estimates by BARTBQ
  t0 <- proc.time()
  BARTresults <- computeBART(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys, save_posterior=TRUE, num_cv=num_cv)
  t1 <- proc.time()
  bartTime <- (t1 - t0)[[1]]
  # population average income estimation by Monte Carlo
  MIresults <- computeMI(trainX.num, trainY, candidateX.num, candidateY, num_iterations=num_new_surveys, seed = num_cv)
  
  # GPBQ
  lengthscale <- optimise_gp_r(as.matrix(trainX), trainY, kernel = "matern32", epochs = 500)
  t0 <- proc.time()
  GPresults <- computeGPBQEmpirical(as.matrix(trainX), trainY, as.matrix(candidateX), candidateY, epochs=num_new_surveys, kernel="matern32", lengthscale=lengthscale)
  t1 <- proc.time()
  GPTime <- (t1 - t0)[[1]]
  
  # store results
  results <- data.frame(
    "epochs" = c(1:num_new_surveys),
    "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
    "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
    "GPMean" = GPresults$meanValueGP, "GPsd" = sqrt(GPresults$varianceGP),
    # "PoptMean" = ground_truths$mi_ground_truths[1], "BpoptMean" = ground_truths$bart_ground_truths[1],
    "runtimeBART" = rep(bartTime, num_new_surveys),
    "runtimeGP" = rep(GPTime, num_new_surveys)
  )
  write.csv(results, file = paste0(resultPath, "results", num_cv, ".csv"), row.names=FALSE)
  results_models <- list("BART"=BARTresults, "MI"=MIresults, "GP"=GPresults)
  save(results_models, file = paste0(plotPath, "results", num_cv, ".RData"))
  
  real <- results$PoptMean[1]
  Breal <- results$BpoptMean[1]
  
  print(c("Real BART-Int: ", Breal))
  print(c("Real MC: ", real))
  print(c("BART-Int: ", results$BARTMean[num_new_surveys]))
  print(c("GP-BQ: ", results$GPMean[num_new_surveys]))
  print(c("MI: ", results$MIMean[num_new_surveys]))
}

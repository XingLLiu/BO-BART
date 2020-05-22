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
num_data <- args[4]   # set to 2000 for this lengthscale
num_design <- args[5]   # set to 2000 for this lengthscale


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
trainData_full <- trainData[complete.cases(trainData),]
candidateData <- candidateData[complete.cases(candidateData),]
candidateData <- candidateData[1:num_data, ]
# compute the real population mean log income
poptMean <- mean(c(trainData$log_Total_person_income, candidateData$log_Total_person_income))
lengthscales <- read.csv("Figures/populationStudy/lengthscales_10000.csv")
ground_truths <- read.csv(paste("results/populationStudy/popt_", num_design,"_", num_data, ".csv", sep=""))

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
    BARTresults <- computeBART(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys)
    t1 <- proc.time()
    bartTime <- (t1 - t0)[[1]]
    # population average income estimation by Monte Carlo
    # MIresults <- computeMI(trainX.num, trainY, candidateX.num, candidateY, num_iterations=num_new_surveys)
    MIresults <- computeMI(trainX.num, trainY, candidateX.num, candidateY, num_iterations=num_new_surveys, seed = num_cv)
    # plot(MIresults$meanValueMI, ylim = c(10.9, 11.1), xlab = "num_iterations", ylab = "mean population")
    # legend("topright", legend=c("MC integration"))
    # abline(h = poptMean, col = "red")
    
    # GPBQ
    lengthscale <- lengthscales$lengthscales[num_cv]

    t0 <- proc.time()
    GPresults <- computeGPBQEmpirical(as.matrix(trainX), trainY, as.matrix(candidateX), candidateY, epochs=num_new_surveys, lengthscale=lengthscale)
    t1 <- proc.time()
    GPTime <- (t1 - t0)[[1]]
    
    # store results
    results <- data.frame(
         "epochs" = c(1:num_new_surveys),
         "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
         "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
         "GPMean" = GPresults$meanValueGP, "GPsd" = GPresults$varianceGP,
         "PoptMean" = ground_truths$mi_ground_truths[1], "BpoptMean" = ground_truths$bart_ground_truths[1],
         "runtimeBART" = rep(bartTime, num_new_surveys),
         "runtimeGP" = rep(GPTime, num_new_surveys)
     )
    # results <- data.frame(
    #     "epochs" = c(1:num_new_surveys),
    #     "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
    #     "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
    #     "BRSMean" = BRSresults$meanValueBRS, "BRSsd" = BRSresults$standardDeviationBRS,
    #     "GPMean" = GPresults$meanValueGP, "GPsd" = GPresults$varianceGP,
    #     "PoptMean" = poptMean
    # )
	  #results <- data.frame("epochs"=c(1:num_new_surveys), "GPMean"=GPresults$meanValueGP, "GPsd"=GPresults$varianceGP)
    write.csv(results, file = paste0(resultPath, "results", num_cv, ".csv"), row.names=FALSE)
    results_models <- list("BART"=BARTresults, "MI"=MIresults, "GP"=GPresults)
    save(results_models, file = paste0(plotPath, "results", num_cv, ".RData"))
    
    #results_models <- list("GP"=GPresults)
  	#save(results_models, file = paste0(plotPath, "gpresults", num_cv, ".RData"))
    
    real <- results$PoptMean[1]
    Breal <- results$BpoptMean[1]
    
    print(c("Real BART-Int: ", Breal))
    print(c("Real MC: ", real))
    print(c("BART-Int: ", results$BARTMean[num_new_surveys]))
    print(c("GP-BQ: ", results$GPMean[num_new_surveys]))
    print(c("MI: ", results$MIMean[num_new_surveys]))
}


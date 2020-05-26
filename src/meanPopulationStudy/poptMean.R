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
num_design <- args[5]

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
poptMean <- mean(c(trainData$log_Total_person_income, candidateData$log_Total_person_income))

for (num_cv in num_cv_start:num_cv_end) {
  
    # set new seed
    set.seed(num_cv)
    print(num_cv)
    trainData <- trainData[sample(c(1:dim(trainData)[1]), num_design),]
    candidateData <- candidateData[1:num_data, ]

    # # linear regression toy example
    # fullData <- rbind(trainData, candidateData)
    # lm.fit <- lm(log_Total_person_income~., data = fullData)
    # res.sd <- sd(lm.fit$residuals)
    # trainData$log_Total_person_income <- lm.fit$fitted.values[1:nrow(trainData)] + rnorm(nrow(trainData), res.sd)
    # candidateData$log_Total_person_income <- lm.fit$fitted.values[-(1:nrow(trainData))] + rnorm(nrow(candidateData), res.sd)
    # poptMean <- mean(c(trainData$log_Total_person_income, candidateData$log_Total_person_income))

    # # reconstruct order of the data
    # set.seed(2020)
    # trainData <- trainData[sample(nrow(trainData)), ]
    # candidateData <- candidateData[sample(nrow(candidateData)), ]

    # extract covariates and response
    cols <- ncol(trainData)
    # trainX <- trainData[1:500, -cols]
    # trainY <- trainData[1:500, cols]


    # candidateX <- candidateData[1:5000, -cols]
    # candidateY <- candidateData[1:5000, cols]

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
    MIresults <- computeMI(trainX.num, trainY, candidateX.num, candidateY, num_iterations=nrow(candidateX.num), seed = num_cv)
    # plot(MIresults$meanValueMI, ylim = c(10.9, 11.1), xlab = "num_iterations", ylab = "mean population")
    # legend("topright", legend=c("MC integration"))
    # abline(h = poptMean, col = "red")
    
    # GPBQ
    lengthscale <- optimise_gp_r(as.matrix(trainX), trainY, kernel = "rbf", epochs = 500)
    t0 <- proc.time()
    GPresults <- computeGPBQEmpirical(as.matrix(trainX), trainY, as.matrix(candidateX), candidateY, epochs=num_new_surveys, lengthscale=lengthscale)
    t1 <- proc.time()
    GPTime <- (t1 - t0)[[1]]
    
    # store results
    results <- data.frame(
         "epochs" = c(1:num_new_surveys),
         "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
         "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
         "GPMean" = GPresults$meanValueGP, "GPsd" = sqrt(GPresults$varianceGP), "PoptMean" = poptMean,
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
    
    real <- results$PoptMean[1]
    
    print(c("Real MC: ", real))
    print(c("BART-Int: ", results$BARTMean[num_new_surveys]))
    print(c("GP-BQ: ", results$GPMean[num_new_surveys]))
    print(c("MI: ", results$MIMean[num_new_surveys]))
}

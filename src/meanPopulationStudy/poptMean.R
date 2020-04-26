library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(caret)
set.seed(0)
source("src/meanPopulationStudy/bartMean.R")
source("src/meanPopulationStudy/gpMean_2.R")
source("src/optimise_gp.R")

# paths to save results and plots
resultPath <- "results/populationStudy/"
plotPath <- "Figures/populationStudy/"

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args

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

# reconstruct order of the data
set.seed(2020)
trainData <- trainData[sample(nrow(trainData)), ]
candidateData <- candidateData[sample(nrow(candidateData)), ]

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

for (num_cv in 1:5) {
    # set new seed
    set.seed(num_cv)

    # compute population average income estimates by BARTBQ
    BARTresults <- computeBART(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys)
    
    # population average income estimation by Monte Carlo
    MIresults <- computeMI(trainX.num, trainY, candidateX.num, candidateY, num_iterations=num_new_surveys)
    
    # GPBQ
    # lengthscale <- optimise_gp_r(as.matrix(one_hot(data.table(trainX))), trainY, kernel = "rbf", epochs = 50000)
    # GPresults <- computeGPBQEmpirical(as.matrix(trainX), trainY, as.matrix(candidateX), candidateY, epochs=num_new_surveys, lengthscale=lengthscale)
    
    # population average income estimation by block random sampling
    BRSresults <- computeBRS(trainX.num, trainY, candidateX.num, candidateY, group = "Race", num_iterations=num_new_surveys)
    
    # store results
    # results <- data.frame(
    #     "epochs" = c(1:num_new_surveys),
    #     "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
    #     "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
    #     "BRSMean" = BRSresults$meanValueBRS, "BRSsd" = BRSresults$standardDeviationBRS, 
    #     "PoptMean" = poptMean, "BpoptMean" = BARTpoptMean
    # )
    results <- data.frame(
        "epochs" = c(1:num_new_surveys),
        "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
        "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
        "BRSMean" = BRSresults$meanValueBRS, "BRSsd" = BRSresults$standardDeviationBRS,
        "PoptMean" = poptMean
    )
    write.csv(results, file = paste0(resultPath, "results", num_cv, ".csv"), row.names=FALSE)
    results_models <- list("BART"=BARTresults, "MI"=MIresults, "BRS"=BRSresults)
    save(results_models, file = paste0(plotPath, "results", num_cv, ".RData"))
    
    real <- results$PoptMean[1]
    # Breal <- results$BpoptMean[1]

    # 1. Open jpeg file
    pdf(paste0(plotPath, "results", num_cv, ".pdf"), width = 8, height = 10)
    par(mfrow = c(1,2), pty = "s")
    # ymax <- max(c(abs(results$BARTMean - real), abs(results$BRSMean - real), abs(results$MIMean - real)))
    ymax <- max(c(abs(results$BARTMean - real), abs(results$GPMean - real)))
    plot(results$epochs,
         abs(results$BARTMean - results$PoptMean),
         ty="l",
         ylab = "Absolute Error",
         xlab = "num_iterations",
         col = "orangered",
         ylim = c(0, ymax)
    )
    abline(h=0)
    
    points(results$epochs, abs(results$MIMean - real), ty="l", col = "chartreuse4")
    points(results$epochs, abs(results$BRSMean - real), ty="l", col = "dodgerblue")

    legend("topright", legend=c("Block sampling", "BART-Int", "Monte Carlo sampling"),
           col=c("dodgerblue", "orangered", "chartreuse4"), cex=0.8, lty = c(1,1,1,1))
    
    # ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$BRSMean - 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
    # ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$BRSMean + 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
    ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$GPMean - 2*results$GPsd))
    ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$GPMean + 2*results$GPsd))    
    
    plot(results$epochs, 
         results$MIMean, 
         ty="l", 
         ylab = "Mean Population",
         xlab = "num_iterations",
         col = "chartreuse4",
         ylim = c(ymin, ymax)
    )
    polygon(c(results$epochs, rev(results$epochs)), 
            c(
              results$MIMean + 2*results$MIsd, 
              rev(results$MIMean - 2*results$MIsd)
            ), 
            col = adjustcolor("chartreuse4", alpha.f = 0.10), 
            border = "chartreuse4", lty = c("dashed", "solid"))
    points(results$epochs, results$BRSMean, ty="l", col = "dodgerblue")
    polygon(c(results$epochs, rev(results$epochs)), 
            c(
              results$BRSMean + 2*results$BRSsd, 
              rev(results$BRSMean - 2*results$BRSsd)
            ), 
            col = adjustcolor("dodgerblue", alpha.f = 0.10), 
            border = "dodgerblue", lty = c("dashed", "solid"))
    points(results$epochs, results$BARTMean, ty="l", col = "orangered")
    polygon(c(results$epochs, rev(results$epochs)), 
            c(
              results$BARTMean + 2*results$BARTsd, 
              rev(results$BARTMean - 2*results$BARTsd)
            ), 
            col = adjustcolor("orangered", alpha.f = 0.10), 
            border = "orangered", lty = c("dashed", "solid"))
    abline(h=results$PoptMean)
    dev.off()
}

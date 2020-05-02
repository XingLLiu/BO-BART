setwd(getwd())

library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(data.table)
library(mltools)
# set seed to enable reproduction

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args


# read in data
trainData <- read.csv("train2.csv")
candidateData <- read.csv("candidate2.csv")
populationData <- read.csv("../../data/full_data.csv")

cols <- ncol(trainData)

# convert num to factor, log income
convert <- function(data) {
  
  log_Total_person_income <- log(data[, cols])
  data <- sapply(data[, 2:(cols-1)], as.factor)
  data <- data.frame(data)
  data <- cbind(data, log_Total_person_income)
}

trainData <- convert(trainData)
candidateData <- convert(candidateData)
populationData <- convert(populationData)


# reconstruct order of the data
trainData <- trainData[sample(nrow(trainData)), ]
candidateData <- candidateData[sample(nrow(candidateData)), ]

# extract covariates and response
trainX <- trainData[, -(cols-1)]
trainY <- trainData[, cols-1]

candidateX <- candidateData[, -(cols-1)]
candidateY <- candidateData[, cols-1]

# compute the real population mean income
poptMean <- mean(populationData$log_Total_person_income)

for (num_cv in 1:5) {

    # compute population average income estimates by BARTBQ
    source("./bartMean.R")
    BARTresults <- computeBART(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys)
    
    # compute population average income estimates by GPBQ
    source("./gpMean.R")
    GPresults <- computeGPBQEmpirical(as.matrix(one_hot(data.table(trainX))), trainY, 38, as.matrix(one_hot(data.table(candidateX))), candidateY, num_new_surveys, kernel="rbf", lengthscale=4.805, sequential=TRUE) 
    
    # store results
    results <- data.frame(
        "epochs" = c(1:num_new_surveys),
        "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
        "GPMean" = GPresults$meanValueGP, "GPsd" = sqrt(GPresults$varianceGP),
        "PoptMean" = poptMean
    )
    write.csv(results, file = paste("catResults", num_cv, ".csv", sep=""), row.names=FALSE)
    
    real <- results$PoptMean[1]

    # 1. Open jpeg file
    pdf(paste("catResults", num_cv, ".pdf", sep=""), width = 8,5, height = 10)
    par(mfrow = c(1,2), pty = "s")
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
    
    points(results$epochs, abs(results$GPMean - real), ty="l", col = "chartreuse4")

    legend("topright", legend=c("BART-Int", "GP-int"),
           col=c("orangered", "chartreuse4"), cex=0.8, lty = c(1,1,1,1))
    
    ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$GPMean - 2*results$GPsd))
    ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$GPMean + 2*results$GPsd))
    
    plot(results$epochs, 
         results$GPMean, 
         ty="l", 
         ylab = "Mean Population",
         xlab = "num_iterations",
         col = "chartreuse4",
         ylim = c(ymin, ymax)
    )
    polygon(c(results$epochs, rev(results$epochs)), 
            c(
              results$GPMean + 2*results$GPsd, 
              rev(results$GPMean - 2*results$GPsd)
            ), 
            col = adjustcolor("chartreuse4", alpha.f = 0.10), 
            border = "chartreuse4", lty = c("dashed", "solid"))
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

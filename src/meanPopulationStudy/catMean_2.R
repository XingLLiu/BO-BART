setwd(getwd())

library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(data.table)
library(mltools)
library(reticulate)
# set seed to enable reproduction

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args


# read in data
populationData <- read.csv("../../data/full_data.csv")

cols <- ncol(populationData)
index <- populationData[, 1]
Log_Total_person_income <- log(populationData[, cols])
populationData <- sapply(populationData[, 2:(cols-1)], as.factor)
populationData <- cbind(index, populationData, Log_Total_person_income)
populationData <- data.frame(populationData)
populationData[, cols] <- Log_Total_person_income

train <- c()
ref <- c()
candidate <- populationData

set.seed(123)


for (i in 2:(cols-1)) {

  level <- levels(populationData[, i])
  nlevel <- length(level)

  for (j in 1:nlevel){

    ss <- candidate[which(candidate[, i] == level[j]), ]
    n <- sample(nrow(ss), 1)
    train <- rbind(train, ss[n, ])
    candidate <- candidate[-which(candidate[, 1] == ss[n, 1]), ]
    }
}


trainData <- train[, -1]
candidateData <- candidate[, -1]



# extract covariates and response
trainX <- trainData[, -(cols-1)]
trainY <- trainData[, cols-1]

candidateX <- candidateData[, -(cols-1)]
candidateY <- candidateData[, cols-1]

# compute the real population mean income
poptMean <- mean(log(populationData$Log_Total_person_income))

for (num_cv in 1:1) {

    # compute population average income estimates by BARTBQ
    source("./bartMean.R")
    BARTresults <- computeBART(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys)
    
    # compute population average income estimates by GPBQ
    source("./gpMean.R")
    source("../optimise_gp.R")
    ls <- optimise_gp_r(as.matrix(one_hot(data.table(trainX))), trainY, kernel = "rbf", epochs = 50000)
    GPresults <- computeGPBQEmpirical(as.matrix(one_hot(data.table(trainX))), trainY, 38, as.matrix(one_hot(data.table(candidateX))), candidateY, num_new_surveys, kernel="rbf", lengthscale=ls, sequential=TRUE) 
    
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
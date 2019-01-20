setwd("~/Documents/GitHub/BO-BART/src/meanPopulationStudy/")

library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args[1]

# read in data
trainData <- read.csv("../../data/collegeTrain.csv")
candidateData <- read.csv("../../data/collegeRemain.csv")
populationData <- read.csv("../../data/collegeData.csv")

#extract covariates and response
cols <- dim(trainData)[2] - 1
trainX <- trainData[,-c(1, (cols+1))]
trainY <- trainData[, (cols+1)]

candidateX <- candidateData[,-c(1, (cols+1))]
candidateY <- candidateData[, (cols+1)]

#trainX <- trainX[1:50, ]
#trainY <- trainY[1:50]

# compute the real population mean income
poptMean <- mean(populationData$Grad.Rate)

# Compute population mean
source("./bartMean.R")
BARTResults <- computePopulationMean(trainX, trainY, candidateX, candidateY, num_iterations = num_new_surveys)

MImean <- c()
MIstandardDeviation <- c()
for (i in 1:num_new_surveys) {
MImean[i] <- mean(c(trainY, candidateY[1:i]))
MIstandardDeviation[i] <- sqrt( var(c(trainY, candidateY[1:i])) ) 
}

BR <- BRcomputeMean(trainX, trainY, candidateX, candidateY, num_iterations = num_new_surveys)

results <- data.frame(
    "epochs" = c(1:num_new_surveys),
    "BARTMean" = BARTResults$meanValueBART, "BARTsd" = BARTResults$standardDeviationBART,
    "MIMean" = MImean, "MIsd" = MIstandardDeviation, "BRmean" = BR$BRmean, "BRsd" = BR$BRstandardDeviation,
    "PoptMean" = poptMean
)
write.csv(results, file = "./population.csv",row.names=FALSE)

cat("BART Population Mean", BARTResults$meanValueBART)
cat("True Population Mean", poptMean)

# 1. Open jpeg file
png("./population.png", width = 700, height = 583)
# 2. Create the plot
par(mfrow = c(1,2), pty = "s")
plot(x = c(1:num_new_surveys), y = MImean,
    pch = 16, type = "l",
    xlab = "Number of Queries", ylab = "Population mean", col = "blue",
    main = NULL,
    lty = 1,
    ylim = c(60,75),
    xaxs="i", yaxs="i"
    )
lines(x = c(1:num_new_surveys), BARTResults$meanValueBART, type = 'l', col = "red", lty = 1)
lines(x = c(1:num_new_surveys), BR$BRmean, type = 'l', col = "green", lty = 1)
abline(a = poptMean, b = 0, lty = 1, col = "black")
legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling", "Actual Mean"),
        col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1))

plot(x = c(1:num_new_surveys), y = MIstandardDeviation,
    pch = 16, type = "l",
    xlab = "Number of Queries", ylab = "Standard deviation", col = "blue",
    main = NULL,
    lty = 1,
    ylim = c(0.8, 1.5),
    xaxs="i", yaxs="i"
    )
lines(x = c(1:num_new_surveys), BARTResults$standardDeviationBART, type = 'l', col = "red", lty = 1)
lines(x = c(1:num_new_surveys), BR$BRstandardDeviation, type = 'l', col = "green", lty = 1)
legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling"),
        col=c("blue", "red", "green"), cex=0.8, lty = c(1,1))

# 3. Close the file
dev.off()


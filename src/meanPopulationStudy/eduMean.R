setwd("~/Documents/GitHub/BO-BART/src/meanPopulationStudy/")

library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)

# read in data
trainData <- read.csv("../../data/train.csv")
candidateData <- read.csv("../../data/candidate.csv")
populationData <- read.csv("../../data/full_data.csv")

studies <- c("Highschool and below", "Beyond highschool")
index <- matrix(c(1,17,16,24), nrow = 2, ncol = 2)

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args[1]

for (edu in 1:2) {

  study <- studies[edu]
  # extract covariates and response
  cols <- dim(trainData)[2] - 1
  trainX <- trainData[,-c(1, (cols+1))]
  trainY <- log(trainData[, (cols+1)])

  candidateX <- candidateData[,-c(1, (cols+1))]
  candidateY <- log(candidateData[, (cols+1)])

  # stratify the male population
  trainY <- trainY[trainX$Education >= index[edu, 1] & trainX$Education <= index[edu, 2]]
  trainX <- trainX[(trainX$Education >= index[edu, 1] & trainX$Education <= index[edu, 2]), -4]

  candidateY <- candidateY[candidateX$Education >= index[edu, 1] & candidateX$Education <= index[edu, 2]]
  candidateX <- candidateX[(candidateX$Education >= index[edu, 1] & candidateX$Education <= index[edu, 2]), -4]

  trainX <- trainX[1:100, ]
  trainY <- trainY[1:100]

  # compute the real population mean income
  poptMean <- mean(log(populationData[(populationData$Education >= index[edu, 1] & populationData$Education <= index[edu, 2]), ]$Total_person_income))

  # run the BART regression model
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
  write.csv(results, file = "./%s.csv" %--% c(study), row.names=FALSE)

  cat("BART %s Mean" %--% c(study), BARTResults$meanValueBART[num_new_surveys])
  cat("True Population Mean", poptMean)

  # 1. Open jpeg file
  png("./%s.png" %--% c(study), width = 700, height = 583)
  # 2. Create the plot
  par(mfrow = c(1,2), pty = "s")
  plot(x = c(1:num_new_surveys), y = MImean,
       pch = 16, type = "l",
       xlab = "Number of Queries", ylab = "Population mean", col = "blue",
       main = NULL,
       lty = 1,
       ylim = c(9.5+(edu-1)/2, 10.5+(edu-1)/2),
       xaxs="i", yaxs="i"
       )
  lines(x = c(1:num_new_surveys), BARTResults$meanValueBART, type = 'l', col = "red", lty = 1)
  lines(x = c(1:num_new_surveys), BR$BRmean, type = 'l', col = "green", lty = 1)
  abline(a = poptMean, b = 0, lty = 1)
  legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling",  "Actual Mean"),
         col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1))
  
  plot(x = c(1:num_new_surveys), y = MIstandardDeviation,
       pch = 16, type = "l",
       xlab = "Number of Queries", ylab = "Standard deviation", col = "blue",
       main = NULL,
       lty = 1,
       ylim = c(1.1-(edu-1)*0.1, 1.6-(edu-1)*0.1),
       xaxs="i", yaxs="i"
       )
  lines(x = c(1:num_new_surveys), BARTResults$standardDeviationBART, type = 'l', col = "red", lty = 1)
  lines(x = c(1:num_new_surveys), BR$BRstandardDeviation, type = 'l', col = "green", lty = 1)
  legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling"),
         col=c("blue", "red", "green"), cex=0.8, lty = c(1,1))
  
  # 3. Close the file
  dev.off()
}

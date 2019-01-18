setwd("~/Documents/GitHub/BO-BART/src/meanPopulationStudy/")

library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)

# read in data
trainData <- read.csv("../../data/train.csv")
candidateData <- read.csv("../../data/candidate.csv")
populationData <- read.csv("../../data/full_data.csv")

studies <- c("White", "NA", "American Indian")

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args[1]

for (race in c(1,3)) {

  study <- studies[race]
  # extract covariates and response
  cols <- dim(trainData)[2] - 1
  trainX <- trainData[,-c(1, (cols+1))]
  trainY <- log(trainData[, (cols+1)])

  candidateX <- candidateData[,-c(1, (cols+1))]
  candidateY <- log(candidateData[, (cols+1)])

  # stratify the male population
  trainY <- trainY[trainX$Race == race]
  trainX <- trainX[trainX$Race == race, -8]

  candidateY <- candidateY[candidateX$Race == race]
  candidateX <- candidateX[candidateX$Race == race, -8]

  if (race == 3 & num_new_surveys > length(candidateY)){
    
    print(c("out of candidate size, maximum posible value has been set as", floor(length(candidateY) * 0.8)))
    
    num_new_surveys <- floor(length(candidateY) * 0.8)
    
  } else {

    trainX <- trainX[1:50, ]
    trainY <- trainY[1:50]

  }

  # compute the real population mean income
  poptMean <- mean(log(populationData[populationData$Race == race, ]$Total_person_income))

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
  jpeg("./%s.jpg" %--% c(study), width = 700, height = 583)
  # 2. Create the plot
  par(mfrow = c(1,2), pty = "s")
  plot(x = c(1:num_new_surveys), y = MImean,
       pch = 16, type = "l",
       xlab = "Number of Queries", ylab = "Population mean", col = "blue",
       main = "Mean income of %s vs N \nusing %s" %--% c(study, num_new_surveys),
       lty = 1,
       ylim = c(9-(race-1)/2, 11-(race-1)/2),
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
       main = "Standard Deviation income of %s vs N \nusing %s" %--% c(study, num_new_surveys),
       lty = 1,
       ylim = c(0.8+(race-1)*0.3, 1.4+(race-1)*0.3),
       xaxs="i", yaxs="i"
       )
  lines(x = c(1:num_new_surveys), BARTResults$standardDeviationBART, type = 'l', col = "red", lty = 1)
  lines(x = c(1:num_new_surveys), BR$BRstandardDeviation, type = 'l', col = "green", lty = 1)
  legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling"),
         col=c("blue", "red", "green"), cex=0.8, lty = c(1,1))
  
  # 3. Close the file
  dev.off()
}

setwd("~/Documents/GitHub/BO-BART/src/meanPopulationStudy/")
source("./bartMean.R")

train <- read.csv("train.csv")
candidate <- read.csv("candidate.csv")
population <- read.csv("full_data.csv")

train <- train[, -1]
candidate <- candidate[, -1]
population <- population[, -1]

par(mfrow = c(2, 2))

for (i in c(50, 200, 500, 1019)){
  
  trainX <- train[1:i, 1:8]
  trainY <- log(train[1:i, 9])
  candidateX <- candidate[, 1:8]
  candidateY <- log(candidate[, 9])
  
  # compute the real population mean income
  poptMean <- mean(log(population$Total_person_income))
  
  # Compute population mean
  BARTResults <- computePopulationMean(trainX, trainY, candidateX, candidateY, num_iterations = 10)
  
  MImean <- c()
  MIstandardDeviation <- c()
  for (j in 1:10) {
    MImean[j] <- mean(c(trainY, candidateY[1:j]))
    MIstandardDeviation[j] <- sqrt( var(c(trainY, candidateY[1:j])) ) 
  }
  
  BR <- BRcomputeMean(trainX, trainY, candidateX, candidateY, num_iterations = 10)
  
  plot(x = c(1:10), y = MImean * 10,
       pch = 16, type = "l",
       xlab = "Number of Queries", ylab = "Population mean", col = "blue",
       main = "Mean population income of vs N \nusing %s",
       lty = 1,
       ylim = c(70,120),
       xaxs="i", yaxs="i")
  lines(x = c(1:10), BARTResults$meanValueBART * 10, type = 'l', col = "red", lty = 1)
  lines(x = c(1:10), BR$BRmean * 10, type = 'l', col = "green", lty = 1)
  abline(a = poptMean * 10, b = 0, lty = 1, col = "black")
  #legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling", "Actual Mean"),
         #col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1))
  
}
setwd("~/Documents/GitHub/BO-BART/src/meanPopulationStudy/")
# read in data
trainData <- read.csv("../../data/train.csv")
candidateData <- read.csv("../../data/candidates.csv")
populationData <- read.csv("../../data/fullData.csv")
studies <- c("Male", "Female")

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args[1]

for (gender in c(1,2)) {

  study <- studies[gender]
  # extract covariates and response
  cols <- dim(trainData)[2] - 1
  trainX <- trainData[,-c(1, (cols+1))]
  trainY <- trainData[, (cols+1)]

  candidateX <- candidateData[,-c(1, (cols+1))]
  candidateY <- candidateData[, (cols+1)]

  # stratify the male population
  trainY <- trainY[trainX$SEX == gender]
  trainX <- trainX[trainX$SEX == gender,-3]

  candidateY <- candidateY[candidateX$SEX == gender]
  candidateX <- candidateX[candidateX$SEX == gender,-3]

  # run the BART regression model
  source("./bartMean.R")
  BARTResults<- computePopulationMean(trainX, trainY, candidateX, candidateY, num_iterations = num_new_surveys)

  # compute the real population mean income
  poptMean <- mean(populationData[populationData$SEX == gender,]$INCOME)

  MImean <- c()
  MIstandardDeviation <- c()
  for (i in 1:num_new_surveys) {
    MImean[i] <- mean(c(trainY, candidateY[1:i]))
    MIstandardDeviation[i] <- sqrt( var(c(trainY, candidateY[1:i])) )
  }

  results <- data.frame(
        "epochs" = c(1:num_new_surveys),
        "BARTMean" = BARTResults$meanValueBART, "BARTsd" = BARTResults$standardDeviationBART,
        "MIMean" = MImean, "MIsd" = MIstandardDeviation,
        "PoptMean" = poptMean
  )
  write.csv(results, file = "./%s.csv" %--% c(study),row.names=FALSE)

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
       ylim = c(0, 60000),
       xaxs="i", yaxs="i"
       )
  lines(x = c(1:num_new_surveys), BARTResults$meanValueBART, type = 'l', col = "red", lty = 1)
  abline(a = poptMean, b = 0, lty = 4)
  legend("topleft", legend=c("Monte Carlo", "BART", "Actual Mean"),
         col=c("blue", "red", "green"), cex=0.8, lty = c(1,1,1))
  
  plot(x = c(1:num_new_surveys), y = MIstandardDeviation,
       pch = 16, type = "l",
       xlab = "Number of Queries", ylab = "Standard deviation", col = "blue",
       main = "Standard Deviation income of %s vs N \nusing %s" %--% c(study, num_new_surveys),
       lty = 1,
       ylim = c(0, 70000),
       xaxs="i", yaxs="i"
       )
  lines(x = c(1:num_new_surveys), BARTResults$standardDeviationBART, type = 'l', col = "red", lty = 1)
  legend("topleft", legend=c("Monte Carlo", "BART"),
         col=c("blue", "red"), cex=0.8, lty = c(1,1))
  
  # 3. Close the file
  dev.off()
}

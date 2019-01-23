setwd("~/Documents/GitHub/BO-BART/src/meanPopulationStudy/")

library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(docstring)
source("./bartMean.R")
# set seed to enable reproduction
set.seed(1223)

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args[1]

# read in data
trainData <- read.csv("../../data/train.csv")
candidateData <- read.csv("../../data/remain.csv")
populationData <- read.csv("../../data/full_data.csv")

# reconstruct order of the data
trainData <- trainData[sample(nrow(trainData)), ]
candidateData <- candidateData[sample(nrow(candidateData)), ]

# extract covariates and response
cols <- dim(trainData)[2] - 1

trainX <- trainData[, -c(1, (cols+1))]
trainY <- trainData[, (cols+1)]/1000

candidateX <- candidateData[, -c(1, (cols+1))]
candidateY <- candidateData[, (cols+1)]/1000

# compute the real population mean income
poptMean <- mean(populationData$Total_person_income/1000)

# compute population average income estimates by BARTBQ
BARTResults <- computeBART(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys)

# population average income estimation by Monte Carlo
MIResults <- computeMI(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys)

# population average income estimation by block random sampling
BRSResults <- computeBRS(trainX, trainY, candidateX, candidateY, group = "Race", num_iterations=num_new_surveys)

# store results
results <- data.frame(
    "epochs" = c(1:num_new_surveys),
    "BARTMean" = BARTResults$meanValueBART, "BARTsd" = BARTResults$standardDeviationBART,
    "MIMean" = MIResults$meanValueMI, "MIsd" = MIResults$standardDeviationMI, 
    "BRSMean" = BRSResults$meanValueBRS, "BRSsd" = BRSResults$standardDeviationBRS, 
    "PoptMean" = poptMean
)
write.csv(results, file = "./population.csv", row.names=FALSE)

# 1. Open jpeg file
png("./population.png", width = 700, height = 583)
# 2. Create the plot
par(mfrow = c(1,2), pty = "s")
plot(x = c(1:num_new_surveys), y = MIResults$meanValueMI,
    pch = 16, type = "l",
    xlab = "Number of candidates added", ylab = "Average income (thousands of US$)", col = "blue",
    main = NULL,
    lty = 1,
    ylim = c(35, 70),
    xaxs="i", yaxs="i"
    )
lines(x = c(1:num_new_surveys), BARTResults$meanValueBART, type = 'l', col = "red", lty = 1)
lines(x = c(1:num_new_surveys), BRSResults$meanValueBRS, type = 'l', col = "green", lty = 1)
abline(a = poptMean, b = 0, lty = 1, col = "black")
legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling", "Actual Mean"),
        col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1))

plot(x = c(1:num_new_surveys), y = MIResults$standardDeviationMI,
    pch = 16, type = "l",
    xlab = "Number of candidates added", ylab = "Standard deviation", col = "blue",
    main = NULL,
    lty = 1,
    ylim = c(0, 15),
    xaxs="i", yaxs="i"
    )
lines(x = c(1:num_new_surveys), BARTResults$standardDeviationBART, type = 'l', col = "red", lty = 1)
lines(x = c(1:num_new_surveys), BRSResults$standardDeviationBRS, type = 'l', col = "green", lty = 1)
legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling"),
        col=c("blue", "red", "green"), cex=0.8, lty = c(1,1))

# 3. Close the file
dev.off()

# data segmentation by education level
studies <- c("Highschool and below", "Beyond highschool")
mat <- matrix(c(1,16,17,24), nrow = 2, ncol = 2)

for (cat in 1:2) {

  study <- studies[cat]
  
  # stratify the population by education
  eduTrainY <- trainY[trainX$Education >= mat[1, cat] & trainX$Education <= mat[2, cat]]
  eduTrainX <- trainX[(trainX$Education >= mat[1, cat] & trainX$Education <= mat[2, cat]), ]

  eduCandidateY <- candidateY[candidateX$Education >= mat[1, cat] & candidateX$Education <= mat[2, cat]]
  eduCandidateX <- candidateX[(candidateX$Education >= mat[1, cat] & candidateX$Education <= mat[2, cat]), ]

  # compute the real population mean income
  eduMean <- mean(populationData[(populationData$Education >= mat[1, cat] & populationData$Education <= mat[2, cat]), ]$Total_person_income/1000)
  
  # stratified population average income estimation by Monte Carlo
  eduMIResults <- computeMI(eduTrainX, eduTrainY, eduCandidateX, eduCandidateY, num_iterations=num_new_surveys)

  # stratified population average income estimation by block random sampling
  eduBRSResults <- computeBRS(eduTrainX, eduTrainY, eduCandidateX, eduCandidateY, "Race", num_iterations=num_new_surveys)

  # store results
  results <- data.frame(
      "epochs" = c(1:num_new_surveys),
      "BARTMean" = BARTResults$eduMeanValueBART[cat, ], "BARTsd" = BARTResults$eduStandardDeviationBART[cat, ],
      "MIMean" = eduMIResults$meanValueMI, "MIsd" = eduMIResults$standardDeviationMI, 
      "BRSMean" = eduBRSResults$meanValueBRS, "BRSsd" = eduBRSResults$standardDeviationBRS,
      "PoptMean" = eduMean
  )
  write.csv(results, file = "./%s.csv" %--% c(study), row.names=FALSE)

  # 1. Open jpeg file
  png("./%s.png" %--% c(study), width = 700, height = 583)
  # 2. Create the plot
  par(mfrow = c(1,2), pty = "s")
  plot(x = c(1:num_new_surveys), y = eduMIResults$meanValueMI, 
       pch = 16, type = "l",
       xlab = "Number of candidates added", ylab = "Average income (thousands of US$)", col = "blue",
       main = NULL,
       lty = 1,
       ylim = c(25+10*(cat-1), 70),
       xaxs="i", yaxs="i"
       )
  lines(x = c(1:num_new_surveys), BARTResults$eduMeanValueBART[cat, ], type = 'l', col = "red", lty = 1)
  lines(x = c(1:num_new_surveys), eduBRSResults$meanValueBRS, type = 'l', col = "green", lty = 1)
  abline(a = eduMean, b = 0, lty = 1)
  legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling",  "Actual Mean"),
         col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1))
  
  plot(x = c(1:num_new_surveys), y = eduMIResults$standardDeviationMI,
       pch = 16, type = "l",
       xlab = "Number of candidates added", ylab = "Standard deviation", col = "blue",
       main = NULL,
       lty = 1,
       ylim = c(0, 15),
       xaxs="i", yaxs="i"
       )
  lines(x = c(1:num_new_surveys), BARTResults$eduStandardDeviationBART[cat, ], type = 'l', col = "red", lty = 1)
  lines(x = c(1:num_new_surveys), eduBRSResults$standardDeviationBRS, type = 'l', col = "green", lty = 1)
  legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling"),
         col=c("blue", "red", "green"), cex=0.8, lty = c(1,1))
  
  # 3. Close the file
  dev.off()
}


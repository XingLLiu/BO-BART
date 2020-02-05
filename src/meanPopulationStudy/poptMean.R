setwd(getwd())

library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
source("./bartMean.R")
# set seed to enable reproduction
set.seed(0)

args <- as.double(commandArgs(TRUE))
num_new_surveys <- 100

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
BARTresults <- computeBART(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys)

# population average income estimation by Monte Carlo
MIresults <- computeMI(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys)

# population average income estimation by block random sampling
BRSresults <- computeBRS(trainX, trainY, candidateX, candidateY, group = "Race", num_iterations=num_new_surveys)

# store results
results <- data.frame(
    "epochs" = c(1:num_new_surveys),
    "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
    "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
    "BRSMean" = BRSresults$meanValueBRS, "BRSsd" = BRSresults$standardDeviationBRS, 
    "PoptMean" = poptMean
)
write.csv(results, file = "./population.csv", row.names=FALSE)

real <- results$PoptMean[1]
# 1. Open jpeg file
plot_points <- seq(0, num_new_surveys, 10)
pdf("results.pdf", width = 8,5, height = 10)
par(mfrow = c(1,2), pty = "s")
ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$BRSMean - 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$BRSMean + 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
plot(results$epochs,
     abs(results$BARTMean - results$PoptMean),
     ty="l",
     ylab = "Absolute Error",
     xlab = "num_iterations",
     col = "chartreuse4",
     ylim = c(0, 40)
)
abline(h=0)
points(results$epochs[plot_points], abs(results$MIMean[plot_points] - real), col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
points(results$epochs[plot_points], abs(results$BARTMean[plot_points] - real), col = "orangered",bg='orangered', pch=21, lwd=3)
points(results$epochs, abs(results$BARTMean - real), ty="l", col = "orangered")
points(results$epochs, abs(results$BRSMean-real), ty="l", col = "dodgerblue")
points(results$epochs[plot_points], abs(results$BRSMean[plot_points] - real), col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
legend("topright", legend=c("Block sampling", "BART-Int", "GP-BQ"),
       col=c("chartreuse4", "orangered", "dodgerblue"), cex=0.8, lty = c(1,1,1,1))

ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$BRSMean - 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$BRSMean + 2*results$BRSsd, results$MIMean[1:num_new_surveys]))

plot(results$epochs, 
     results$MIMean, 
     ty="l", 
     ylab = "Mean Population",
     xlab = "num_iterations",
     col = "chartreuse4",
     ylim = c(ymin, ymax)
)
points(results$epochs[plot_points], results$MIMean[plot_points], col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
points(results$epochs, results$BRSMean, ty="l", col = "dodgerblue")
points(results$epochs[plot_points], results$BRSMean[plot_points], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
polygon(c(results$epochs, rev(results$epochs)), 
        c(
          results$BRSMean + 2*results$BRSsd, 
          rev(results$BRSMean - 2*results$BRSsd)
        ), 
        col = adjustcolor("dodgerblue", alpha.f = 0.10), 
        border = "dodgerblue", lty = c("dashed", "solid"))
points(results$epochs, results$BARTMean, ty="l", col = "orangered")
points(results$epochs[plot_points], results$BARTMean[plot_points], col = "orangered",bg='orangered', pch=21, lwd=3)
polygon(c(results$epochs, rev(results$epochs)), 
        c(
          results$BARTMean + 2*results$BARTsd, 
          rev(results$BARTMean - 2*results$BARTsd)
        ), 
        col = adjustcolor("orangered", alpha.f = 0.10), 
        border = "orangered", lty = c("dashed", "solid"))
abline(h=results$PoptMean)
dev.off()

# data segmentation by education level
# studies <- c("Highschool and below", "Beyond highschool")
# mat <- matrix(c(1,16,17,24), nrow = 2, ncol = 2)
# 
# for (cat in 1:2) {
# 
#   study <- studies[cat]
#   
#   # stratify the population by education
#   eduTrainY <- trainY[trainX$Education >= mat[1, cat] & trainX$Education <= mat[2, cat]]
#   eduTrainX <- trainX[(trainX$Education >= mat[1, cat] & trainX$Education <= mat[2, cat]), ]
# 
#   eduCandidateY <- candidateY[candidateX$Education >= mat[1, cat] & candidateX$Education <= mat[2, cat]]
#   eduCandidateX <- candidateX[(candidateX$Education >= mat[1, cat] & candidateX$Education <= mat[2, cat]), ]
# 
#   # compute the real population mean income
#   eduMean <- mean(populationData[(populationData$Education >= mat[1, cat] & populationData$Education <= mat[2, cat]), ]$Total_person_income/1000)
#   
#   # stratified population average income estimation by Monte Carlo
#   eduMIresults <- computeMI(eduTrainX, eduTrainY, eduCandidateX, eduCandidateY, num_iterations=num_new_surveys)
# 
#   # stratified population average income estimation by block random sampling
#   eduBRSresults <- computeBRS(eduTrainX, eduTrainY, eduCandidateX, eduCandidateY, "Race", num_iterations=num_new_surveys)
# 
#   # store results
#   results <- data.frame(
#       "epochs" = c(1:num_new_surveys),
#       "BARTMean" = BARTresults$eduMeanValueBART[cat, ], "BARTsd" = BARTresults$eduStandardDeviationBART[cat, ],
#       "MIMean" = eduMIresults$meanValueMI, "MIsd" = eduMIresults$standardDeviationMI, 
#       "BRSMean" = eduBRSresults$meanValueBRS, "BRSsd" = eduBRSresults$standardDeviationBRS,
#       "PoptMean" = eduMean
#   )
#   write.csv(results, file = "./%s.csv" %--% c(study), row.names=FALSE)
# 
#   # 1. Open jpeg file
#   png("./%s.png" %--% c(study), width = 700, height = 583)
#   # 2. Create the plot
#   par(mfrow = c(1,2), pty = "s")
#   plot(x = c(1:num_new_surveys), y = eduMIresults$meanValueMI, 
#        pch = 16, type = "l",
#        xlab = "Number of candidates added", ylab = "Average income (thousands of US$)", col = "blue",
#        main = NULL,
#        lty = 1,
#        ylim = c(25+10*(cat-1), 70),
#        xaxs="i", yaxs="i"
#        )
#   lines(x = c(1:num_new_surveys), BARTresults$eduMeanValueBART[cat, ], type = 'l', col = "red", lty = 1)
#   lines(x = c(1:num_new_surveys), eduBRSresults$meanValueBRS, type = 'l', col = "green", lty = 1)
#   abline(a = eduMean, b = 0, lty = 1)
#   legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling",  "Actual Mean"),
#          col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1))
#   
#   plot(x = c(1:num_new_surveys), y = eduMIresults$standardDeviationMI,
#        pch = 16, type = "l",
#        xlab = "Number of candidates added", ylab = "Standard deviation", col = "blue",
#        main = NULL,
#        lty = 1,
#        ylim = c(0, 15),
#        xaxs="i", yaxs="i"
#        )
#   lines(x = c(1:num_new_surveys), BARTresults$eduStandardDeviationBART[cat, ], type = 'l', col = "red", lty = 1)
#   lines(x = c(1:num_new_surveys), eduBRSresults$standardDeviationBRS, type = 'l', col = "green", lty = 1)
#   legend("topleft", legend=c("Monte Carlo", "BART", "Block sampling"),
#          col=c("blue", "red", "green"), cex=0.8, lty = c(1,1))
#   
#   # 3. Close the file
#   dev.off()
# }


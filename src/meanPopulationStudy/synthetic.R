# library(dbarts)
# library(data.tree)
# mu <- data.frame(matrix(0, ncol = 7, nrow = 7))
# grid <- data.frame(matrix(0, ncol = 2, nrow = 0))
# sd <- 0.1
# samples <- data.frame(matrix(0, ncol = 3, nrow = 0))
# colnames(samples) <- c("i", "j", "count")
# n_samples <- 50
# for (i in 1:7) {
#   for (j in 1:7) {
#     mu[i, j] = 2
#     grid <- rbind(grid, c(i, j))
#     df <- data.frame(
#       "i"= rep(i, n_samples),
#       "j"= rep(j, n_samples),
#       "count"= rnorm(n_samples, mean=mu[i, j], sd=1)
#     )
#     samples <- rbind(samples, df)
#   }
# }
# model <- bart(samples[, c(1, 2)], samples$count, keeptrees=TRUE, keepevery=3L, 
#      nskip=500, ndpost=5000, ntree=60, k=3, usequant=FALSE)      
# 
# fValues <- predict(model, as.matrix(grid))
# means <- colMeans(fValues)
# means
# dim(fValues)
library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(caret)
source("src/meanPopulationStudy/bartMean.R")
source("src/meanPopulationStudy/gpMean_2.R")
source("src/optimise_gp.R")
# paths to save results and plots
resultPath <- "results/"
plotPath <- "Figures/"

set.seed(0)
num_new_surveys <- 100
num_variables <- 10
num_data <- 2500
df <- data.frame(matrix(0, ncol = num_variables, nrow = num_data))
colnames(df) <- c(1:num_variables)
p = 10
for (variable in 1:num_variables) {
  # num_classes <- sample(c(2:10), 1)
  num_classes <- 2
  print(num_classes)
  data <- sample.int(num_classes, size=num_data, replace=TRUE) - 1
  df[, variable] <- as.factor(data)
}


## generate gaussian via one-hot encoding
dummyFullData <- dummyVars("~.", data = df)
trainX <- data.frame(predict(dummyFullData, newdata = df[1:500,]))
candidateX <- data.frame(predict(dummyFullData, newdata = df[501:2500,]))
# trainX <- df[1:500,]
# candidateX <- df[501:2500,]
coeffs <- rnorm(p, mean=0, 1)
trainY <- rowSums(trainX * coeffs)#  + rnorm(dim(trainX)[1], 0,1) ;X\beta + N(0, 1) 
candidateY <- rowSums(candidateX * coeffs) # + rnorm(dim(candidateX)[1], 0,1); X\beta + N(0, 1)
real <- 2^p / (2^p - 1) * sum(coeffs) # 2^p / (2^p - 1) \sum_{c=1}^p \beta_p

for (num_cv in 2:20) {
  # set new seed
  set.seed(num_cv)
  print(num_cv)
  # compute population average income estimates by BARTBQ
  BARTresults <- computeBART(trainX, trainY, as.matrix(candidateX), candidateY, num_iterations=num_new_surveys)
  
  # population average income estimation by Monte Carlo
  MIresults <- computeMI(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys, seed=num_cv)
  
  # GPBQ
  if (num_cv == 1) {
    lengthscale <- optimise_gp_r(as.matrix(trainX), trainY, kernel = "rbf", epochs = 500)
  }
  else {
    lengthscale=lengthscale # change it
  }
  GPresults <- computeGPBQEmpirical(as.matrix(trainX), trainY, as.matrix(candidateX), candidateY, epochs=num_new_surveys, lengthscale=lengthscale)
  
  plot(BARTresults$meanValueBART, ty="l", ylim=c(1.5, 1.8)); abline(h=real)
  plot(MIresults$meanValueMI, ty="l"); abline(h=real)
  plot(GPresults$meanValueGP, ty="l"); abline(h=real)
  # population average income estimation by block random sampling
  # BRSresults <- computeBRS(trainX.num, trainY, candidateX.num, candidateY, group = "Race", num_iterations=num_new_surveys)
  
  # store results
  # results <- data.frame(
  #   "epochs" = c(1:num_new_surveys),
  #   "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
  #   "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
  #   "BRSMean" = BRSresults$meanValueBRS, "BRSsd" = BRSresults$standardDeviationBRS, 
  #   "PoptMean" = poptMean, "BpoptMean" = BARTpoptMean
  # )
  results <- data.frame(
    "epochs" = c(1:num_new_surveys),
    "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
    "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
    "GPMean" = GPresults$meanValueGP, "GPsd" = GPresults$varianceGP,
    "real" = real
  )
  write.csv(results, file = paste0(resultPath, "syntheticCategorical", num_cv, ".csv"), row.names=FALSE)
  results_models <- list("BART"=BARTresults, "MI"=MIresults, "GP"=GPresults)
  save(results_models, file = paste0(plotPath, "syntheticCategorical", num_cv, ".RData"))
  
  print(c("Real", real))
  print(c("BART-Int", BARTresults$meanValueBART[num_new_surveys]))
  print(c("MI", MIresults$meanValueMI[num_new_surveys]))
  print(c("GP-BQ", GPresults$meanValueGP[num_new_surveys]))
  
  # 1. Open jpeg file
  # pdf(paste0(plotPath, "syntheticCategorical", num_cv, ".pdf"), width = 8, height = 10)
  # par(mfrow = c(1,2), pty = "s")
  # # ymax <- max(c(abs(results$BARTMean - real), abs(results$BRSMean - real), abs(results$MIMean - real)))
  # ymax <- max(c(abs(results$BARTMean - real), abs(results$GPMean - real)))
  # plot(results$epochs,
  #      abs(results$BARTMean - results$PoptMean),
  #      ty="l",
  #      ylab = "Absolute Error",
  #      xlab = "num_iterations",
  #      col = "orangered",
  #      ylim = c(0, ymax)
  # )
  # abline(h=0)
  # 
  # points(results$epochs, abs(results$MIMean - real), ty="l", col = "chartreuse4")
  # points(results$epochs, abs(results$GPMean - real), ty="l", col = "darkgoldenrod")
  # 
  # legend("topright", legend=c("BART-Int", "MI", "GPBQ"),
  #        col=c("orangered", "chartreuse4", "darkgoldenrod"), cex=0.8, lty = c(1,1,1,1))
  # 
  # # ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$BRSMean - 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
  # # ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$BRSMean + 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
  # ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$GPMean - 2*results$GPsd))
  # ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$GPMean + 2*results$GPsd))    
  # 
  # plot(results$epochs, 
  #      results$MIMean, 
  #      ty="l", 
  #      ylab = "Mean",
  #      xlab = "num_iterations",
  #      col = "chartreuse4",
  #      ylim = c(ymin, ymax)
  # )
  # polygon(c(results$epochs, rev(results$epochs)), 
  #         c(
  #           results$MIMean + 2*results$MIsd, 
  #           rev(results$MIMean - 2*results$MIsd)
  #         ), 
  #         col = adjustcolor("chartreuse4", alpha.f = 0.10), 
  #         border = "chartreuse4", lty = c("dashed", "solid"))
  # points(results$epochs, results$BRSMean, ty="l", col = "dodgerblue")
  # polygon(c(results$epochs, rev(results$epochs)), 
  #         c(
  #           results$BRSMean + 2*results$BRSsd, 
  #           rev(results$BRSMean - 2*results$BRSsd)
  #         ), 
  #         col = adjustcolor("dodgerblue", alpha.f = 0.10), 
  #         border = "dodgerblue", lty = c("dashed", "solid"))
  # points(results$epochs, results$BARTMean, ty="l", col = "orangered")
  # polygon(c(results$epochs, rev(results$epochs)), 
  #         c(
  #           results$BARTMean + 2*results$BARTsd, 
  #           rev(results$BARTMean - 2*results$BARTsd)
  #         ), 
  #         col = adjustcolor("orangered", alpha.f = 0.10), 
  #         border = "orangered", lty = c("dashed", "solid"))
  # 
  # points(results$epochs, results$GPMean, ty="l", col = "darkgoldenrod")
  # polygon(c(results$epochs, rev(results$epochs)),
  #         c(
  #           results$BARTMean + 2*results$GPsd, 
  #           rev(results$BARTMean - 2*results$GPsd)
  #         ), 
  #         col = adjustcolor("darkgoldenrod", alpha.f = 0.10), 
  #         border = "darkgoldenrod", lty = c("dashed", "solid"))
  # 
  # abline(h=results$PoptMean)
  # dev.off()
}

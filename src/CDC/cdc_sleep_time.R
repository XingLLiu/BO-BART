library(foreign)
df <- read.xport("data/brfss/LLCP2017.XPT ")
df <- df[!is.na(df$SLEPTIM1), ]
df <- df[df$SLEPTIM1<=24, ]

# Employment","Sex","Education","Disability","Health_insurance","Own_child","Race"
cols <- c("SEX", "MARITAL", "EDUCA", "VETERAN3", "EMPLOY1", "DEAF", "BLIND", "SMOKE100", "SLEPTIM1")
sleep_time <- df[, "SLEPTIM1"]
cdc_data <- df[ ,cols]
for (i in 1:8){ cdc_data[,i] <- as.factor(cdc_data[, i]) }
colnames(cdc_data) <- cols
cdc_data <- cdc_data[complete.cases(cdc_data),]

library(lhs)
library(dbarts)
library(data.tree)
library(matrixStats)
library(caret)
set.seed(0)

source("src/meanPopulationStudy/bartMean2.R")
source("src/meanPopulationStudy/gpMean_2.R")
source("src/optimise_gp.R")

# paths to save results and plots
resultPath <- "results/cdc/"
plotPath <- "Figures/cdc/"

args <- as.double(commandArgs(TRUE))
num_new_surveys <- args[1]
num_cv_start <- args[2]
num_cv_end <- args[3]
num_data <- args[4]   # set to 2000 for this lengthscale
num_design <- args[5]

num_training <- 10000
set.seed(123)
training_indices <- sample.int(nrow(cdc_data), num_training)

trainData_full <- cdc_data[training_indices,]
candidateData_full <- cdc_data[-training_indices,]

# compute the real population mean sleep time
set.seed(1123)
data <- cdc_data[sample.int(nrow(cdc_data), 10000),]
dim <- ncol(data)
model <- bart(data[,1:(dim-1)], data[, dim], keeptrees=TRUE, keepevery=3L,
                nskip=1000, ndpost=10000, ntree=50, k=3, usequant=FALSE)
poptMean <- mean(model$yhat.train.mean)

for (num_cv in num_cv_start:num_cv_end) {
  
    # set new seed
    set.seed(num_cv)
    print(num_cv)
    trainData <- trainData_full[sample(c(1:dim(trainData_full)[1]), num_design), ]
    candidateData <- candidateData_full[1:num_data, ]
    
    cols <- ncol(trainData)
    
    trainX <- trainData[, -cols]
    trainY <- trainData[, cols]
    candidateX <- candidateData[, -cols]
    candidateY <- candidateData[, cols]

    # one-hot encoding
    trainX.num <- trainX
    candidateX.num <- candidateX
    dummyFullData <- dummyVars("~.", data = rbind(trainX, candidateX))
    trainX <- data.frame(predict(dummyFullData, newdata = trainX))
    candidateX <- data.frame(predict(dummyFullData, newdata = candidateX))
    # save(trainX, trainX.num, trainY, candidateX, candidateX.num, candidateY, file = "data/survey_data.RData")
    # load(file = "data/survey_data.RData")
    # compute population average income estimates by BARTBQ
    t0 <- proc.time()
    BARTresults <- computeBART(trainX, trainY, candidateX, candidateY, num_iterations=num_new_surveys)
    t1 <- proc.time()
    bartTime <- (t1 - t0)[[1]]
    # population average income estimation by Monte Carlo
    # MIresults <- computeMI(trainX.num, trainY, candidateX.num, candidateY, num_iterations=num_new_surveys)
    MIresults <- computeMI(trainX.num, trainY, candidateX.num, candidateY, num_iterations=num_new_surveys, seed = num_cv)
    # plot(MIresults$meanValueMI, ylim = c(10.9, 11.1), xlab = "num_iterations", ylab = "mean population")
    # legend("topright", legend=c("MC integration"))
    # abline(h = poptMean, col = "red")
    
    # GPBQ
    lengthscale <- optimise_gp_r(as.matrix(trainX), trainY, kernel = "rbf", epochs = 500)
    # lengthscale <- 3.35
    t0 <- proc.time()
    GPresults <- computeGPBQEmpirical(as.matrix(trainX), trainY, as.matrix(candidateX), candidateY, epochs=num_new_surveys, lengthscale=lengthscale)
    t1 <- proc.time()
    GPTime <- (t1 - t0)[[1]]
    
    # store results
    results <- data.frame(
         "epochs" = c(1:num_new_surveys),
         "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
         "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
         "GPMean" = GPresults$meanValueGP, "GPsd" = GPresults$varianceGP, "PoptMean" = poptMean,
         "runtimeBART" = rep(bartTime, num_new_surveys),
         "runtimeGP" = rep(GPTime, num_new_surveys)
     )
    # results <- data.frame(
    #     "epochs" = c(1:num_new_surveys),
    #     "BARTMean" = BARTresults$meanValueBART, "BARTsd" = BARTresults$standardDeviationBART,
    #     "MIMean" = MIresults$meanValueMI, "MIsd" = MIresults$standardDeviationMI, 
    #     "BRSMean" = BRSresults$meanValueBRS, "BRSsd" = BRSresults$standardDeviationBRS,
    #     "GPMean" = GPresults$meanValueGP, "GPsd" = GPresults$varianceGP,
    #     "PoptMean" = poptMean
    # )
	  #results <- data.frame("epochs"=c(1:num_new_surveys), "GPMean"=GPresults$meanValueGP, "GPsd"=GPresults$varianceGP)
    write.csv(results, file = paste0(resultPath, "results", num_cv, ".csv"), row.names=FALSE)
    results_models <- list("BART"=BARTresults, "MI"=MIresults, "GP"=GPresults)
    save(results_models, file = paste0(plotPath, "results", num_cv, ".RData"))
    
    real <- results$PoptMean[1]
    
    print(c("Real MC: ", real))
    print(c("BART-Int: ", results$BARTMean[num_new_surveys]))
    print(c("GP-BQ: ", results$GPMean[num_new_surveys]))
    print(c("MI: ", results$MIMean[num_new_surveys]))

    # 1. Open jpeg file
    pdf(paste0(plotPath, "result", num_cv, ".pdf"), width = 8, height = 10)
    par(mfrow = c(1,2), pty = "s")
    # ymax <- max(c(abs(results$BARTMean - real), abs(results$BRSMean - real), abs(results$MIMean - real)))
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

    points(results$epochs, abs(results$MIMean - real), ty="l", col = "chartreuse4")
    points(results$epochs, abs(results$GPMean - real), ty="l", col = "dodgerblue")
    
    legend("topright", legend=c("BART-Int", "MI", "GPBQ"),
            col=c("orangered", "chartreuse4", "dodgerblue"), cex=0.8, lty = c(1,1,1,1))
    
    # ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$BRSMean - 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
    # ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$BRSMean + 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
    ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$MIMean - 2*results$MIsd))
    ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$MIMean + 2*results$MIsd))    

    plot(results$epochs, 
        results$MIMean, 
        ty="l", 
        ylab = "Mean Sleep Time",
        xlab = "num_iterations",
        col = "chartreuse4",
        ylim = c(ymin, ymax)
    )
    # polygon(c(results$epochs, rev(results$epochs)), 
    #         c(
    #             results$MIMean + 2*results$MIsd, 
    #             rev(results$MIMean - 2*results$MIsd)
    #         ), 
    #         col = adjustcolor("chartreuse4", alpha.f = 0.10), 
            # border = "chartreuse4", lty = c("dashed", "solid"))
    points(results$epochs, results$BARTMean, ty="l", col = "orangered")
    polygon(c(results$epochs, rev(results$epochs)), 
            c(
                results$BARTMean + 2*results$BARTsd, 
                rev(results$BARTMean - 2*results$BARTsd)
            ), 
            col = adjustcolor("orangered", alpha.f = 0.10), 
            border = "orangered", lty = c("dashed", "solid"))
    points(results$epochs, results$GPMean, ty="l", col = "dodgerblue")
    polygon(c(results$epochs, rev(results$epochs)),
            c(
                results$GPMean + 2*results$GPsd, 
                rev(results$GPMean - 2*results$GPsd)
            ), 
            col = adjustcolor("dodgerblue", alpha.f = 0.10), 
            border = "dodgerblue", lty = c("dashed", "solid"))
    
    abline(h=results$PoptMean)
    dev.off()
}
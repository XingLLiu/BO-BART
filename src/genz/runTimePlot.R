# Plot run time for each dimension
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To see results:
#     Plots stored in {ROOT}/report/Figures for plots


methods <- c("BART", "MI", "GP")
bestMethod <- matrix(NA, ncol = 6, nrow = 6)
rownames(bestMethod) <- c("cont", "copeak", "disc", "gaussian", "oscil", "prpeak")
colnames(bestMethod) <- c("1", "2", "3", "5", "10", "20")

# Initialize BART run time
BARTRunTimeAll <- matrix(NA, nrow = 1, ncol = 6)

# Retrieve GP run time
GPRunTimeAll <- read.csv("../results/genz/GPRunTime/GPTime.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

for (i in 1:6){
    for (j in 1:6){

        # global parameters: dimension
        dimensionsList <- c(1,2,3,5,10,20)
        args <- cbind(dimensionsList, c(1:6))
        dim <- args[j, 1]
        whichGenz <- args[i, 2]
        num_iterations <- 500
        
        # Skil if dim = 1 for discontinuous
        if (whichGenz == 3 & dim == 1) { next } 

        # Find Genz function
        if (whichGenz == 1) { genzFunctionName <-  "cont" }
        if (whichGenz == 2) { genzFunctionName <-  "copeak" }
        if (whichGenz == 3) { genzFunctionName <-  "disc" }
        if (whichGenz == 4) { genzFunctionName <-  "gaussian" }
        if (whichGenz == 5) { genzFunctionName <-  "oscil" }
        if (whichGenz == 6) { genzFunctionName <-  "prpeak" }

        # Set path for estimated integral values
        fileName <- paste('results', toString(whichGenz), 'dim', toString(dim), '.csv', sep='')
        filePath <- paste('../results/genz', toString(whichGenz), fileName, sep='/')

        # Retrieve BART run time
        integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
        BARTRunTime <- integrals[1, 9]

        # Retrieve GPBQ run time
        GPRunTime <- GPRunTimeAll[1, (8 * i + (j - 1))]
    
        # 1. Open jpeg file
        plotName <- paste("convergenceMean", toString(whichGenz), toString(dim), "Dimensions", ".eps", sep = "")
        plotPath <- gsub(paste("../report/Figures/", toString(whichGenz), "/", plotName, sep=""), pattern = "csv", replacement="jpg")
        postscript(plotPath, width = 2100, height = 1794)

        # 2. Create the plot
        # Set y limits
        yLimMean <- cbind(quantile(predictionBART$meanValueBART, probs=c(0.1, 1), na.rm=TRUE),
                       quantile(predictionMonteCarlo$meanValueMonteCarlo, probs=c(0.1, 1), na.rm=TRUE))
        
        yLowLimMean <- min(yLimMean[, 1])
        if (yLowLimMean < 0) {scalingFactor <- 2}
        else {scalingFactor <- -1}
        yLowLimMean <-  yLowLimMean * scalingFactor

        yHighLimMean <- max(yLimMean[, 2])
        if (yHighLimMean > 0) {scalingFactor <- 2}
        else {scalingFactor <- -1}
        yHighLimMean <-  yHighLimMean * scalingFactor
        yHighLimMean <- max(yLimMean[, 2]) * scalingFactor
        
        par(mfrow = c(1,2), pty = "s", mar=c(2, 4, 4, 7))
        plot(x = c(1:num_iterations), y = predictionMonteCarlo$meanValueMonteCarlo,
            pch = 16, type = "l",
            xlab = "Number of epochs N", ylab = "Integral approximation", col = "blue",
            main = paste("Convergence of methods: mean vs N \nusing", genzFunctionName, "with", num_iterations, "epochs in", dim, "dim"),
            # ylim = c(-real, real + real),
            ylim = c(yLowLimMean, yHighLimMean), 
            lty = 1,
            xaxs="i", yaxs="i"
            )
        lines(x = c(1:num_iterations), predictionBART$meanValueBART, type = 'l', col = "red", lty = 1)
        lines(x = c(1:num_iterations), predictionGPBQ$meanValueGP, type = 'l', col = "green", lty = 1)
        abline(a = real, b = 0, lty = 4)
        legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ", "Actual"),
            col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1,1))

        # 2. Create the plot
        # Set y limits
        yLimSd <- cbind(quantile(log(predictionBART$standardDeviationBART), probs=c(0.1, 0.9), na.rm=TRUE),
                       quantile(log(predictionMonteCarlo$standardDeviationMonteCarlo), probs=c(0.1, 0.9), na.rm=TRUE))
        
        yLowLimSd <- min(yLimSd[, 1])
        if (yLowLimSd < 0) {scalingFactor <- 1.5}
        else {scalingFactor <- 0.2}
        yLowLimSd <-  yLowLimSd * scalingFactor

        yHighLimSd <- max(yLimSd[, 2])
        if (yHighLimSd > 0) {scalingFactor <- 2}
        else {scalingFactor <- 0}
        yHighLimSd <-  yHighLimSd * scalingFactor
        yHighLimSd <- max(yLimSd[, 2]) * scalingFactor

        plot(x = log(c(2:num_iterations)), y = log(predictionMonteCarlo$standardDeviationMonteCarlo[-1]),
            pch = 16, type = "l",
            xlab = "Log number of epochs N", ylab = "Log standard deviation", col = "blue",
            main = paste("Convergence of methods: log(sd) vs log(N) \nusing", genzFunctionName, "with", num_iterations, "epochs in", dim, "dim"),
            lty = 1,
            ylim = c(yLowLimSd, yHighLimSd),
            xaxs="i", yaxs="i")
        lines(x = log(c(2:num_iterations)), log(predictionBART$standardDeviationBART[-1]), type = 'l', col = "red", lty = 1)
        # lines(x = log(c(2:num_iterations)), log(predictionGPBQ$standardDeviationGP[-1]), type = 'l', col = "green", lty = 1)

        # legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ"),
        #     col=c("blue", "red", "green"), cex=0.8, lty = c(1,1,1,1))

        legend("topleft", legend=c("MC Integration", "BART BQ"),
              col=c("blue", "red"), cex=0.8, lty = c(1,1,1,1))   
        # 3. Close the file
        dev.off()

        # Compute abs error
        epoch <- nrow(integrals)
        absDifference <- abs(c(predictionBART$meanValueBART[epoch], predictionMonteCarlo$meanValueMonteCarlo[epoch], 
                           predictionGPBQ$meanValueGP[epoch]) - real)
        bestMethod[i, j] <- methods[which(absDifference == min(absDifference))[1]]
    }
    cat(genzFunctionName, "done", '\n')
}

write.csv(bestMethod, file = "../results/genz/bestMethods.csv")

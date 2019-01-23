# Plot results for each Genz integral after 500 epochs.
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To manually adjust y limits of the plots:
#     Uncomment ylim = ylims[1:2] and ylim = ylims[3:3] in the two plot functions, 
#     and input the four limits when running the file
# To see results:
#     Plots stored in {ROOT}/report/Figures for plots

ylims <- as.double(commandArgs(TRUE))

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

        # Retrieve estimated integral values
        integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
        predictionBART <- data.frame("meanValueBART" = integrals[, 2], "standardDeviationBART" = integrals[, 3])
        predictionMonteCarlo <- data.frame("meanValueMonteCarlo" = integrals[, 4], "standardDeviationMonteCarlo" = integrals[, 5])
        predictionGPBQ <-  data.frame("meanValueGP" = integrals[, 6], "standardDeviationGP" = integrals[, 7])

        # Retrieve analytical integral values
        whichDimension <- which(dim == dimensionsList)
        analyticalIntegrals <- read.csv("./genz/integrals.csv", header = FALSE)
        real <- analyticalIntegrals[whichGenz, whichDimension]

        # 1. Open eps file
        # plotName <- paste("convergenceMean", toString(whichGenz), toString(dim), "Dimensions", sep = "")
        plotName <- paste("convergenceMean", toString(whichGenz), toString(dim), "Dimensions", ".eps", sep = "")
        plotPath <- gsub(paste("../report/Figures/", toString(whichGenz), "/", plotName, sep=""), pattern = "csv", replacement="jpg")
        postscript(plotPath, width = 2800, height = 1794)

        # 2. Create the plot
        # Set y limits
        yLimMean <- cbind(quantile(predictionBART$meanValueBART, probs=c(0, 1), na.rm=TRUE),
                       quantile(predictionMonteCarlo$meanValueMonteCarlo, probs=c(0, 1), na.rm=TRUE))
        
        yLowLimMean <- min(yLimMean[, 1])
        if (yLowLimMean < 0) {scalingFactor <- 1.1}
        else {scalingFactor <- 0.5}
        yLowLimMean <-  yLowLimMean * scalingFactor

        yHighLimMean <- max(yLimMean[, 2])
        if (yHighLimMean > 0) {scalingFactor <- 0.8}
        else {scalingFactor <- 0.9}
        yHighLimMean <-  yHighLimMean * scalingFactor
        yHighLimMean <- yHighLimMean * scalingFactor
        
        par(mfrow = c(1,2), pty = "s", mar=c(2, 5, 2, 5))
        plot(x = c(1:num_iterations), y = predictionMonteCarlo$meanValueMonteCarlo,
            pch = 16, type = "l",
            xlab = "Number of epochs N", ylab = "Integral estimates", col = "blue",
            ylim = c(yLowLimMean, yHighLimMean), 
            # ylim = ylims[1:2],
            xaxs="i", yaxs="i",
            cex.lab = 2.0,
            cex.axis = 1.5
            )
        lines(x = c(1:num_iterations), predictionBART$meanValueBART, type = 'l', col = "red", lty = 1)
        lines(x = c(1:num_iterations), predictionGPBQ$meanValueGP, type = 'l', col = "green", lty = 1)
        abline(a = real, b = 0, lty = 4)
        axis(4, at = signif(real, digits=4), las = 2, cex.axis = 1.2)
        legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ", "Actual"),
            col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1,1), bg="white")

        # 2. Create the plot
        # Set y limits
        yLimSd <- cbind(quantile(log(predictionBART$standardDeviationBART), probs=c(0.1, 0.9), na.rm=TRUE),
                       quantile(log(predictionMonteCarlo$standardDeviationMonteCarlo), probs=c(0.1, 0.9), na.rm=TRUE))
        
        yLowLimSd <- min(yLimSd[, 1])
        if (yLowLimSd < 0) {scalingFactor <- 1}
        else {scalingFactor <- -0.2}
        yLowLimSd <-  yLowLimSd * scalingFactor

        yHighLimSd <- max(yLimSd[, 2])
        if (yHighLimSd > 0) {scalingFactor <- 1}
        else {scalingFactor <- -1}
        yHighLimSd <-  yHighLimSd * scalingFactor
        yHighLimSd <- yHighLimSd * scalingFactor

        plot(x = log(c(2:num_iterations)), y = log(predictionMonteCarlo$standardDeviationMonteCarlo[-1]),
            pch = 16, type = "l",
            xlab = "Log number of epochs N", ylab = "Log standard deviation", col = "blue",
            lty = 1,
            ylim = c(yLowLimSd, yHighLimSd),
            # ylim = ylims[3:4],
            xaxs="i", yaxs="i",
            cex.lab = 2.0,
            cex.axis = 1.5
            )
        lines(x = log(c(2:num_iterations)), log(predictionBART$standardDeviationBART[-1]), type = 'l', col = "red", lty = 1)
        legend("topleft", legend=c("MC Integration", "BART BQ"), bg="white",
              col=c("blue", "red"), cex=0.8, lty = c(1,1,1,1), pt.cex = 1.5)   
        # 3. Close the file
        dev.off()

    }
    cat(genzFunctionName, "done", '\n')
}

# Plot results for Genz integral estimates compared with GPBQ with less sequential design samples.
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To see results:
#     Plots stored in {ROOT}/report/Figures for plots


methods <- c("BART", "MI", "GP")
bestMethod <- matrix(NA, ncol = 6, nrow = 6)
rownames(bestMethod) <- c("cont", "copeak", "disc", "gaussian", "oscil", "prpeak")
colnames(bestMethod) <- c("1", "2", "3", "5", "10", "20")

for (i in 5:5){
    for (j in 5:5){


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
    
        # Retrieve GP-without-sequential-design estimates
        GPBQpath <- paste('../results/genz/GPBQSequentialDesign/results', toString(whichGenz), "dim", toString(dim), "gp.csv", sep="")
        results <- (read.csv(GPBQpath, header=TRUE, stringsAsFactors = FALSE))[1:num_iterations, ]
        predictionGPBQNoSeqDesign <- data.frame("meanValueGPNoSeq" = results[, 1], "standardDeviationGPNoSeq" = results[, 2])

        # Retrieve analytical integral values
        whichDimension <- which(dim == dimensionsList)
        analyticalIntegrals <- read.csv("./genz/integrals.csv", header = FALSE)
        real <- analyticalIntegrals[whichGenz, whichDimension]
        # real <- analyticalIntegrals[5, 5]

        # 1. Open jpeg file
        plotName <- paste("convergenceMean", toString(whichGenz), toString(dim), "Dimensions", "NoSeqDes", ".eps", sep = "")
        plotPath <- gsub(paste("../report/Figures/", toString(whichGenz), "/", plotName, sep=""), pattern = "csv", replacement="jpg")
        postscript(plotPath, width = 2100, height = 1794)
        # jpeg(plotPath, width = 2100, height = 1794, res=200)
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
        
        par(mar=c(7, 6, 6, 10))
        # plot(x = c(1:num_iterations), y = predictionMonteCarlo$meanValueMonteCarlo,
        #     pch = 16, type = "l",
        #     xlab = "Number of epochs N", ylab = "Integral approximation", col = "blue",
        #     ylim = c(0.5, 0.77), 
        #     lty = 1,
        #     xaxs="i", yaxs="i",
        #      cex.lab = 2.2,
        #      cex.axis = 1.8
        #     )
        # lines(x = c(1:num_iterations), predictionBART$meanValueBART, type = 'l', col = "red", lty = 1)
        # lines(x = c(1:num_iterations), predictionGPBQ$meanValueGP, type = 'l', col = "green", lty = 1)
        # lines(x = c(1:num_iterations), predictionGPBQNoSeqDesign$meanValueGPNoSeq, type = 'l', col = "chocolate", lty = 1)
        # lines(x = c(1:num_iterations), results[, 2], type = 'l', col = "cadetblue", lty = 1)
        # lines(x = c(1:num_iterations), results[, 3], type = 'l', col = "gold", lty = 1)
        # lines(x = c(1:num_iterations), results[, 4], type = 'l', col = "magenta", lty = 1)
        # legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ", "Actual", "GP 10", "GP 30", "GP 50", "GP 70"),
        #     col=c("blue", "red", "green", "chocolate", "black", "cadetblue", "gold", "magenta"), cex=0.8, lty = c(1,1,1,1))
        # axis(4, at = signif(real, digits=4), las = 2, cex.axis=1.8)

        plot(x = c(1:num_iterations), y = results[, 6],
            pch = 16, type = "l",
            xlab = "Number of epochs N", ylab = "Integral estimate", col = "blue",
            ylim = c(0.142, 0.175), 
            lty = 1,
            xaxs="i", yaxs="i",
            cex.lab = 2.2,
            cex.axis = 1.8
        )
        lines(x = c(1:num_iterations), results[, 2], type = 'l', col = "chocolate", lty = 1)
        lines(x = c(1:num_iterations), results[, 3], type = 'l', col = "cadetblue", lty = 1)
        lines(x = c(1:num_iterations), results[, 4], type = 'l', col = "gold", lty = 1)
        lines(x = c(1:num_iterations), results[, 5], type = 'l', col = "magenta", lty = 1)
        abline(a = real, b = 0, lty = 4, col = "black")
        legend("topright", legend=c("GP 10", "GP 30", "GP 50", "GP 70", "GP 100", "Actual"),
            col=c("chocolate", "cadetblue", "gold", "magenta", "blue", "black"), cex=0.8, lty = c(1,1,1,1))
        axis(4, at = signif(real, digits=4), las = 2, cex.axis=1.8)

        dev.off()

    }
    cat(genzFunctionName, "done", '\n')
}


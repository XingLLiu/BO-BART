# Plot results for each Genz integral after 500 epochs.
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To see results:
#     Plots stored in {ROOT}/report/Figures for plots

# Load required packages
source("./packages/requiredPackages.R")
requiredPackages()


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

        # Source the correct file
        fileName <- paste('results', toString(whichGenz), 'dim', toString(dim), '.csv', sep='')
        filePath <- paste('./genz/results', toString(whichGenz), fileName, sep='/')

        # Retrieve integral values
        integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
        predictionBART <- data.frame("BARTMean" = integrals[, 2], "standardDeviationBART" = integrals[, 3])
        predictionMonteCarlo <- data.frame("meanValueMonteCarlo" = integrals[, 4], "standardDeviationMonteCarlo" = integrals[, 5])
        predictionGPBQ <-  data.frame("meanValueGP" = integrals[, 6], "standardDeviationGP" = integrals[, 7])

        # Retrieve analytical integrals
        whichDimension <- which(dim == dimensionsList)
        analyticalIntegrals <- read.csv("./genz/integrals.csv", header = FALSE)
        real <- analyticalIntegrals[whichGenz, whichDimension]

        # 1. Open jpeg file
        plotPath <- gsub(paste("../report/Figures/", toString(whichGenz), "/", fileName, sep=""), pattern = "csv", replacement="jpg")
        jpeg(plotPath, width = 2100, height = 1794, res=200)
        # 2. Create the plot
        par(mfrow = c(1,2), pty = "s")
        plot(x = c(1:num_iterations), y = predictionMonteCarlo$meanValueMonteCarlo,
            pch = 16, type = "l",
            xlab = "Number of epochs N", ylab = "Integral approximation", col = "blue",
            main = paste("Convergence of methods: mean vs N \nusing", genzFunctionName, "with", num_iterations, "epochs in", dim, "dim"),
            ylim = c(-real, real + real), 
            lty = 1,
            xaxs="i", yaxs="i"
            )
        lines(x = c(1:num_iterations), predictionBART$meanValueBART, type = 'l', col = "red", lty = 1)
        lines(x = c(1:num_iterations), predictionGPBQ$meanValueGP, type = 'l', col = "green", lty = 1)
        abline(a = real, b = 0, lty = 4)
        legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ", "Actual"),
            col=c("blue", "red", "green", "black"), cex=0.8, lty = c(1,1,1,1))

        # 2. Create the plot
        plot(x = log(c(2:num_iterations)), y = log(predictionMonteCarlo$standardDeviationMonteCarlo[-1]),
            pch = 16, type = "l",
            xlab = "Number of epochs N", ylab = "Log standard deviation", col = "blue",
            main = paste("Convergence of methods: log(sd) vs log(N) \nusing", genzFunctionName, "with", num_iterations, "epochs in", dim, "dim"),
            lty = 1,
            xaxs="i", yaxs="i")
        lines(x = log(c(2:num_iterations)), log(predictionBART$standardDeviationBART[-1]), type = 'l', col = "red", lty = 1)
        lines(x = log(c(2:num_iterations)), log(predictionGPBQ$standardDeviationGP[-1]), type = 'l', col = "green", lty = 1)
        legend("topleft", legend=c("MC Integration", "BART BQ", "GP BQ"),
            col=c("blue", "red", "green"), cex=0.8, lty = c(1,1,1,1))
        # 3. Close the file
        dev.off()

    }
    cat(genzFunctionName, "done", '\n')
}
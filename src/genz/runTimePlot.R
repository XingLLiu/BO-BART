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

# Initialize BART and GP run time
BARTRunTimeAll <- matrix(NA, nrow = 6, ncol = 6)
GPRunTimeAll <- matrix(NA, nrow = 6, ncol = 6)
GPRunTimeRearranged <- matrix(NA, nrow = 6, ncol = 6)
BARTRunTime <- matrix(NA, nrow = 1, ncol = 6)

# Retrieve GP run time
GPRunTimeAll <- read.csv("../results/genz/GPRunTime/GPTime.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

for (j in 1:6){
    for (i in 1:6){

        # global parameters: dimension
        dimensionsList <- c(1,2,3,5,10,20)
        args <- cbind(dimensionsList, c(1:6))
        dim <- args[j, 1]
        whichGenz <- args[i, 2]
        num_iterations <- 500
        
        # Skil if dim = 1 for discontinuous
        if (whichGenz == 3 & dim == 1) { BARTRunTimeAll[i, j] <- NA; next} 

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
        BARTRunTimeAll[i, j] <- integrals[1, 9]

        # Rearrange GP run time
        GPRunTimeRearranged[i, j] <- GPRunTimeAll[1, (1 + 6 * i + j)]
    }

    cat("Dim", dim, "done", '\n')
}

# Compute mean time of all genz 
GPRunTime <- colMeans(GPRunTimeRearranged, na.rm = TRUE)
BARTRunTime <- colMeans(BARTRunTimeAll, na.rm = TRUE)


# 1. Open jpeg file
# plotName <- paste("convergenceMean", toString(whichGenz), toString(dim), "Dimensions", ".eps", sep = "")
# plotPath <- gsub(paste("../report/Figures/", toString(whichGenz), "/", plotName, sep=""), pattern = "csv", replacement="jpg")

plotPath <- "../report/Figures/runTimePlot.eps"
postscript(plotPath, width = 2100, height = 1794)

# 2. Create the plot
# Set y limits

par(mar=c(7, 7, 4, 4))
plot(x = 1:6, y = BARTRunTime,
    pch = 16, type = "b",
    xlab = "Dimension", ylab = "Time", col = "black",
    ylim = c(0, 20000),
    xlim = c(0.8, 6.2),
    cex = 1.5,
    lty = 2,
    xaxs="i", yaxs="i",
    xaxt = "n",
    cex.lab = 2.2,
    cex.axis = 1.8
    )
lines(GPRunTime, col = "blue", type = "o", lty = 2, pch = 19, cex = 1.8)
axis(1, at=1:6, labels=c(1, 2, 3, 5, 10, 20), cex.axis=1.8)
legend("topleft", legend=c("BART BQ", "GP BQ"), col=c("black", "blue"), lty=rep(2, 2), cex=1.2, lwd=rep(2,2))

# 3. Close the file
dev.off()

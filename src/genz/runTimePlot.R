# Plot run time for each dimension
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To see results:
#     Plots stored in {ROOT}/report/Figures for plots


methods <- c("BART", "MI", "GP")
bestMethod <- matrix(NA, ncol = 3, nrow = 6)
rownames(bestMethod) <- c("cont", "copeak", "gaussian", "oscil", "prpeak", "step")
colnames(bestMethod) <- c("1", "10", "20")

# Initialize BART and GP run time
BARTRunTimeAll <- matrix(NA, nrow = 7, ncol = 3)
GPRunTimeAll <- matrix(NA, nrow = 7, ncol = 3)
GPRunTimeRearranged <- matrix(NA, nrow = 7, ncol = 3)
BARTRunTime <- matrix(NA, nrow = 1, ncol = 3)

# Retrieve GP run time
# GPRunTimeAll <- read.csv("../../results/genz/GPRunTime/GPTime.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)

for (j in c(1, 2, 3)){
    for (i in c(1, 2, 4, 5, 6, 7)){

        # global parameters: dimension
        dimensionsList <- c(1, 10,20)
        dim <- dimensionsList[j]
        whichGenz <- i
        
        # Skil if dim = 1 for discontinuous
        if (whichGenz == 3 & dim == 1) { BARTRunTimeAll[i, j] <- NA; next} 

        # Find Genz function
        if (whichGenz == 1) { genzFunctionName <-  "cont" }
        if (whichGenz == 2) { genzFunctionName <-  "copeak" }
        if (whichGenz == 3) { genzFunctionName <-  "disc" }
        if (whichGenz == 4) { genzFunctionName <-  "gaussian" }
        if (whichGenz == 5) { genzFunctionName <-  "oscil" }
        if (whichGenz == 6) { genzFunctionName <-  "prpeak" }
        if (whichGenz == 7) { genzFunctionName <-  "step" }

        for (num_cv in 1:5) {
            # Set path for estimated integral values
            fileName <- paste(toString(genzFunctionName), 'Dim', toString(dim), "", "Uniform", "_", toString(num_cv),  '.csv', sep='')
            filePath <- paste('../../mlbox/results_3_matern_k2/results/genz', toString(whichGenz), fileName, sep='/')
    
            # Retrieve BART run time
            integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
            BARTRunTimeAll[i, j] <- integrals[1, 9]
    
            # Rearrange GP run time
            GPRunTimeAll[i, j] <- integrals[1, 11]
        }
    }
    cat("Dim", dim, "done", '\n')
}
for (j in c(1, 2, 4, 5, 6, 7)) {
    pdf(paste("../../Figures_matern32_k2", "/", toString(j), "runtime",  ".pdf", sep = ""), width = 5, height = 5.5)
    par(pty = "s")
    plot(c(1, 10, 20), 
         BARTRunTimeAll[j, ], 
         ty = "l", 
         col = "orangered", 
         ylab = "Time in seconds",
         xlab = "dimension",
         ylim = c(0, max(BARTRunTimeAll[j, ], GPRunTimeAll[j, ])),
         cex.lab = 1.5,
         cex.axis = 1.5
    )
    points(c(1, 10, 20), BARTRunTimeAll[j, ], col = "orangered", bg='orangered', pch=21, lwd=3)
    points(c(1, 10, 20), GPRunTimeAll[j, ], col = "dodgerblue", ty = "l")
    points(c(1, 10, 20), GPRunTimeAll[j, ], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
    dev.off()
}


# # 1. Open jpeg file
# # plotName <- paste("convergenceMean", toString(whichGenz), toString(dim), "Dimensions", ".eps", sep = "")
# # plotPath <- gsub(paste("../report/Figures/", toString(whichGenz), "/", plotName, sep=""), pattern = "csv", replacement="jpg")
# 
# plotPath <- "../report/Figures/runTimePlot.eps"
# postscript(plotPath, width = 2100, height = 1794)
# 
# # 2. Create the plot
# # Set y limits
# 
# par(mar=c(7, 7, 4, 4))
# plot(x = 1:6, y = BARTRunTime,
#     pch = 16, type = "b",
#     xlab = "Dimension", ylab = "Time", col = "black",
#     ylim = c(0, 20000),
#     xlim = c(0.8, 6.2),
#     cex = 1.5,
#     lty = 2,
#     xaxs="i", yaxs="i",
#     xaxt = "n",
#     cex.lab = 2.2,
#     cex.axis = 1.8
#     )
# lines(GPRunTime, col = "blue", type = "o", lty = 2, pch = 19, cex = 1.8)
# axis(1, at=1:6, labels=c(1, 2, 3, 5, 10, 20), cex.axis=1.8)
# legend("topleft", legend=c("BART BQ", "GP BQ"), col=c("black", "blue"), lty=rep(2, 2), cex=1.2, lwd=rep(2,2))
# 
# # 3. Close the file
# dev.off()

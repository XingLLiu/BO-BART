# Set path for estimated integral values
fileName <- paste(toString("step"), 'Dim', toString(dim), "Gaussian",'.csv', sep='')
filePath <- paste('../../mlbox/results_2/results/genz', toString(whichGenz), fileName, sep='/')

# Retrieve estimated integral values
integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
predictionBART <- data.frame("meanValueBART" = integrals[, 2], "standardDeviationBART" = integrals[, 3])
predictionMonteCarlo <- data.frame("meanValueMonteCarlo" = integrals[, 4], "standardDeviationMonteCarlo" = integrals[, 5])
predictionGPBQ <-  data.frame("meanValueGP" = integrals[, 6], "standardDeviationGP" = integrals[, 7])

# Retrieve analytical integral values
whichDimension <- which(dim == dimensionsList)
analyticalIntegrals <- read.csv("../../results/genz/integrals.csv", header = FALSE)
real <- analyticalIntegrals[whichGenz, whichDimension]

# 1. Open eps file
# plotName <- paste("convergenceMean", toString(whichGenz), toString(dim), "Dimensions", sep = "")
plotName <- "convergenceMean1DimensionsGaussian.pdf"
plotPath <- gsub(paste("../../Figures_final/", toString(whichGenz), "/", plotName, sep=""), pattern = "csv", replacement="jpg")

integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
ymin <- min(c(predictionBART$meanValueBART, predictionGPBQ$meanValueGP, predictionMonteCarlo$meanValueMonteCarlo[5:50]))
ymax <- max(c(predictionBART$meanValueBART, predictionGPBQ$meanValueGP, predictionMonteCarlo$meanValueMonteCarlo[5:50]))

# posterior plots
pdf(plotPath, width = 8,5, height = 10)
plot_points <- seq(0, 50, 5)
plot(integrals$epochs, 
     integrals$MIMean, 
     ty="l", 
     main = paste("Posterior Convergence", toString(genzFunctionName), "function"),
     ylab = "Integral",
     xlab = "num_iterations",
     col = "chartreuse4",
     ylim = c(ymin, ymax)
)
points(integrals$epochs[plot_points], integrals$MIMean[plot_points], col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
points(integrals$epochs, integrals$GPMean, ty="l", col = "dodgerblue")
points(integrals$epochs[plot_points], integrals$GPMean[plot_points], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
polygon(c(integrals$epochs, rev(integrals$epochs)), 
        c(
          integrals$GPMean + 2*integrals$GPsd, 
          rev(integrals$GPMean - 2*integrals$GPsd)
        ), 
        col = adjustcolor("dodgerblue", alpha.f = 0.10), 
        border = "dodgerblue", lty = c("dashed", "solid"))
points(integrals$epochs, integrals$BARTMean, ty="l", col = "orangered")
points(integrals$epochs[plot_points], integrals$BARTMean[plot_points], col = "orangered",bg='orangered', pch=21, lwd=3)
polygon(c(integrals$epochs, rev(integrals$epochs)), 
        c(
          integrals$BARTMean + 2*integrals$BARTsd, 
          rev(integrals$BARTMean - 2*integrals$BARTsd)
        ), 
        col = adjustcolor("orangered", alpha.f = 0.10), 
        border = "orangered", lty = c("dashed", "solid"))
abline(h=integrals$actual)
legend("topright", legend=c("MI", "BART-Int", "GP-BQ", "Actual"),
       col=c("chartreuse4", "orangered", "dodgerblue", "black"), cex=1, lty = c(1,1,1,1))
dev.off()

# Plot results for each Genz integral after 500 epochs.
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To manually adjust y limits of the plots:
#     Uncomment ylim = ylims[1:2] and ylim = ylims[3:3] in the two plot functions, 
#     and input the four limits when running the file
# To see results:
#     Plots stored in {ROOT}/report/Figures for plots
# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}
ylims <- as.double(commandArgs(TRUE))
mape_error <- read.csv("figures_code/mapeValues.csv")
for (i in c(7)){
  for (j in c(1)){
    for (sequential in c("")){
      # global parameters: dimension
      dimensionsList <- c(1,10,20)
      dim <- dimensionsList[j]
      whichGenz <- i
      
      # Skil if dim = 1 for discontinuous
      if (whichGenz == 3 & dim == 1) { break } 
      
      # Find Genz function
      if (whichGenz == 1) { genzFunctionName <-  "cont" }
      if (whichGenz == 2) { genzFunctionName <-  "copeak" }
      if (whichGenz == 3) { genzFunctionName <-  "disc" }
      if (whichGenz == 4) { genzFunctionName <-  "gaussian" }
      if (whichGenz == 5) { genzFunctionName <-  "oscil" }
      if (whichGenz == 6) { genzFunctionName <-  "prpeak" }
      if (whichGenz == 7) { genzFunctionName <-  "step" }
      if (whichGenz == 8) { genzFunctionName <-  "additive_gaussian" }
      
      for (num_cv in c(14)) {
        # Set path for estimated integral values
        fileName <- paste(toString(genzFunctionName), 'Dim', toString(dim), "Uniform_", toString(num_cv),'.csv', sep='')
        filePath <- paste('results/newest_genz_runs', toString(whichGenz), fileName, sep='/')
        
        # Retrieve estimated integral values
        integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
        predictionBART <- data.frame("meanValueBART" = integrals[, 2], "standardDeviationBART" = integrals[, 3])
        predictionMonteCarlo <- data.frame("meanValueMonteCarlo" = integrals[, 4], "standardDeviationMonteCarlo" = integrals[, 5])
        predictionGPBQ <-  data.frame("meanValueGP" = integrals[, 6], "standardDeviationGP" = integrals[, 7])
        
        # Retrieve analytical integral values
        whichDimension <- which(dim == dimensionsList)
        analyticalIntegrals <- read.csv("results/genz/integrals.csv", header = FALSE)
        if (whichGenz==7) {
          real = 0.5
        } else {
          real <- analyticalIntegrals[whichGenz, whichDimension]
        }
        
        
        # 1. Open eps file
        # plotName <- paste("convergenceMean", toString(whichGenz), toString(dim), "Dimensions", sep = "")
        plotName <- paste("convergenceMean", toString(whichGenz), sequential, toString(dim), "Dimensions_", toString(num_cv),".pdf", sep = "")
        plotPath <- gsub(paste("Figures/", toString(whichGenz), "/", plotName, sep=""), pattern = "csv", replacement="jpg")
        plotName_error <- paste("convergenceMean", toString(whichGenz), sequential, toString(dim), "Dimensions_error_", toString(num_cv),".pdf", sep = "")
        plotPath_error <- gsub(paste("Figures/", toString(whichGenz), "/", plotName_error, sep=""), pattern = "csv", replacement="jpg")
        
        integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
        
        ymin <- min(c(predictionBART$meanValueBART - 2*predictionBART$standardDeviationBART, 
                      predictionGPBQ$meanValueGP - 2*predictionGPBQ$standardDeviationGP,
                      predictionMonteCarlo$meanValueMonteCarlo[1:20]))
        ymax <- max(c(predictionBART$meanValueBART, 
                      predictionGPBQ$meanValueGP + 2*predictionGPBQ$standardDeviationGP,
                      predictionMonteCarlo$meanValueMonteCarlo[1:20]))
        
        # posterior plots
        pdf(plotPath, width = 8.5, height = 8)
        par(pty = "s")
        plot_points <- seq(0, 50, 5)
        plot(integrals$epochs, 
             integrals$MIMean, 
             ty="l", 
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
        # 
        # print(plotPath)
        # pdf(plotPath, width = 10, height = 5)
        # # 2. Create the plot
        # # Set y limits
        # yLimMean <- cbind(quantile(predictionBART$meanValueBART, probs=c(0, 1), na.rm=TRUE),
        #                quantile(predictionMonteCarlo$meanValueMonteCarlo, probs=c(0, 1), na.rm=TRUE))
        # 
        # yLowLimMean <- min(yLimMean[, 1])
        # if (yLowLimMean < 0) {
        #     scalingFactor <- 1.1
        # } else {
        #     scalingFactor <- 0.5
        # }
        # yLowLimMean <-  yLowLimMean * scalingFactor
        # 
        # yHighLimMean <- max(yLimMean[, 2])
        # if (yHighLimMean > 0) {
        #     scalingFactor <- 0.8
        # } else {
        #     scalingFactor <- 0.9
        # }
        # yHighLimMean <-  yHighLimMean * scalingFactor
        # 
        # ymin <- min(c(predictionBART$meanValueBART, predictionGPBQ$meanValueGP, predictionMonteCarlo$meanValueMonteCarlo[5:500]))
        # ymax <- max(c(predictionBART$meanValueBART, predictionGPBQ$meanValueGP, predictionMonteCarlo$meanValueMonteCarlo[5:500]))
        # 
        # 
        # par(mfrow = c(1,2), pty = "s", mar=c(2, 5, 2, 5))
        # plot(x = c(1:num_iterations), y = predictionMonteCarlo$meanValueMonteCarlo,
        #     pch = 16, type = "l",
        #     xlab = "Number of epochs N", ylab = "Integral estimates", col = "blue",
        #     ylim = c(ymin, ymax), 
        #     # ylim = ylims[1:2yHighLimMean
        #     xaxs="i", yaxs="i",
        #     cex.lab = 1.3,
        #     cex.axis = 1.3
        #     )
        # lines(x = c(1:num_iterations), predictionBART$meanValueBART, type = 'l', col = "red", lty = 1)
        # lines(x = c(1:num_iterations), predictionGPBQ$meanValueGP, type = 'l', col = "green", lty = 1)
        # abline(a = real, b = 0, lty = 4)
        # axis(4, at = signif(real, digits=4), las = 2, cex.axis = 1.2)
        # legend("topright", legend=c("MC Integration", "BART BQ", "GP BQ", "Actual"),
        #     col=c("blue", "red", "green", "black"), cex=0.7, lty = c(1,1,1,1), bg="white")
        # 
        # # 2. Create the plot
        # # Set y limits
        # yLimSd <- cbind(quantile(log(predictionBART$standardDeviationBART), probs=c(0.1, 0.9), na.rm=TRUE),
        #                quantile(log(predictionMonteCarlo$standardDeviationMonteCarlo), probs=c(0.1, 0.9), na.rm=TRUE))
        # 
        # yLowLimSd <- min(yLimSd[, 1])
        # if (yLowLimSd < 0) {
        #     scalingFactor <- 1
        # } else {
        #     scalingFactor <- -0.2
        # }
        # yLowLimSd <-  yLowLimSd * scalingFactor
        # 
        # yHighLimSd <- max(yLimSd[, 2])
        # if (yHighLimSd > 0) {
        #     scalingFactor <- 1
        # } else {
        #     scalingFactor <- -1
        # }
        # yHighLimSd <-  yHighLimSd * scalingFactor
        # yHighLimSd <- yHighLimSd * scalingFactor
        # 
        # ymin <- min(c(predictionBART$standardDeviationBART, predictionGPBQ$standardDeviationGP, predictionMonteCarlo$standardDeviationMonteCarlo[5:500]))
        # ymax <- max(c(predictionBART$standardDeviationBART, predictionGPBQ$standardDeviationGP, predictionMonteCarlo$standardDeviationMonteCarlo[5:500]))
        # 
        # plot(x = log(c(2:num_iterations)), y = log(predictionMonteCarlo$standardDeviationMonteCarlo[-1]/sqrt(c(1:(num_iterations-1)))),
        #     pch = 16, type = "l",
        #     xlab = "Log number of epochs N", ylab = "Log standard deviation", col = "blue",
        #     lty = 1,
        #     ylim = c(yLowLimSd, yHighLimSd),
        #     # ylim = ylims[3:4],
        #     xaxs="i", yaxs="i",
        #     cex.lab = 1.3,
        #     cex.axis = 1.3
        #     )
        # lines(x = log(c(2:num_iterations)), log(predictionBART$standardDeviationBART[-1]), type = 'l', col = "red", lty = 1)
        # legend("topright", legend=c("MC Integration", "BART BQ"), bg="white",
        #       col=c("blue", "red"), cex=0.7, lty = c(1,1,1,1), pt.cex = 1.5)   
        # 3. Close the file
        dev.off()
        
        pdf(plotPath_error, width = 8,5, height = 10)
        plot(integrals$epochs,
             abs(integrals$MIMean - real),
             ty="l",
             ylab = "Absolute Error",
             xlab = "num_iterations",
             col = "chartreuse4"
        )
        abline(h=0)
        points(integrals$epochs[plot_points], abs(integrals$MIMean[plot_points] - real), col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
        points(integrals$epochs, abs(integrals$GPMean-real), ty="l", col = "dodgerblue")
        points(integrals$epochs[plot_points], abs(integrals$GPMean[plot_points] - real), col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
        points(integrals$epochs, abs(integrals$BARTMean - real), ty="l", col = "orangered")
        points(integrals$epochs[plot_points], abs(integrals$BARTMean[plot_points] - real), col = "orangered",bg='orangered', pch=21, lwd=3)
        dev.off()
        cat(genzFunctionName, "done", '\n')
        
        # pdf(plotPath_error, width = 8,5, height = 10)
        # plot(integrals$epochs, 
        #      log(abs(integrals$MIMean - real)), 
        #      ty="l", 
        #      ylab = "Absolute Error",
        #      xlab = "num_iterations",
        #      col = "chartreuse4"
        # )
        # abline(h=0)
        # points(integrals$epochs[plot_points], log(abs(integrals$MIMean[plot_points] - real)), col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
        # points(integrals$epochs, log(abs(integrals$GPMean-real)), ty="l", col = "dodgerblue")
        # points(integrals$epochs[plot_points], log(abs(integrals$GPMean[plot_points] - real)), col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
        # points(integrals$epochs, log(abs(integrals$BARTMean - real)), ty="l", col = "orangered")
        # points(integrals$epochs[plot_points], log(abs(integrals$BARTMean[plot_points] - real)), col = "orangered",bg='orangered', pch=21, lwd=3)
        # dev.off()
        # cat(genzFunctionName, "done", '\n')
      }
    }
    pdf(paste("Figures", "/combined_", toString(whichGenz), ".pdf", sep = ""), width = 9, height = 5)
    par(mfrow = c(1,2), pty = "s")
    plot(integrals$epochs[abs(integrals$MIMean - real)!=0]+50,
        abs(integrals$MIMean - real)[abs(integrals$MIMean - real)!=0],
        ty="l",
        ylab = "Absolute Error",
        xlab = expression("n'" + n[seq]),
        col = "chartreuse4",
        cex.lab = 2,
        cex.axis = 2,
        log = "y",
        ylim = c(0.00001, 0.1)
    )
    abline(h=0)
    points(integrals$epochs[plot_points]+50, abs(integrals$MIMean[plot_points] - real), col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
    points(integrals$epochs+50, abs(integrals$GPMean-real), ty="l", col = "dodgerblue")
    points(integrals$epochs[plot_points]+50, abs(integrals$GPMean[plot_points] - real), col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
    points(integrals$epochs+50, abs(integrals$BARTMean - real), ty="l", col = "orangered")
    points(integrals$epochs[plot_points]+50, abs(integrals$BARTMean[plot_points] - real), col = "orangered",bg='orangered', pch=21, lwd=3)
    legend("bottomleft", legend=c("BART-Int", "GP-BQ", "MI"),
          col=c("orangered", "dodgerblue", "chartreuse4"), cex=1.6, lty = c(1,1,1,1), bty="n")
    
    plot_points <- seq(0, 50, 5)
    
    plot(integrals$epochs+50, 
         integrals$MIMean, 
         ty="l", 
         ylab = "Integral",
         xlab = expression("n'"+n[seq]),
         col = "chartreuse4",
         ylim = c(ymin, ymax),
         cex.lab = 1.7,
         cex.axis = 1.7
    )
    points(integrals$epochs[plot_points]+50, integrals$MIMean[plot_points], col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
    points(integrals$epochs+50, integrals$GPMean, ty="l", col = "dodgerblue")
    points(integrals$epochs[plot_points]+50, integrals$GPMean[plot_points], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
    polygon(c(integrals$epochs+50, rev(integrals$epochs+50)), 
            c(
              integrals$GPMean + 2*integrals$GPsd, 
              rev(integrals$GPMean - 2*integrals$GPsd)
            ), 
            col = adjustcolor("dodgerblue", alpha.f = 0.10), 
            border = "dodgerblue", lty = c("dashed", "solid"))
    points(integrals$epochs+50, integrals$BARTMean, ty="l", col = "orangered")
    points(integrals$epochs[plot_points]+50, integrals$BARTMean[plot_points], col = "orangered",bg='orangered', pch=21, lwd=3)
    polygon(c(integrals$epochs+50, rev(integrals$epochs+50)), 
            c(
              integrals$BARTMean + 2*integrals$BARTsd, 
              rev(integrals$BARTMean - 2*integrals$BARTsd)
            ), 
            col = adjustcolor("orangered", alpha.f = 0.10), 
            border = "orangered", lty = c("dashed", "solid"))
    for (n_seq in 1:20) {
      bart_posterior <- load("results/genz/posterior_BART_Dim1_step_14_%s.RData" %--% c(n_seq))
      if (n_seq == 1) {
        num_posterior_samples <- length(posterior_samples$posterior_samples)
        points(rep(n_seq+50, num_posterior_samples), 
             posterior_samples$posterior_samples,
             col = "orangered",bg='orangered',
             cex=0.1)
      } else {
        points(rep(n_seq+50, num_posterior_samples), posterior_samples$posterior_samples, 
               col = "orangered", bg='orangered', cex=0.1)
      }
    }
    abline(h=integrals$actual)
    dev.off()
  }
}
pdf(paste("Figures", "/combined_", toString(whichGenz), "BART.pdf", sep = ""), width = 5, height = 5)
par(pty="s")
for (n_seq in 1:20) {
  bart_posterior <- load("results/genz/posterior_BART_Dim1_step_14_%s.RData" %--% c(n_seq))
  if (n_seq == 1) {
    num_posterior_samples <- length(posterior_samples$posterior_samples)
    plot(rep(n_seq, num_posterior_samples), 
         posterior_samples$posterior_samples,
         xlim=c(1,20), 
         ylim=c(0.4, 0.6), col = "orangered",bg='orangered',
         cex=0.1)
  } else {
    points(rep(n_seq, num_posterior_samples), posterior_samples$posterior_samples, 
           col = "orangered", bg='orangered', cex=0.1)
  }
}
abline(h=0.5)
dev.off()

num_new_surveys = 1000
for (num_cv in 1:5) {
  # compute population average income estimates by BARTBQ
  
  results <- read.csv(paste("meanPopulationStudy/", "results", num_cv, ".csv", sep = ""))  
  real <- results$PoptMean[1]
  # 1. Open jpeg file
  plot_points <- seq(0, 1000, 100)
  pdf(paste("results", num_cv, ".pdf", sep=""), width = 8.5, height = 5)
  par(mfrow = c(1,2), pty = "s")
  ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$BRSMean - 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
  ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$BRSMean + 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
  plot(results$epochs,
       abs(results$BARTMean - results$PoptMean),
       ty="l",
       ylab = "Absolute Error",
       xlab = "num_iterations",
       col = "orangered",
       ylim = c(0, 40),
       cex.lab = 1.5,
       cex.axis = 1.5
  )
  abline(h=0)
  points(results$epochs[plot_points], abs(results$BARTMean[plot_points] - real), col = "orangered",bg='orangered', pch=21, lwd=3)
  points(results$epochs[plot_points], abs(results$MIMean[plot_points] - real), col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
  points(results$epochs, abs(results$MIMean - real), ty="l", col = "chartreuse4")
  points(results$epochs, abs(results$BRSMean-real), ty="l", col = "dodgerblue")
  points(results$epochs[plot_points], abs(results$BRSMean[plot_points] - real), col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
  legend("topright", legend=c("Block sampling", "BART-Int", "BR", "Actual"),
         col=c("chartreuse4", "orangered", "dodgerblue", "black"), cex=0.5, lty = c(1,1,1,1))
  
  ymin <- min(results$BARTMean - 2*results$BARTsd, results$BRSMean[-num_new_surveys] - 2*results$BRSsd[-num_new_surveys], results$MIMean)
  ymax <- max(results$BARTMean + 2*results$BARTsd, results$BRSMean[-num_new_surveys] + 2*results$BRSsd[-num_new_surveys], results$MIMean)
  
  plot(results$epochs, 
       results$MIMean, 
       ty="l", 
       ylab = "Mean Population",
       xlab = "num_iterations",
       col = "chartreuse4",
       ylim = c(ymin, ymax),
       cex.lab = 1.5,
       cex.axis = 1.5
  )
  points(results$epochs[plot_points], results$MIMean[plot_points], col = "chartreuse4", bg='chartreuse4', pch=21, lwd=3)
  points(results$epochs, results$BRSMean, ty="l", col = "dodgerblue")
  points(results$epochs[plot_points], results$BRSMean[plot_points], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
  polygon(c(results$epochs, rev(results$epochs)), 
          c(
            results$BRSMean + 2*results$BRSsd, 
            rev(results$BRSMean - 2*results$BRSsd)
          ), 
          col = adjustcolor("dodgerblue", alpha.f = 0.10), 
          border = "dodgerblue", lty = c("dashed", "solid"))
  points(results$epochs, results$BARTMean, ty="l", col = "orangered")
  points(results$epochs[plot_points], results$BARTMean[plot_points], col = "orangered",bg='orangered', pch=21, lwd=3)
  polygon(c(results$epochs, rev(results$epochs)), 
          c(
            results$BARTMean + 2*results$BARTsd, 
            rev(results$BARTMean - 2*results$BARTsd)
          ), 
          col = adjustcolor("orangered", alpha.f = 0.10), 
          border = "orangered", lty = c("dashed", "solid"))
  abline(h=results$PoptMean)
  dev.off()
}

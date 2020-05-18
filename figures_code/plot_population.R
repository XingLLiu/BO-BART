
# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

# c(c(1:8), c(10), c(13:50))
for (num_cv in c(c(1:8), c(10), c(13:50))) {
  results <- read.csv("results/populationStudy/results%s.csv" %--% num_cv)
  # gpresults <- read.csv("results/populationStudy/gpresults%s.csv" %--% num_cv)
  # results$GPMean <- gpresults$GPMean
  # results$GPsd <- gpresults$GPsd
  real <- results$BpoptMean[1]
  real <- mean(c(trainY, candidateY))
  # 1. Open jpeg file
  pdf(paste0(plotPath, "results", num_cv, ".pdf"), width = 8, height = 10)
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
  points(results$epochs, abs(results$BRSMean - real), ty="l", col = "dodgerblue")
  points(results$epochs, abs(results$GPMean - real), ty="l", col = "darkgoldenrod")
  
  legend("topright", legend=c("Block sampling", "BART-Int", "Monte Carlo sampling", "GPBQ"),
         col=c("dodgerblue", "orangered", "chartreuse4", "darkgoldenrod"), cex=0.8, lty = c(1,1,1,1))
  
  # ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$BRSMean - 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
  # ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$BRSMean + 2*results$BRSsd, results$MIMean[1:num_new_surveys]))
  ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$GPMean - 2*results$GPsd))
  ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$GPMean + 2*results$GPsd))    
  
  plot(results$epochs, 
       results$MIMean, 
       ty="l", 
       ylab = "Mean Population",
       xlab = "num_iterations",
       col = "chartreuse4",
       ylim = c(ymin, ymax)
  )
  polygon(c(results$epochs, rev(results$epochs)), 
          c(
            results$MIMean + 2*results$MIsd, 
            rev(results$MIMean - 2*results$MIsd)
          ), 
          col = adjustcolor("chartreuse4", alpha.f = 0.10), 
          border = "chartreuse4", lty = c("dashed", "solid"))
  points(results$epochs, results$BRSMean, ty="l", col = "dodgerblue")
  polygon(c(results$epochs, rev(results$epochs)), 
          c(
            results$BRSMean + 2*results$BRSsd, 
            rev(results$BRSMean - 2*results$BRSsd)
          ), 
          col = adjustcolor("dodgerblue", alpha.f = 0.10), 
          border = "dodgerblue", lty = c("dashed", "solid"))
  points(results$epochs, results$BARTMean, ty="l", col = "orangered")
  polygon(c(results$epochs, rev(results$epochs)), 
          c(
            results$BARTMean + 2*results$BARTsd, 
            rev(results$BARTMean - 2*results$BARTsd)
          ), 
          col = adjustcolor("orangered", alpha.f = 0.10), 
          border = "orangered", lty = c("dashed", "solid"))
  
  points(results$epochs, results$GPMean, ty="l", col = "darkgoldenrod")
  polygon(c(results$epochs, rev(results$epochs)), 
          c(
            results$GPMean + 2*results$GPsd, 
            rev(results$GPMean - 2*results$GPsd)
          ), 
          col = adjustcolor("darkgoldenrod", alpha.f = 0.10), 
          border = "darkgoldenrod", lty = c("dashed", "solid"))
  
  abline(h=results$PoptMean)
  dev.off()
}

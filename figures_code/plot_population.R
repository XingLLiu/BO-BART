
# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}
# paths to save results and plots
plotPath <- "Figures/populationStudy/"
num_design=500
bart = c()
gp = c()
mi = c()
    
# c(c(1:8), c(10), c(13:50))
for (num_cv in c(11)) {
  results <- read.csv("results/populationStudy/populationStudy/results%s.csv" %--% num_cv)
  gpresults <- read.csv("hpc/gpresults%s.csv" %--% num_cv)
  # results$GPMean <- gpresults$GPMean
  # results$GPsd <- gpresults$GPsd
  real <- results$BpoptMean[1]
  bart <- c(bart, abs((results$BARTMean[1500] - real)/ real))
  gp <- c(gp, abs((gpresults$GPMean[1500] - real) / real))
  mi <- c(mi, abs((results$MIMean[1500] - real) / real))
    
  # 1. Open jpeg file
  pdf(paste0(plotPath, "results", num_cv, ".pdf"), width = 9, height = 5)
  par(mfrow = c(1,2), pty = "s")
  ymax <- max(c(abs(results$BARTMean - real), abs(results$GPMean - real)))
  plot(results$epochs+num_design,
       abs(results$BARTMean - results$PoptMean),
       ty="l",
       xlab = expression(n[seq]),
       ylab = "Absolute Error",
       col = "orangered",
       ylim = c(0.001, 1),
       log="y",
       cex.lab=1.6,
       cex.axis=2,
  )
  
  points(results$epochs+num_design, abs(results$MIMean - real), ty="l", col = "chartreuse4")
  points(results$epochs+num_design, abs(gpresults$GPMean - real), ty="l", col = "dodgerblue")
  
  ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$GPMean - 2*sqrt(gpresults$GPsd)))
  ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$GPMean + 2*sqrt(gpresults$GPsd)))    
  
  plot(results$epochs+num_design, 
        results$MIMean, 
        ty="l", 
        ylab = "Mean Log Income",
        xlab = expression(n[seq]),
        col = "chartreuse4",
        ylim = c(real-0.4, real+0.4),
        cex.lab=1.5,
        cex.axis=2,
  )

  points(results$epochs+num_design, results$BARTMean, ty="l", col = "orangered")
  polygon(c(results$epochs+num_design, rev(results$epochs+num_design)), 
          c(
            results$BARTMean + 2*sqrt(results$BARTsd), 
            rev(results$BARTMean - 2*sqrt(results$BARTsd))
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
  points(gpresults$epochs+num_design, gpresults$GPMean, ty="l", col = "dodgerblue")
  polygon(c(gpresults$epochs+num_design, rev(gpresults$epochs+num_design)), 
          c(
            gpresults$GPMean + 2*sqrt(gpresults$GPsd), 
            rev(gpresults$GPMean - 2*sqrt(gpresults$GPsd))
          ), 
          col = adjustcolor("dodgerblue", alpha.f = 0.10), 
          border = "dodgerblue", lty = c("dashed", "solid"))
  legend("bottomleft", legend=c("Monte Carlo", "BART-Int", "GPBQ"),
          col=c("chartreuse4", "orangered", "dodgerblue"), cex=1.6, lty = c(1,1,1), bty="n")
  
  abline(h=real)
  dev.off()
}
bart_mape = signif(mean(bart), 7)
bart_sd = signif(sd(bart) / sqrt(20), 7)
gp_mape = signif(mean(abs(gp)), 7)
gp_sd = signif(sd(gp) / sqrt(20), 7)
mi_mape = signif(mean(abs(mi)), 7)
mi_sd = signif(sd(mi)/sqrt(20), 7)

jpeg(paste0(plotPath, "BARTresultsFull", num_cv, ".pdf"), width = 9, height = 5)
par(pty = "s")
for (num_cv in c(1:20)) {
    results <- read.csv("results/populationStudy/populationStudy/results%s.csv" %--% num_cv)
    # gpresults <- read.csv("results/populationStudy/gpresults%s.csv" %--% num_cv)
      # results$GPMean <- gpresults$GPMean
      # results$GPsd <- gpresults$GPsd
      real <- c(results$BpoptMean[1], real)
      # 1. Open jpeg file
        # ymax <- max(c(abs(results$BARTMean - real), abs(results$BRSMean - real), abs(results$MIMean - real)))
        if (num_cv == 1) {
            ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$GPMean - 2*results$GPsd))
            ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$GPMean + 2*results$GPsd))    
            plot(results$epochs+num_design, 
                           results$BARTMean, 
                           ty="l", 
                           ylab = "Mean Log Income",
                           xlab = expression(n[seq]),
                           col = "orangered",
                           ylim = c(real[1]-0.4, real[1]+0.4),
                           cex.lab=1.5,
                           cex.axis=2,
                           lwd=0.1
                      )
          } else {
              points(results$epochs+num_design, results$BARTMean, ty="l", col = "orangered", lwd =0.1)
            }
    }
abline(h=mean(real))
dev.off()

jpeg(paste0(plotPath, "GPresultsFull", num_cv, ".pdf"), width = 9, height = 5)
par(pty = "s")
for (num_cv in c(1:20)) {
  results <- read.csv("results/populationStudy/populationStudy/results%s.csv" %--% num_cv)
  # gpresults <- read.csv("results/populationStudy/gpresults%s.csv" %--% num_cv)
  # results$GPMean <- gpresults$GPMean
  # results$GPsd <- gpresults$GPsd
  real <- c(results$BpoptMean[1], real)
  # 1. Open jpeg file
  # ymax <- max(c(abs(results$BARTMean - real), abs(results$BRSMean - real), abs(results$MIMean - real)))
  if (num_cv == 1) {
    ymin <- min(c(results$BARTMean - 2*results$BARTsd, results$GPMean - 2*results$GPsd))
    ymax <- max(c(results$BARTMean + 2*results$BARTsd, results$GPMean + 2*results$GPsd))    
    plot(results$epochs+num_design, 
      results$GPMean, 
      ty="l", 
      ylab = "Mean Log Income",
      xlab = expression(n[seq]),
      col = "chartreuse4",
      ylim = c(real[1]-0.4, real[1]+0.4),
      cex.lab=1.5,
      cex.axis=2,
      lwd=1
    )
  } else {
      points(results$epochs+num_design, results$GPMean, ty="l", col = "chartreuse4", lwd =1)
  }
}
abline(h=mean(real))
dev.off()
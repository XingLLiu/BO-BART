bart_error <- 0
MI_error <- 0
BR_error <- 0
for (num_cv in 1:5) {
  results <- read.csv(paste("meanPopulationStudy/", "results", num_cv, ".csv", sep = ""))
  epochs <- nrow(results)
  
  bart_error <- bart_error + abs(results$BARTMean[epochs] - results$PoptMean[epochs])
  MI_error <- MI_error + abs(results$MIMean[epochs] - results$PoptMean[epochs])
  BR_error <- BR_error + abs(results$BRSMean[epochs] - results$PoptMean[epochs])
}

cat("BART ", bart_error/5, "\n")
cat("MI ", MI_error/5,"\n")
cat("BR ", BR_error/5,"\n")
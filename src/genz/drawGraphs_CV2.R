# Plot results for each Genz integral after 500 epochs.
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To manually adjust y limits of the plots:
#     Uncomment ylim = ylims[1:2] and ylim = ylims[3:3] in the two plot functions, 
#     and input the four limits when running the file
# To see results:
#     Plots stored in {ROOT}/report/Figures for plots


# ylims <- as.double(commandArgs(TRUE))
# # abs_error <- read.csv("../../Figures_matern32_k2/rmseValues.csv")
# dimensionsList <- c(1,2,3,5,10,20)
# measure <- "Uniform"   # or "Gaussian"
# num_cv_total <- 5

# for (whichGenz in c(9)){
#   for (dim in c(10)){
#     for (sequential in c("")){

#       # Skip if dim = 1 for discontinuous
#       if (whichGenz == 3 & dim == 1) { break } 
      
#       # Find Genz function
#       if (whichGenz == 1) { genzFunctionName <-  "cont" }
#       if (whichGenz == 2) { genzFunctionName <-  "copeak" }
#       if (whichGenz == 3) { genzFunctionName <-  "disc" }
#       if (whichGenz == 4) { genzFunctionName <-  "gaussian" }
#       if (whichGenz == 5) { genzFunctionName <-  "oscil" }
#       if (whichGenz == 6) { genzFunctionName <-  "prpeak" }
#       if (whichGenz == 7) { genzFunctionName <-  "step" }
#       if (whichGenz == 9) { genzFunctionName <-  "additive_gaussian" }
      
#       for (num_cv in 1:num_cv_total) {
#         # Set path for estimated integral values
#         fileName <- paste(toString(genzFunctionName), 'Dim', toString(dim), sequential, measure, "_", toString(num_cv),'.csv', sep='')
#         filePath <- paste('../../results/genz', toString(whichGenz), fileName, sep='/')
        
#         # Retrieve estimated integral values
#         integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
#         predictionBART <- data.frame("meanValueBART" = integrals[, 2], "standardDeviationBART" = integrals[, 3])
#         predictionMonteCarlo <- data.frame("meanValueMonteCarlo" = integrals[, 4], "standardDeviationMonteCarlo" = integrals[, 5])
#         predictionGPBQ <-  data.frame("meanValueGP" = integrals[, 6], "standardDeviationGP" = integrals[, 7])
        
#         if (num_cv == 1){
#           predictionComb <- data.frame("meanValueBART" = predictionBART$meanValueBART,
#                                        "standardDeviationBART" = predictionBART$standardDeviationBART,
#                                        "meanValueMonteCarlo" = predictionMonteCarlo$meanValueMonteCarlo,
#                                        "standardDeviationMonteCarlo" = predictionMonteCarlo$standardDeviationMonteCarlo,
#                                        "meanValueGP" = predictionGPBQ$meanValueGP,
#                                        "standardDeviationGP" = predictionGPBQ$standardDeviationGP)
#         } else{
#           predictionComb$meanValueBART <- cbind(predictionComb$meanValueBART, predictionBART$meanValueBART)
#           predictionComb$standardDeviationBART <- cbind(predictionComb$standardDeviationBART, predictionBART$standardDeviationBART)
#           predictionComb$meanValueMonteCarlo <- cbind(predictionComb$meanValueMonteCarlo, predictionMonteCarlo$meanValueMonteCarlo)
#           predictionComb$standardDeviationMonteCarlo <- cbind(predictionComb$standardDeviationMonteCarlo, predictionMonteCarlo$standardDeviationMonteCarlo)
#           predictionComb$meanValueGP <- cbind(predictionComb$meanValueGP, predictionGPBQ$meanValueGP)
#           predictionComb$standardDeviationGP <- cbind(predictionComb$standardDeviationGP, predictionGPBQ$standardDeviationGP)
#         }

#         # Retrieve analytical integral values
#         real <- integrals$actual[1]
        
#         # 1. Open eps file
#         plotName <- paste("convergenceMean", toString(whichGenz), sequential, toString(dim), "Dimensions_", toString(num_cv),".pdf", sep = "")
#         plotRootPath <- paste("../../Figures/", toString(whichGenz), "/", plotName, sep="")
#         plotPath <- gsub(plotRootPath, pattern = "csv", replacement="jpg")
#         plotName_error <- paste("convergenceMean", toString(whichGenz), sequential, toString(dim), "Dimensions_error_", toString(num_cv),".pdf", sep = "")
#         plotPath_error <- gsub(plotRootPath, pattern = "csv", replacement="jpg")
        
#         ymin <- min( apply(cbind(predictionBART$meanValueBART, predictionMonteCarlo$meanValueMonteCarlo, predictionGPBQ$meanValueGP), 2, FUN = min) )
#         ymax <- max( apply(cbind(predictionBART$meanValueBART, predictionMonteCarlo$meanValueMonteCarlo, predictionGPBQ$meanValueGP), 2, FUN = max) )
        
#         # Posterior plots
#         pdf(plotPath, width = 8.5, height = 8)
#         par(pty = "s")

#         plot(integrals$epochs, 
#             integrals$MIMean, 
#             ty="l", 
#             ylab = "Integral",
#             xlab = "num_iterations",
#             col = "chartreuse4",
#             ylim = c(ymin, ymax),
#             lwd = 1.5
#         )
#         polygon(c(integrals$epochs, rev(integrals$epochs)), 
#                 c(
#                   integrals$GPMean + 2*integrals$GPsd, 
#                   rev(integrals$GPMean - 2*integrals$GPsd)
#                 ), 
#                 col = adjustcolor("dodgerblue", alpha.f = 0.2), 
#                 border = adjustcolor("dodgerblue", alpha.f = 0.2),
#                 lty = "solid")
#         polygon(c(integrals$epochs, rev(integrals$epochs)), 
#                 c(
#                   integrals$BARTMean + 2*integrals$BARTsd, 
#                   rev(integrals$BARTMean - 2*integrals$BARTsd)
#                 ), 
#                 col = adjustcolor("orangered", alpha.f = 0.2), 
#                 border = adjustcolor("orangered", alpha.f = 0.2),
#                 lty = "solid")
#         points(integrals$epochs, integrals$GPMean, ty="l", col = "dodgerblue", lwd = 1.5)
#         points(integrals$epochs, integrals$BARTMean, ty="l", col = "orangered", lwd = 1.5)
#         abline(h=integrals$actual)
#         legend("topright", legend=c("MI", "BART-Int", "GP-BQ", "Actual"),
#               col=c("chartreuse4", "orangered", "dodgerblue", "black"), cex=1, lty = c(1,1,1,1))
   
#         # Close the file
#         dev.off()
        
#         pdf(plotPath_error, width = 8,5, height = 10)
#         plot(integrals$epochs,
#             abs(integrals$MIMean - real),
#             ty="l",
#             ylab = "Absolute Error",
#             xlab = "num_iterations",
#             col = "chartreuse4",
#             lwd = 1.5
#         )
#         abline(h=0)
#         points(integrals$epochs, abs(integrals$GPMean-real), ty="l", col = "dodgerblue", lwd = 1.5)
#         points(integrals$epochs, abs(integrals$BARTMean - real), ty="l", col = "orangered", lwd = 1.5)
#         dev.off()

#         cat(genzFunctionName, "done", '\n')
#       }
#     }

#     # Combined plots for all cv runs
#     combinedGPsd <- apply(predictionComb$meanValueGP, 1, sd) / sqrt(num_cv_total)
#     combinedBARTsd <- apply(predictionComb$meanValueBART, 1, sd) / sqrt(num_cv_total)
#     combinedMeanMC <- apply(predictionComb$meanValueMonteCarlo, 1, mean)
#     combinedMeanGP <- apply(predictionComb$meanValueGP, 1, mean)
#     combinedMeanBART <- apply(predictionComb$meanValueBART, 1, mean)

#     pdf(paste("../../Figures", "/combined_", toString(whichGenz), ".pdf", sep = ""), width = 8.5, height = 5)
#     par(mfrow = c(1,2), pty = "s")
#     plot(integrals$epochs,
#          abs(combinedMeanMC - real),
#          ty="l",
#          ylab = "Absolute Error",
#          xlab = "num_iterations",
#          col = "chartreuse4",
#          lwd = 1.5
#     )
#     abline(h=0)
#     points(integrals$epochs, abs(combinedMeanGP - real), ty = "l", col = "dodgerblue", lwd = 1.5)
#     points(integrals$epochs, abs(combinedMeanBART - real), ty="l", col = "orangered", lwd = 1.5)
#     legend("topright", legend=c("BART-Int", "GP-BQ", "MI"),
#            col=c("orangered", "dodgerblue", "chartreuse4"), cex=1.2, lty = c(1,1,1,1))
    
#     plot(integrals$epochs, 
#          combinedMeanMC, 
#          ty="l", 
#          ylab = "Integral",
#          xlab = "num_iterations",
#          col = "chartreuse4",
#          ylim = c(ymin, ymax),
#          lwd = 1.5
#     )
#     polygon(c(integrals$epochs, rev(integrals$epochs)), 
#             c(
#               combinedMeanGP + 2*combinedGPsd, 
#               rev(combinedMeanGP - 2*combinedGPsd)
#             ), 
#             col = adjustcolor("dodgerblue", alpha.f = 0.20), 
#             border = adjustcolor("dodgerblue", alpha.f = 0.20),
#             lty = "solid")
#     polygon(c(integrals$epochs, rev(integrals$epochs)), 
#             c(
#               combinedMeanBART + 2*combinedBARTsd,
#               rev(combinedMeanBART - 2*combinedBARTsd)
#             ), 
#             col = adjustcolor("orangered", alpha.f = 0.20), 
#             border = adjustcolor("orangered", alpha.f = 0.20),
#             lty = "solid")
#     points(integrals$epochs, combinedMeanGP, ty="l", col = "dodgerblue", lwd = 1.5)
#     points(integrals$epochs, combinedMeanBART, ty="l", col = "orangered", lwd = 1.5)
#     abline(h=integrals$actual)
#     dev.off()
#   }
  
# }




plot_results <- function(args)
{  
  dims_list <- args$dims_list
  genz_list <- args$genz_list
  sequential_list <- args$sequential_list
  measure <- args$measure
  num_cv_total <- args$num_cv

  for (whichGenz in genz_list) {
    for (dim in dims_list) {
      for (sequential in sequential_list) {

        # Skip if dim = 1 for discontinuous
        if (whichGenz == 3 & dim == 1) { break } 
        
        # Find Genz function
        if (whichGenz == 1) { genzFunctionName <-  "cont" }
        if (whichGenz == 2) { genzFunctionName <-  "copeak" }
        if (whichGenz == 3) { genzFunctionName <-  "disc" }
        if (whichGenz == 4) { genzFunctionName <-  "gaussian" }
        if (whichGenz == 5) { genzFunctionName <-  "oscil" }
        if (whichGenz == 6) { genzFunctionName <-  "prpeak" }
        if (whichGenz == 7) { genzFunctionName <-  "step" }
        if (whichGenz == 9) { genzFunctionName <-  "additive_gaussian" }
        
        for (num_cv in 1:num_cv_total) {
          # Set path for estimated integral values
          fileName <- paste(toString(genzFunctionName), 'Dim', toString(dim), sequential, tools::toTitleCase(measure), "_", toString(num_cv),'.csv', sep='')
          filePath <- paste('results/genz', toString(whichGenz), fileName, sep='/')
          
          # Retrieve estimated integral values
          integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
          integrals$GPsd <- replace(integrals$GPsd, is.na(integrals$GPsd), 0)

          predictionBART <- data.frame("meanValueBART" = integrals$BARTMean, "standardDeviationBART" = integrals$BARTsd)
          predictionMonteCarlo <- data.frame("meanValueMonteCarlo" = integrals$MIMean, "standardDeviationMonteCarlo" = integrals$MIsd)
          predictionGPBQ <-  data.frame("meanValueGP" = integrals$GPMean, "standardDeviationGP" = integrals$GPsd)
          
          if (num_cv == 1){
            predictionComb <- data.frame("meanValueBART" = predictionBART$meanValueBART,
                                        "standardDeviationBART" = predictionBART$standardDeviationBART,
                                        "meanValueMonteCarlo" = predictionMonteCarlo$meanValueMonteCarlo,
                                        "standardDeviationMonteCarlo" = predictionMonteCarlo$standardDeviationMonteCarlo,
                                        "meanValueGP" = predictionGPBQ$meanValueGP,
                                        "standardDeviationGP" = predictionGPBQ$standardDeviationGP,
                                        "timeBART" = integrals$runtimeBART,
                                        "timeMC" = integrals$runtimeMI,
                                        "timeGP" = integrals$runtimeGP)

          } else{
            predictionComb$meanValueBART <- cbind(predictionComb$meanValueBART, predictionBART$meanValueBART)
            predictionComb$standardDeviationBART <- cbind(predictionComb$standardDeviationBART, predictionBART$standardDeviationBART)
            predictionComb$meanValueMonteCarlo <- cbind(predictionComb$meanValueMonteCarlo, predictionMonteCarlo$meanValueMonteCarlo)
            predictionComb$standardDeviationMonteCarlo <- cbind(predictionComb$standardDeviationMonteCarlo, predictionMonteCarlo$standardDeviationMonteCarlo)
            predictionComb$meanValueGP <- cbind(predictionComb$meanValueGP, predictionGPBQ$meanValueGP)
            predictionComb$standardDeviationGP <- cbind(predictionComb$standardDeviationGP, predictionGPBQ$standardDeviationGP)
            predictionComb$timeBART <- cbind(predictionComb$timeBART, integrals$runtimeBART)
            predictionComb$timeGP <- cbind(predictionComb$timeGP, integrals$runtimeGP)
          }

          # Retrieve analytical integral values
          real <- integrals$actual[1]
          
          # 1. Open eps file
          print(paste0("dim: ", dim))
          plotName <- paste("convergenceMean", toString(whichGenz), sequential, toString(dim), "Dimensions_", toString(num_cv),".pdf", sep = "")
          plotRootPath <- paste("Figures/", toString(whichGenz), "/", sep="")
          plotPath <- gsub(paste(plotRootPath, plotName, sep=""), pattern = "csv", replacement="jpg")
          plotName_error <- paste("convergenceMean", toString(whichGenz), sequential, toString(dim), "Dimensions_error_", toString(num_cv),".pdf", sep = "")
          plotPath_error <- gsub(paste(plotRootPath, plotName_error, sep=""), pattern = "csv", replacement="jpg")
          
          ymin <- min( apply(cbind(predictionBART$meanValueBART, predictionMonteCarlo$meanValueMonteCarlo, predictionGPBQ$meanValueGP), 2, FUN = min) )
          ymax <- max( apply(cbind(predictionBART$meanValueBART, predictionMonteCarlo$meanValueMonteCarlo, predictionGPBQ$meanValueGP), 2, FUN = max) )
          
          # Posterior plots
          pdf(plotPath, width = 8.5, height = 8)
          par(pty = "s")

          plot(integrals$epochs, 
              integrals$MIMean, 
              ty="l", 
              ylab = "Integral",
              xlab = "Number of training points",
              col = "chartreuse4",
              ylim = c(ymin, ymax),
              lwd = 1.5
          )
          polygon(c(integrals$epochs, rev(integrals$epochs)), 
                  c(
                    integrals$GPMean + 2*integrals$GPsd, 
                    rev(integrals$GPMean - 2*integrals$GPsd)
                  ), 
                  col = adjustcolor("dodgerblue", alpha.f = 0.2), 
                  border = adjustcolor("dodgerblue", alpha.f = 0.2),
                  lty = "solid")
          polygon(c(integrals$epochs, rev(integrals$epochs)), 
                  c(
                    integrals$BARTMean + 2*integrals$BARTsd, 
                    rev(integrals$BARTMean - 2*integrals$BARTsd)
                  ), 
                  col = adjustcolor("orangered", alpha.f = 0.2), 
                  border = adjustcolor("orangered", alpha.f = 0.2),
                  lty = "solid")
          points(integrals$epochs, integrals$GPMean, ty="l", col = "dodgerblue", lwd = 1.5)
          points(integrals$epochs, integrals$BARTMean, ty="l", col = "orangered", lwd = 1.5)
          abline(h=integrals$actual)
          legend("topright", legend=c("MI", "BART-Int", "GP-BQ", "Actual"),
                col=c("chartreuse4", "orangered", "dodgerblue", "black"), cex=1, lty = c(1,1,1,1))
    
          # Close the file
          dev.off()
          
          pdf(plotPath_error, width = 8,5, height = 10)
          plot(integrals$epochs,
              abs(integrals$MIMean - real),
              ty="l",
              ylab = "Absolute Error",
              xlab = "Number of training points",
              col = "chartreuse4",
              lwd = 1.5
          )
          abline(h=0)
          points(integrals$epochs, abs(integrals$GPMean-real), ty="l", col = "dodgerblue", lwd = 1.5)
          points(integrals$epochs, abs(integrals$BARTMean - real), ty="l", col = "orangered", lwd = 1.5)
          dev.off()

          cat(genzFunctionName, "done", '\n')
        }

        # Combined plots for all cv runs
        if (num_cv_total > 1){
          combinedGPsd <- apply(predictionComb$meanValueGP, 1, sd)  
          combinedBARTsd <- apply(predictionComb$meanValueBART, 1, sd)
          combinedMeanMC <- apply(predictionComb$meanValueMonteCarlo, 1, mean)
          combinedMeanGP <- apply(predictionComb$meanValueGP, 1, mean)
          combinedMeanBART <- apply(predictionComb$meanValueBART, 1, mean)
          combinedTimeGP <- apply(predictionComb$timeGP, 1, mean)
          combinedTimeBART <- apply(predictionComb$timeBART, 1, mean)

          ymin <- min( apply(cbind(combinedMeanBART, combinedMeanMC, combinedMeanGP), 2, FUN = min) )
          ymax <- max( apply(cbind(combinedMeanBART, combinedMeanMC, combinedMeanGP), 2, FUN = max) )
          ymin_err <- min( apply(cbind(abs(combinedMeanBART - real), abs(combinedMeanMC - real), abs(combinedMeanGP - real)), 2, FUN = min) )
          ymax_err <- max( apply(cbind(abs(combinedMeanBART - real), abs(combinedMeanMC - real), abs(combinedMeanGP - real)), 2, FUN = max) )
          ymax_time <- max( apply(cbind(combinedTimeBART, combinedTimeGP), 2, FUN = max) )

          pdf(paste("Figures", "/combined_", toString(whichGenz), "Dim", dim, ".pdf", sep = ""), width = 8.5, height = 5)
          par(mfrow = c(1,2), pty = "s")
          plot(integrals$epochs,
              abs(combinedMeanMC - real),
              ty="l",
              ylab = "Absolute Error",
              xlab = "Number of training points",
              col = "chartreuse4",
              lwd = 1.5,
              ylim = c(ymin_err, ymax_err)
          )
          abline(h=0)
          points(integrals$epochs, abs(combinedMeanGP - real), ty = "l", col = "dodgerblue", lwd = 1.5)
          points(integrals$epochs, abs(combinedMeanBART - real), ty="l", col = "orangered", lwd = 1.5)
          legend("topright", legend=c("BART-Int", "GP-BQ", "MI"),
                col=c("orangered", "dodgerblue", "chartreuse4"), cex=1.2, lty = c(1,1,1,1))
          
          plot(integrals$epochs, 
              combinedMeanMC, 
              ty="l", 
              ylab = "Integral",
              xlab = "Number of training points",
              col = "chartreuse4",
              ylim = c(ymin, ymax),
              lwd = 1.5
          )
          polygon(c(integrals$epochs, rev(integrals$epochs)), 
                  c(
                    combinedMeanGP + 2*combinedGPsd, 
                    rev(combinedMeanGP - 2*combinedGPsd)
                  ), 
                  col = adjustcolor("dodgerblue", alpha.f = 0.20), 
                  border = adjustcolor("dodgerblue", alpha.f = 0.20),
                  lty = "solid")
          polygon(c(integrals$epochs, rev(integrals$epochs)), 
                  c(
                    combinedMeanBART + 2*combinedBARTsd,
                    rev(combinedMeanBART - 2*combinedBARTsd)
                  ), 
                  col = adjustcolor("orangered", alpha.f = 0.20), 
                  border = adjustcolor("orangered", alpha.f = 0.20),
                  lty = "solid")
          abline(h=integrals$actual)
          points(integrals$epochs, combinedMeanGP, ty="l", col = "dodgerblue", lwd = 1.5)
          points(integrals$epochs, combinedMeanBART, ty="l", col = "orangered", lwd = 1.5)

          # Close file
          dev.off()

          # Runtime plot
          pdf(paste("Figures", "/combined_time_", toString(whichGenz), "Dim", dim, ".pdf", sep = ""), width = 8.5, height = 5)
          plot(integrals$epochs,
              combinedTimeGP,
              ty="b",
              ylab = "Runtime",
              xlab = "Number of training points",
              col = "dodgerblue",
              # lwd = 1.5,
              ylim = c(0, ymax_time * 2)
          )
          points(integrals$epochs, combinedTimeBART, ty="b", col = "orangered")
          legend("topright", legend=c("BART-Int", "GP-BQ"),
                col=c("orangered", "dodgerblue"), cex=1.2, lty = c(1,1,1,1))

          # Close file
          dev.off()
        }
      }
    } 
  }

  cat(paste("Done plotting:", plotRootPath, "\n"))
}
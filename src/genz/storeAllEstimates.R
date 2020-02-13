# Store all estimates and best methods in separate csv's.
# To run:
#     In Unix terminal: 
#     Rscript {ROOT}/src/genz/drawGraphs.R
# To see results:
#     Results stored in {ROOT}/report/genz/resultsAll.csv and bestmethods.csv

methods <- c("BART", "MI", "GP")
bestMethod <- matrix(NA, ncol = 6, nrow = 7)
rownames(bestMethod) <- c("cont", "copeak", "disc", "gaussian", "oscil", "prpeak", "step")
colnames(bestMethod) <- c("1", "2", "3", "5", "10", "20")

# Store all results in a single csv file
resultsAll <- matrix(NA, ncol = 7, nrow = 28)
rownames(resultsAll) <- matrix(c(rep("cont",4), rep("copeak", 4), rep("disc", 4), rep("gaussian", 4), 
                          rep("oscil", 4), rep("prpeak", 4), rep("step", 4)), ncol = 1)
colnames(resultsAll) <- c("Methods", "1", "2", "3", "5", "10", "20")   
resultsAll[, 1] <- c("True", "BART", "MI", "GP")  

rmseValues <- matrix(NA, ncol = 7, nrow = 21)
rownames(rmseValues) <- matrix(c(rep("cont",3), rep("copeak", 3), rep("disc", 3), rep("gaussian", 3), 
                          rep("oscil", 3), rep("prpeak", 3),  rep("step", 3)), ncol = 1)
colnames(rmseValues) <- c("Methods", "1", "2", "3", "5", "10", "20") 
rmseValues[, 1] <- c("BART", "MI", "GP") 

# RMSE function
RMSE <- function(x, y){
  return( sqrt(mean((x - y)^2)) )
}

for (i in c(1, 2, 4, 5, 6, 7)){
    for (j in c(1, 5, 6)){

        # global parameters: dimension
        dimensionsList <- c(1,2,3,5,10,20)
        dim <- dimensionsList[j]
        whichGenz <- i
        
        # Skip if dim = 1 or 2 for discontinuous
        # if (whichGenz == 3 & dim == 1 ) { 
        #   next 
        # } else if (whichGenz == 3 & dim == 2) {
        #   next
        # } else if (whichGenz == 3 & dim == 3) {
        #   next
        # }
        if (whichGenz == 3) { 
            next
        }
        # Find Genz function
        if (whichGenz == 1) { genzFunctionName <-  "cont" }
        if (whichGenz == 2) { genzFunctionName <-  "copeak" }
        if (whichGenz == 3) { genzFunctionName <-  "disc" }
        if (whichGenz == 4) { genzFunctionName <-  "gaussian" }
        if (whichGenz == 5) { genzFunctionName <-  "oscil" }
        if (whichGenz == 6) { genzFunctionName <-  "prpeak" }
        if (whichGenz == 7) { genzFunctionName <-  "step" }
        if (whichGenz == 8) { genzFunctionName <-  "mix" }

        # Set path for estimated integral values
        fileName <- paste(toString(genzFunctionName), 'Dim', toString(dim), "", "Uniform",'.csv', sep='')
        filePath <- paste('../../mlbox/results_3_matern_k2/results/genz', toString(whichGenz), fileName, sep='/')

        # Retrieve estimated integral values
        integrals <- read.csv(filePath, header=TRUE, sep=",", stringsAsFactors = FALSE)
        predictionBART <- data.frame("meanValueBART" = integrals[, 2], "standardDeviationBART" = integrals[, 3])
        predictionMonteCarlo <- data.frame("meanValueMonteCarlo" = integrals[, 4], "standardDeviationMonteCarlo" = integrals[, 5])
        predictionGPBQ <-  data.frame("meanValueGP" = integrals[, 6], "standardDeviationGP" = integrals[, 7])

        # Retrieve analytical integral values
        whichDimension <- which(dim == dimensionsList)
        analyticalIntegrals <- read.csv("../../results/genz/integrals.csv", header = FALSE)
        real <- analyticalIntegrals[whichGenz, whichDimension]

        # Compute abs error
        epoch <- nrow(integrals)
        absDifference <- abs(c(predictionBART$meanValueBART[epoch], predictionMonteCarlo$meanValueMonteCarlo[epoch], 
                           predictionGPBQ$meanValueGP[epoch]) - real)
        bestMethod[i, j] <- methods[which(absDifference == min(absDifference))[1]]

        # Store all estimates in a single csv
        resultsAll[( (4*i - 3) : (4*i) ), (j + 1)] <- signif(c(real, predictionBART$meanValueBART[epoch], predictionMonteCarlo$meanValueMonteCarlo[epoch], 
                                          predictionGPBQ$meanValueGP[epoch]), 3)

        # Store all RMSE
        rmse <- c(RMSE(predictionBART$meanValueBART[epoch], real),
                  RMSE(predictionMonteCarlo$meanValueMonteCarlo[epoch], real),
                  RMSE(predictionGPBQ$meanValueGP[epoch], real)
                  )
        rmseValues[( (3*i - 2) : (3*i) ), (j + 1)] <- signif(rmse, 3)
    }
    cat(genzFunctionName, "done", '\n')
}

write.csv(bestMethod, file = "../../results/genz/bestMethods.csv")
write.csv(resultsAll, file = "../../results/genz/allEstimates.csv")
write.csv(rmseValues, file = "../../results/genz/rmseValues.csv")

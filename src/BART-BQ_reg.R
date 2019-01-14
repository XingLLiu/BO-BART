library(dbarts)
library(data.tree)
source("/Users/christine798/Documents/GitHub/BO-BART/src/BARTBQ.R")

rescale <- function(myData) {
# rescale variables in train set to apply pdf of uniform distribution

    for (i in 1:ncol(myData)) {

        myData[, i] <- (myData[, i] - min(myData[, i])) / (max(myData[, i]) - min(myData[, i]))
        #myData[, i] <- (myData[, i] - min(myData[, i])) / (max(myData[, i]) - min(myData[, i])) - 0.5
        #myData[, i] <- pnorm(myData[, i])
        
    }

    return(myData)

}

singleTreePrediction <- function(tree, x, cutPoints) {
# compute a prediction given by one tree in one posterior draw

    # return prediction when reach terminal node
    if (is.null(tree$leftChild) == FALSE) {
        
        # decision rule on specific element
        decisionRule <- cutPoints[[tree$splitVar]][tree$splitIndex]

        # locate datapoint to next child
        if (x[tree$splitVar] <= decisionRule) {

            singleTreePrediction(tree$leftChild, x, cutPoints)

        } else {
            
            singleTreePrediction(tree$rightChild, x, cutPoints)
       
        }

    } else {
        
        # return single tree prediction
        return(tree$terminal_probability * tree$mu)

    }

}  


singlePosteriorPrediction <- function(drawNum, model, x, nTree) {
# compute prediction given by a posterior draw (first quadrature)

    cutPoints <- dbarts:::createCutPoints(model$fit)
    cut <- array(c(0, 1), c(2, dim))

    posteriorPrediction <- NULL

    for (i in 1:nTree) {
  
        treeList <- getTree(model$fit, 1, drawNum, i)
  
        tree <- FromListSimple(treeList) 
  
        #Modify tree by the functions written above
        tree <- fillProbabilityForNode(tree, cutPoints, cut)
        tree <- terminalProbabilityStore(tree)

        singlePosteriorPrediction <- c(posteriorPrediction, singleTreePrediction(tree, x, cutPoints))

    }
    
    return(mean(singlePosteriorPrediction))

}

posteriorPrediction <- function(testNum, model, test, nTree) {

    x <- test[testNum, ]
    
    nDraw <- dim(model$fit$state[[1]]@savedTreeFits)[3]
    drawNum <- seq(1, nDraw, length.out=nDraw)

    #Extra variables
    var <- list(model, x, nTree)
  
    #Calculate integration over all trees in the draw by mapply
    posteriorPrediction <- mapply(singlePosteriorPrediction, drawNum, MoreArgs=var)

    return(sum(posteriorPrediction))

}



bartBQPredict <- function(model, test) {
# predict response to a new datapoint with BART-BQ

    nTree <- dim(model$fit$state[[1]]@savedTreeFits)[2]
    
    nTest <- nrow(test)
    testNum <- seq(1, nTest, length.out=nTest) 
  
    #Extra variables
    var <- list(model, test, nTree)
  
    #Calculate integration over all trees in the draw by mapply
    prediction <- mapply(posteriorPrediction, testNum, MoreArgs=var)

    return(prediction)

}




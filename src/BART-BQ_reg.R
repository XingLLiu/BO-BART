library(dbarts)
library(data.tree)
library(Rfast)
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

scale <- function(mu, max, min) {

    mu <- (mu + 0.5) * (max - min) + min

    return(mu)

}

fillProbabilityForNode_update <- function(oneTree, cutPoints, cut) 
# Drop data set into the tree and assign them to different nodes 
{
  if ( !is.null(oneTree$leftChild) ) {
    
    decisionRule <- cutPoints[[oneTree$splitVar]][oneTree$splitIndex]

    oneTree$leftChild$probability <- (decisionRule - cut[1, oneTree$splitVar]) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
    #cat("left prob", oneTree$leftChild$probability, "\n")
    oneTree$rightChild$probability <- (cut[2, oneTree$splitVar] - decisionRule) / (cut[2, oneTree$splitVar] - cut[1, oneTree$splitVar])
    #cat("right prob", oneTree$rightChild$probability, "\n")

    range <- cut[, oneTree$splitVar]
    cut[, oneTree$splitVar] <- c(range[1], decisionRule)
    
    # if (cut[1, oneTree$splitVar] == cut[2, oneTree$splitVar]) {
        
    #     print(c("left", cut[, oneTree$splitVar]))

    # }
    
    fillProbabilityForNode_update(oneTree$leftChild, cutPoints, cut)
    
    cut[, oneTree$splitVar] <- c(decisionRule, range[2])

    # if (cut[1, oneTree$splitVar] == cut[2, oneTree$splitVar]) {
        
    #     print(c("right", cut[, oneTree$splitVar]))

    # }
    
    fillProbabilityForNode_update(oneTree$rightChild, cutPoints, cut)
    
  } else if( is.null(oneTree$probability) ) {

    oneTree$probability <- 1

  }
  
  return (oneTree)
}

singleTreePrediction <- function(tree, x, cutPoints, max, min) {
# compute a prediction given by one tree in one posterior draw

    # return prediction when reach terminal node
    if (is.null(tree$leftChild) == FALSE) {
        
        # decision rule on specific element
        decisionRule <- cutPoints[[tree$splitVar]][tree$splitIndex]
        
        # locate datapoint to next child
        if (x[tree$splitVar] <= decisionRule) {

            singleTreePrediction(tree$leftChild, x, cutPoints, max, min)
           
        } else {
            
            singleTreePrediction(tree$rightChild, x, cutPoints, max, min)
            
        }

    } else {

        #cat("prob", tree$terminal_probability, "\n")
        #cat("mu", tree$mu, "\n")
        #return single tree prediction
        #return(tree$terminal_probability * scale(tree$mu, max, min))
        return(scale(tree$mu, max, min))

    }

}  

singlePosteriorPrediction <- function(drawNum, model, x, nTree, cut, max, min) {
# compute prediction given by a posterior draw (first quadrature)

    cutPoints <- dbarts:::createCutPoints(model$fit)

    singlePosteriorPrediction <- rep(NULL, nTree)

    for (i in 1:nTree) {
  
        treeList <- getTree(model$fit, 1, drawNum, i)
        
        tree <- FromListSimple(treeList) 
        
        #Modify tree by the functions written above
        tree <- fillProbabilityForNode_update(tree, cutPoints, cut)
        tree <- terminalProbabilityStore(tree)
    
        singlePosteriorPrediction[i] <- singleTreePrediction(tree, x, cutPoints, max, min)
    
    }
    
    return(mean(singlePosteriorPrediction))

}

posteriorPrediction <- function(testNum, model, test, nTree, cut, max, min) {
    
    x <- test[testNum, ]
    
    nDraw <- dim(model$fit$state[[1]]@savedTreeFits)[3]
    drawNum <- seq(1, nDraw, length.out = nDraw)

    #Extra variables
    var <- list(model, x, nTree, cut, max, min)
    
    #Calculate integration over all trees in the draw by mapply
    posteriorPrediction <- mapply(singlePosteriorPrediction, drawNum, MoreArgs = var)

    return(mean(posteriorPrediction))

}

bartBQPredict <- function(model, train, test, max, min) {
# predict response to a new datapoint with BART-BQ

    cut_origin <- rbind(colMins(as.matrix(train[, 1:(ncol(train)-1)]), value = TRUE), colMaxs(as.matrix(train[, 1:(ncol(train)-1)]), value = TRUE))
    #cut_origin <- array(c(0, 1), c(2, 2))
    
    #print(c("cut", cut))
    nTree <- dim(model$fit$state[[1]]@savedTreeFits)[2]
    
    nTest <- nrow(test)
    testNum <- seq(1, nTest, length.out = nTest) 

    #Extra variables
    var <- list(model, test, nTree, cut_origin, max, min)

    #Calculate integration over all trees in the draw by mapply
    prediction <- mapply(posteriorPrediction, testNum, MoreArgs = var)

    return(prediction)

}




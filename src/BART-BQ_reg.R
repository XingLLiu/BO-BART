library(dbarts)
library(data.tree)
library(hydroGOF)
source("/Users/christine798/Documents/GitHub/BO-BART/src/BARTBQ.R")

netTrain <- read.csv("~/Documents/2017-2018/UROP/report/network/net_train.csv")
netTrain <- net_train[-which(net_train$V27 == 0), ]
netTest <- read.csv("~/Documents/2017-2018/UROP/report/network/net_test.csv")
netTest <- net_test[-which(net_test$V27 == 0), ]

rescale <- function(myData) {
# rescale variables in train set to apply pdf of uniform distribution

    for (i in 1:ncol(myData)) {

        myData[, i] <- (myData[, i] - min(myData[, i])) / (max(myData[, i]) - min(myData[, i]))
        
    }

    return(my_data)

}

singleTreePrediction <- function(tree, x, cutPoints) {
# compute a prediction given by one tree in one posterior draw

    # return prediction when reach terminal node
    if (is.null(tree$leftChild) == FALSE) {
        
        # decision rule on specific element
        decisionRule <- cutPoints[[tree$splitVar]][tree$splitIndex]

        # locate datapoint to next child
        if (x[tree$splitVar] <= decisionRule) {

            singleTreePrediction(tree$leftChild, x)

        } else {
            
            singleTreePrediction(tree$rightChild, x)
       
        }

    } else {
        
        # return single tree prediction
        return(tree$terminal_probability * tree$mu)

    }

}  


posteriorPrediction <- function(model, drawNum, nTree, test, testNum) {
# compute prediction given by a posterior draw (first quadrature)

    x <- test[testNum, ]

    cutPoints<-dbarts:::createCutPoints(model$fit)

    posteriorPrediction <- NULL

    for (i in 1:ntree) {
  
        treeList <- getTree(model$fit, 1, drawNum, i)
  
        tree <- FromListSimple(treeList) 
  
        #Modify tree by the functions written above
        tree <- fillProbabilityForNode(tree, cutPoints, cut)
        tree <- terminalProbabilityStore(tree)

        posteriorPrediction <- c(posteriorPrediction, singleTreePrediction(tree, x, cutPoints))

    }
    
    return(mean(posteriorPrediction))

}

predict <- function(model, test) {
# predict response to a new datapoint with BART-BQ

    nTree <- dim(model$fit$state[[1]]@savedTreeFits)[2]
    drawNum <- dim(model$fit$state[[1]]@savedTreeFits)[3]
    dim <- ncol(test)

    nTest <- nrow(test)

    testNum <- seq(1, nTest, length.out=nTest)
  
    #Extra variables
    var <- list(model, drawNum, nTree, test)
  
    #Calculate integration over all trees in the draw by mapply
    posteriorPrediction <- sum(unlist(mapply(singleTreePrediction, testNum, MoreArgs=var) ) )

}

pred <- colMeans(predict(bart_model, test))


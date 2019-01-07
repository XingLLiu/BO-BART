############
#Bayesian Quadrature with Gaussian Process
############
library('kernlab')

# Function to be learned by BART
# ! not well-behaved when dim is large
sampleTrainData <- function(dim=2, N=10, FUN)  
# Sample train data X and computes corresponding observed (binary) values Y 
# input:
#     dim: 
#     N:
#     FUN: function to be learned
{    
    # X <- rmvnorm(n = N, mean = rep(0, dim), sigma = diag(dim))
    # probability <- apply(X, 1, pmvnorm, diag(dim))
    X <- randomLHS(N, dim)
    probability <- rep(NA, N)
    for (i in 1:N){
      # probability[i] <- pmvnorm(upper = X[i, ], sigma = diag(dim))
      probability[i] <- sum(X[i, ])
    }
    Y <- as.numeric(probability < 1)   # Binary output; arbitrarily chosen threshold value

    return(data.frame(Y, X))
}


# GP classification model
trainData <- sampleTrainData(dim = 3, N = 120)
gpFit <- gausspr(factor(Y)~., data = trainData[1:100, ], type = 'classification', kernel = 'rbfdot', kpar = list('sigma' = 1))

ytest <- as.numeric(predict(gpFit, trainData[101:120, -1])) - 1
# ytest <- predict(gpFit, trainData[1501:2000, -1], type = 'probabilities')[1,] > 0.5
plot(trainData[101:120,1])
points(ytest, col = 'red')

sum(as.numeric(ytest) == trainData[101:120,1])/20


# Existing dataset
data(iris)
test <- gausspr(Species~.,data=iris[1:120, ],var=2)
alpha(test)

# predict on the training set
plotPred <- hist(as.numeric(predict(test,iris[121:150, -5])))
plltTrue <- hist(as.numeric(iris[121:150, 5]))

sum(as.numeric(predict(test,iris[121:150, -5])) == as.numeric(iris[121:150, 5]))/30



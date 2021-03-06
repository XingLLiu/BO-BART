# !/usr/bin/env R
# Load required packages
library(MASS)
library(lhs)
library(data.tree)
library(dbarts)
library(matrixStats)
library(mvtnorm)

# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

# global parameters: dimension

# global parameters: dimension
args <- as.double(commandArgs(TRUE))
dim <- args[1]
num_iterations <- args[2]
whichGenz <- args[3]
cat("Sequential: ", args[4])
# turn on/off sequential design
# 1 denotes TRUE to sequential
# 0 denotes FALSE to sequential
cat("\nBegin testing:\n")
if (args[4] == 1 | is.na(args[4])) {
  sequential <- TRUE
  print("Sequantial design set to TRUE.")
} else {
  sequential <- FALSE
  print("Sequantial design set to FALSE.")
}

if (num_iterations == 1) { stop("NEED MORE THAN 1 ITERATION") }

print(c(dim, num_iterations, whichGenz))
source("src/genz/genz.R") # genz function to test

if (whichGenz < 1 | whichGenz > 8) { stop("undefined genz function. Change 3rd argument to 1-8") }
if (whichGenz == 3 & dim == 1) { stop("incorrect dimension. Discrete Genz function only defined for dimension >= 2") } 

if (whichGenz == 1) { genz <- cont; genzFunctionName <-  deparse(substitute(cont)) }
if (whichGenz == 2) { genz <- copeak; genzFunctionName <-  deparse(substitute(copeak)) }
if (whichGenz == 3) { genz <- disc; genzFunctionName <-  deparse(substitute(disc)) }
if (whichGenz == 4) { genz <- gaussian; genzFunctionName <-  deparse(substitute(gaussian)) }
if (whichGenz == 5) { genz <- oscil; genzFunctionName <-  deparse(substitute(oscil)) }
if (whichGenz == 6) { genz <- prpeak; genzFunctionName <-  deparse(substitute(prpeak)) }
if (whichGenz == 7) { genz <- step; genzFunctionName <-  deparse(substitute(step)) }
if (whichGenz == 8) { genz <- mix; genzFunctionName <-  deparse(substitute(mix)) }

print("Testing with: %s" %--% genzFunctionName)

# prepare training dataset
trainX <- replicate(dim, runif(500))
if (whichGenz == 7) {
  trainY <- genz(trainX)
} else {
  trainY <- genz(trainX)
}

# Bayesian Quadrature method
# set number of new query points using sequential design
source("src/BARTBQ.R")
t0 <- proc.time()
posterior_model <- BART_posterior(dim, trainX, trainY, num_iterations, FUN = genz, sequential)
t1 <- proc.time()
bartTime <- (t1 - t0)[[1]]

x_plot <- replicate(dim, runif(500))
x_order <- order(x_plot)
y_pred <- predict(posterior_model,x_plot)
y_pred_mean <- colMeans(y_pred)
y_pred_sd <- sqrt(colVars(y_pred))

plot(trainX, trainY)
polygon(c(x_plot[x_order], rev(x_plot[x_order])),
        c(y_pred_mean[x_order] - y_pred_sd[x_order], rev(y_pred_mean[x_order] + y_pred_sd[x_order])),
        col = rgb(0, 0, 1, 0.6),
        border = FALSE)
points(x_plot, y_pred_mean, col = "red", cex=0.2)

if (!sequential){
  figName <- "Figures/%s/drawBART%s%sDimNoSequential.pdf" %--% c(whichGenz, genzFunctionName, dim)
  csvName <- "Figures/%s/drawBART%s%sDimNoSequential.csv" %--% c(whichGenz, genzFunctionName, dim)
  groundTruthName <- "Figures/%s/trainDrawBart%s%sDimNoSequential.csv" %--% c(whichGenz, genzFunctionName, dim)
} else {
  figName <- "Figures/%s/drawBART%s%sDim.pdf" %--% c(whichGenz, genzFunctionName, dim)
  csvName <- "Figures/%s/drawBART%s%sDim.csv" %--% c(whichGenz, genzFunctionName, dim)
  groundTruthName <- "Figures/%s/trainDrawBart%s%sDim.csv" %--% c(whichGenz, genzFunctionName, dim)
}

results <- data.frame(
  "x_plot" = x_plot,
  "y_pred" = y_pred
)
groundTruth <- data.frame(
  "trainX" = trainX,
  "trainY" = trainY
)
write.csv(results, file = csvName, row.names=FALSE)
write.csv(groundTruth, file = groundTruthName, row.names=FALSE)

print("Begin Plots")
# 1. Open jpeg file
pdf(figName, width = 10, height = 11)
# 2. Create the plot
par(mfrow = c(1,2), pty = "s")
plot(trainX, trainY)
points(x_plot, y_pred, col = "red", cex=0.2)

legend("topleft", legend=c("BART BQ", "GP BQ", "Actual"),
       col=c("red", "green", "black"), cex=0.8, lty = c(1,1,1,1))

# 3. Close the file
dev.off()

print("Please check {ROOT}/%s for plots" %--% figName)



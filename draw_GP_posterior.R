# !/usr/bin/env R
# Load required packages
library(MASS)
library(cubature)
library(lhs)
library(matrixStats)
library(mvtnorm)
library(doParallel)
library(kernlab)

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
trainX <- randomLHS(500, dim)
trainY <- genz(trainX)
# plot(trainX, trainY)

# Bayesian Quadrature with Gaussian Process
print("Begin Gaussian Process Integration")
library(reticulate)
source("src/optimise_gp.R")
lengthscale <- optimise_gp_r(trainX, trainY, kernel="rbf", epochs=500)

source("src/GPBQ.R")
t0 <- proc.time()
# need to add in function to optimise the hyperparameters
GPResults <- computeGPBQ(trainX, trainY, dim, epochs = num_iterations-1, N=500, FUN = genz, lengthscale,sequential)  
t1 <- proc.time()
GPTime <- (t1 - t0)[[1]]

K <- GPResults$K
X <- GPResults$X
Y <- GPResults$Y
x_plot <- randomLHS(1000, dim)
k_xstar_x <- kernelMatrix(rbfdot(.5/lengthscale^2), matrix(x_plot, ncol=1), X)
k_xstar_xstar <- kernelMatrix(rbfdot(.5/lengthscale^2), 
                              matrix(x_plot, ncol=1), 
                              matrix(x_plot, ncol=1))
jitter = 1e-6
K_inv <- solve(K + diag(jitter, nrow(K)))

gp_post_mean <- k_xstar_x %*% K_inv %*% Y
gp_post_cov <- k_xstar_xstar - k_xstar_x %*% K_inv %*% t(k_xstar_x)

plot(trainX, trainY, ylim=c(-0.5, 1))
points(x_plot, gp_post_mean, col = "red", cex=0.2)

if (!sequential){
  figName <- "Figures/%s/drawGP%s%sDimNoSequential.pdf" %--% c(whichGenz, genzFunctionName, dim)
  csvName <- "Figures/%s/drawGP%s%sDimNoSequential.csv" %--% c(whichGenz, genzFunctionName, dim)
  groundTruthName <- "Figures/%s/trainDrawGP%s%sDimNoSequential.csv" %--% c(whichGenz, genzFunctionName, dim)
} else {
  figName <- "Figures/%s/drawGP%s%sDim.pdf" %--% c(whichGenz, genzFunctionName, dim)
  csvName <- "Figures/%s/drawGP%s%sDim.csv" %--% c(whichGenz, genzFunctionName, dim)
  groundTruthName <- "Figures/%s/trainDrawGP%s%sDim.csv" %--% c(whichGenz, genzFunctionName, dim)
}

results <- data.frame(
  "x_plot" = x_plot,
  "y_pred" = gp_post_mean
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
points(x_plot, gp_post_mean, col = "red", cex=0.2)

legend("topleft", legend=c("BART BQ", "GP BQ", "Actual"),
       col=c("red", "green", "black"), cex=0.8, lty = c(1,1,1,1))

# 3. Close the file
dev.off()

print("Please check {ROOT}/Figures/%s for plots" %--% figName)



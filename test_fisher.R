source("src/genz/fisher_integrands.R")
dim <- 10
# C <- replicate(dim, runif(9, 0.1, 0.9))
C <- rep(0.45958094818045914, dim)
R <- rep(0.13332051229486602, dim)
H <- rep(1.5020584867022158, dim)
F <- rep(4.624266954512351, dim)
P <- rep(1, dim)


# trainX <- matrix(seq(0, 1, 1/1000))
fisher_function <- create_fisher_function(C, R, H, F, P, dim)
# trainY <- fisher_function(trainX)

# plot(trainX[order(trainX),], trainY[order(trainX)], ty= "l")

groundX <- replicate(dim, seq(0, 1, 1/1e7))
groundY <- fisher_function(groundX)
real <- mean(groundY)
print(real)

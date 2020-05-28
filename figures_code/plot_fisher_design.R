load("results/fisher_function/PaperDim1Uniform_20_10.RData")

# trainX <- results_models$BART$trainData$trainX[1:20]
# trainY <- results_models$BART$trainData$trainY[1:20]
# newX <- results_models$BART$trainData$trainX[21:40]
# newY <- results_models$BART$trainData$trainY[21:40]
# gp <- ""

trainX <- results_models$GP$X[1:20]
trainY <- results_models$GP$Y[1:20]
newX <- results_models$GP$X[21:40]
newY <- results_models$GP$Y[21:40]
gp <- "GP"

dim <- 1
num_data <- 20
num_iterations <- 20*dim
measure <- "uniform"
source("src/genz/fisher_integrands.R")
C <- c(0.1706093, 0.5319923, 0.7117816)
R <- c(0.6786221, 0.7207544, 0.5120249)
H <- c(3.029867, 3.065427, 3.114357)
F <- c(4.607562, 4.526243, 2.221768)
P <- c(1, 0, 1)

fisher_function <- create_fisher_function(C[1:dim], R[1:dim], H[1:dim], F[1:dim], P[1:dim], dim)
plotX <- as.matrix(seq(0, 1, 1/1000))

# plot 5
pdf(paste0("figures_code/", "fisher_design_", gp, 5, ".pdf"), width = 5, height = 5)
par(pty="s")
plot(plotX, fisher_function(plotX), ty="l", ylim = c(-1, 3.5), xlab="x", ylab = "y",
     cex.axis = 2, cex.lab = 2)
points(trainX[order(trainX)], trainY[order(trainX)], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
points(newX[order(newX)][1:5], newY[order(newX)][1:5], col = "orangered", bg='orangered', pch=21, lwd=3)
legend("bottomright", legend=c("Design Points", "Sequential Points"),
       col=c("dodgerblue", "orangered"), cex=1.8, bty="n",pch = c(19, 19))
dev.off()


# plot 10
pdf(paste0("figures_code/", "fisher_design_", gp, 10, ".pdf"), width = 5, height = 5)
par(pty="s")
plot(plotX, fisher_function(plotX), ty="l", ylim = c(-1, 3.5), xlab="x", ylab = "y",
     cex.axis = 2, cex.lab = 2)
points(trainX[order(trainX)], trainY[order(trainX)], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
points(newX[order(newX)][1:10], newY[order(newX)][1:10], col = "orangered", bg='orangered', pch=21, lwd=3)
legend("bottomright", legend=c("Design Points", "Sequential Points"),
       col=c("dodgerblue", "orangered"), cex=1.8, bty="n",pch = c(19, 19))
dev.off()

# plot 20
pdf(paste0("figures_code/", "fisher_design_", gp, 20, ".pdf"), width = 5, height = 5)
par(pty="s")
plot(plotX, fisher_function(plotX), ty="l", ylim = c(-1, 3.5), xlab="x", ylab = "y",
     cex.axis = 2, cex.lab = 2)
points(trainX[order(trainX)], trainY[order(trainX)], col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3)
points(newX[order(newX)][1:20], newY[order(newX)][1:20], col = "orangered", bg='orangered', pch=21, lwd=3)
legend("bottomright", legend=c("Design Points", "Sequential Points"),
       col=c("dodgerblue", "orangered"), cex=1.8, bty="n",pch = c(19, 19))
dev.off()



####### 2 Dimensions
load("results/fisher_function/PaperDim2Uniform_20_16.RData")

# trainX1 <- results_models$BART$trainData$X1[1:20]
# trainX2 <- results_models$BART$trainData$X2[1:20]
# newX1 <- results_models$BART$trainData$X1[21:60]
# newX2 <- results_models$BART$trainData$X2[21:60]
dim=2
trainX1 <- results_models$GP$X[1:20, 1]
trainX2 <- results_models$GP$X[1:20, 2]
newX1 <- results_models$GP$X[21:59,1]
newX2 <- results_models$GP$X[21:59,2]
fisher_function <- create_fisher_function(C[1:dim], R[1:dim], H[1:dim], F[1:dim], P[1:dim], dim)

# plot 5
pdf(paste0("figures_code/", "fisher_designDim2_", gp, 10, ".pdf"), width = 5, height = 5)
par(pty="s")
plot(trainX1, trainX2, col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3, ylim=c(0,1.5), xlim=c(0,1.5))
abline(h=1)
abline(v=1)
points(newX1[1:10], newX2[1:10], col = "orangered", bg='orangered', pch=21, lwd=3)
legend("topleft", legend=c("Design Points", "Sequential Points"),
       col=c("dodgerblue", "orangered"), cex=1.8, bty="n",pch = c(19, 19))
dev.off()


# plot 10
pdf(paste0("figures_code/", "fisher_designDim2_", gp, 20, ".pdf"), width = 5, height = 5)
par(pty="s")
plot(trainX1, trainX2, col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3, ylim=c(0,1.5), xlim=c(0,1.5))
abline(h=1)
abline(v=1)
points(newX1[1:20], newX2[1:20], col = "orangered", bg='orangered', pch=21, lwd=3)
legend("topleft", legend=c("Design Points", "Sequential Points"),
       col=c("dodgerblue", "orangered"), cex=1.8, bty="n",pch = c(19, 19))
dev.off()


# plot 20
pdf(paste0("figures_code/", "fisher_design_Dim2", gp, 40, ".pdf"), width = 5, height = 5)
par(pty="s")
plot(trainX1, trainX2, col = "dodgerblue", bg='dodgerblue', pch=21, lwd=3, ylim=c(0,1.5), xlim=c(0,1.5))
abline(h=1)
abline(v=1)
points(newX1[1:40], newX2[1:40], col = "orangered", bg='orangered', pch=21, lwd=3)
legend("topleft", legend=c("Design Points", "Sequential Points"),
       col=c("dodgerblue", "orangered"), cex=1.8, bty="n",pch = c(19, 19))
dev.off()


### heat map
library(ggplot2)
library(graphics)
jet.colors <- colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan", "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))
plotX <- as.matrix(seq(0, 1, 1/100))

x <- plotX
y <- plotX
df <- expand.grid(X=x, Y=y)
df$Z <- fisher_function(df)

pdf(paste0("figures_code/", "fisher_design_Dim2_heatmap", gp, 10, ".pdf"), width = 5, height = 5)
par(pty="s")
df_train <- data.frame(X1=trainX1, X2=trainX2)
df_seq <- data.frame(X1=newX1[1:10], X2=newX2[1:10])
ggplot() + 
        geom_tile(data=df, aes(x=X, y=Y, fill=Z)) +
        geom_point(data=df_train, aes(x=X1, y=X2), colour="dodgerblue") + 
        geom_point(data=df_seq, aes(x=X1, y=X2), colour="orangered") 
dev.off()

pdf(paste0("figures_code/", "fisher_design_Dim2_heatmap", gp, 20, ".pdf"), width = 5, height = 5)
par(pty="s")
df_train <- data.frame(X1=trainX1, X2=trainX2)
df_seq <- data.frame(X1=newX1[1:20], X2=newX2[1:20])
ggplot() + 
        geom_tile(data=df, aes(x=X, y=Y, fill=Z)) +
        geom_point(data=df_train, aes(x=X1, y=X2), colour="dodgerblue") + 
        geom_point(data=df_seq, aes(x=X1, y=X2), colour="orangered") 
dev.off()


pdf(paste0("figures_code/", "fisher_design_Dim2_heatmap", gp, 40, ".pdf"), width = 5, height = 5)
par(pty="s")
df_train <- data.frame(X1=trainX1, X2=trainX2)
df_seq <- data.frame(X1=newX1[1:40], X2=newX2[1:40])
ggplot() + 
        geom_tile(data=df, aes(x=X, y=Y, fill=Z)) +
        geom_point(data=df_train, aes(x=X1, y=X2), colour="dodgerblue") + 
        geom_point(data=df_seq, aes(x=X1, y=X2), colour="orangered") 
dev.off()


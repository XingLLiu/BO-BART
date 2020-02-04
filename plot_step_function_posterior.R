bartres <- read.csv("mlbox/7/drawBARTstep1Dim.csv")
gpres <- read.csv("mlbox/7/drawGPstep1Dim.csv")
actual <- read.csv("mlbox/7/trainDrawBartstep1Dim.csv")
actualgp <- read.csv("mlbox/7/trainDrawGPstep1Dim.csv")

actual <- actual[order(actual$trainX),]
bartres <- bartres[order(bartres$x_plot),]
gpres <- gpres[order(gpres$x_plot),]
plot(bartres$x_plot, 
     bartres$y_pred_mean, 
     ty="l",
     col="red", 
     xlab = "x", 
     ylab = "y",
     main = "Step function 1 dim estimation",
     ylim=c(-0.2, 1.2))
points(gpres$x_plot, gpres$y_pred_mean, ty="l", col="blue")
points(actual$trainX, actual$trainY, col="black", cex=0.2)
legend("topleft", legend=c("BART BQ", "GP BQ", "Actual"),
       col=c("red", "blue", "black"), cex=2, lty = c(1,1,1))
points(actualgp$trainX, actualgp$trainY)

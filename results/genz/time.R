times <- read.csv("GPTime.csv")
times <- times/3600


par(mfrow = c(2,3))
plot(c(1,2,3,5,10,20), c(times[1,]$X1, times[1,]$X2, times[1,]$X3, times[1,]$X5, times[1,]$X10, times[1,]$X20), type="l",
              ylab = "time s", xlab = "dimension")

for (i in 2:6) {
  plot(c(1,2,3,5,10,20), c(times[i,]$X1, times[i,]$X2, times[i,]$X3, times[i,]$X5, times[i,]$X10, times[i,]$X20),
     ylab = "time s", xlab = "dimension", type = "l")
}


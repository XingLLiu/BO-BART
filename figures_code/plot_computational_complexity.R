# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}

nums <- c(10,20,100,500,1000,2000,3000,5000,7000,10000)
bart_times <- c()
mi_times <- c()
gp_times <- c()

i=1
for (num in nums) {
  filename <- "computational_complexity/computational_complexity_stepDim1NoSequentialUniform_num%s.csv" %--% num
  df <- read.csv(filename)
  bart_times[i] <- df$runtimeBART
  mi_times[i] <- df$runtimeMI
  gp_times[i] <- df$runtimeGP
  i <- i+1
}

par(mfrow=c(1,2), pty="s")
plot(nums, log(gp_times), xlab="N", ylab="log(Time) (s)", main="GPBQ", cex.lab=1.4, cex.main=1.9, cex.lab=1.4)
m <- lm(log(gp_times) ~ nums)
coeffs <-as.numeric(m$coefficients)
abline(a=coeffs[1], b=coeffs[2])

plot(nums, bart_times, xlab="N", ylab="Time (s)", main="BART-Int", cex.lab=1.4, cex.main=1.9, cex.lab=1.4)
m <- lm(bart_times~ nums)
coeffs <-as.numeric(m$coefficients)
abline(a=coeffs[1], b=coeffs[2])
# plot(nums, bart_times)

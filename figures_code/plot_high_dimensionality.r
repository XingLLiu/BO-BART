# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}


results <- data.frame(matrix(0, ncol = 6, nrow = 4))
colnames(results) <- c("BART-Int", "BART-Int-se", "GPBQ", "GPBQ-se", "MI", "MI-se")
dims = c(1, 10, 20, 100)
l = 1
for (dim in dims) {
  bart = c()
  gp = c()
  mi = c()
  for (num_cv in 1:5) {
    df <- read.csv("high_dimensionality/9/additive_gaussianDim%sUniform_%s.csv" %--% c(dim, num_cv))
    bart = c(bart, df$BARTMean)
    gp = c(gp, df$GPMean)
    mi = c(mi, df$MIMean)
  }
  results[l, 1] = mean(abs(bart - df$actual))
  results[l, 2] = sd(abs(bart - df$actual)) / sqrt(5)
  results[l, 3] = mean(abs(gp - df$actual))
  results[l, 4] = sd(abs(gp - df$actual)) / sqrt(5)
  results[l, 5] = mean(abs(mi - df$actual))
  results[l, 6] = sd(abs(mi - df$actual)) / sqrt(5)
  l = l+1
}
write.csv(results, "figures_code//high_dimensionality.csv")


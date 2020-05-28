# define string formatting
`%--%` <- function(x, y) 
  # from stack exchange:
  # https://stackoverflow.com/questions/46085274/is-there-a-string-formatting-operator-in-r-similar-to-pythons
{
  do.call(sprintf, c(list(x), y))
}


results <- data.frame(matrix(0, ncol = 8, nrow = 2))
colnames(results) <- c("BART-Int", "BART-Int-se", "GPBQ", "GPBQ-se", "MI", "MI-se", "dim", "n")
dims = c(1, 2, 3)
l = 1
for (num_data in c(20)) {
  for (dim in c(1)) {
    bart = c()
    gp = c()
    mi = c()
    for (num_cv in 1:20) {
      df <- read.csv("results/fisher_function/fisher_function/PaperDim%sUniform_%s_%s.csv" 
                     %--% c(dim, num_data, num_cv))
      bart = c(bart, df$BARTMean)
      gp = c(gp, df$GPMean)
      mi = c(mi, df$MIMean)
    }
  }
  results[l, 1] = mean(abs(bart - df$actual))
  results[l, 2] = sd(abs(bart - df$actual)) / sqrt(20)
  results[l, 3] = mean(abs(gp - df$actual))
  results[l, 4] = sd(abs(gp - df$actual)) / sqrt(20)
  results[l, 5] = mean(abs(mi - df$actual))
  results[l, 6] = sd(abs(mi - df$actual)) / sqrt(20)
  results[l, 7] = dim
  results[l, 8] = num_data*dim
  l = l + 1
}
for (num_data in c(20)) {
  for (dim in c(2)) {
    bart = c()
    gp = c()
    mi = c()
    for (num_cv in 1:20) {
      df <- read.csv("results/fisher_function/fisher_function/PaperDim%sUniform_%s_%s.csv" 
                     %--% c(dim, num_data, num_cv))
      bart = c(bart, df$BARTMean)
      gp = c(gp, df$GPMean)
      mi = c(mi, df$MIMean)
    }
  }
  results[l, 1] = mean(abs(bart - df$actual))
  results[l, 2] = sd(abs(bart - df$actual)) / sqrt(20)
  results[l, 3] = mean(abs(gp - df$actual))
  results[l, 4] = sd(abs(gp - df$actual)) / sqrt(20)
  results[l, 5] = mean(abs(mi - df$actual))
  results[l, 6] = sd(abs(mi - df$actual)) / sqrt(20)
  results[l, 7] = dim
  results[l, 8] = num_data*dim
  l = l + 1
}

for (num_data in c(20)) {
  for (dim in c(3)) {
    bart = c()
    gp = c()
    mi = c()
    for (num_cv in 1:20) {
      df <- read.csv("results/fisher_function/fisher_function/PaperDim%sUniform_%s_%s.csv" 
                     %--% c(dim, num_data, num_cv))
      bart = c(bart, df$BARTMean)
      gp = c(gp, df$GPMean)
      mi = c(mi, df$MIMean)
    }
  }
  results[l, 1] = mean(abs(bart - df$actual))
  results[l, 2] = sd(abs(bart - df$actual)) / sqrt(20)
  results[l, 3] = mean(abs(gp - df$actual))
  results[l, 4] = sd(abs(gp - df$actual)) / sqrt(20)
  results[l, 5] = mean(abs(mi - df$actual))
  results[l, 6] = sd(abs(mi - df$actual)) / sqrt(20)
  results[l, 7] = dim
  results[l, 8] = num_data*dim
  l = l + 1
}

write.csv(results, "figures_code/non_stationarity.csv" %--% c(dim, num_data))



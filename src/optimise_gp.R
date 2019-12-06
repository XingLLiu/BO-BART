library(reticulate)


optimise_gp <- function(trainX, trainY, init_lengthscale) 
{
  source_python("python/gp_tune.py")
  test(trainX)
}
#virtualenv_create()
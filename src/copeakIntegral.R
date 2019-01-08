copeakIntegral <- function(dim = NA)
# Computes the true integral of the d-dimensional copeak function
# See the appendix of the documentation for its derivation
# Input:
#   dim: dimension of the integral
# Output:
#   Value of integral over [0, 1]^d
{
  integral <- 0
  for (n in 1:dim){
    integral <- ( n * (n + 1) )^(-1) * choose((dim - 1), (n - 1)) * (-1)^(n + dim + 1) + integral
  }
  return(integral / factorial(dim))
}


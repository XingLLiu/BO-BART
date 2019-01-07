copeakIntegral <- function(dim = NA)
# Computes the true integral of the d-dimensional copeak function
# See the appendix of the documentation for its derivation
# Input:
#   dim: dimension of the integral
# Output:
{
    integral <- 0
    for (n in dim:1){
      integral <- (n + 1)^(-1) * choose((dim - 1), (n - 1)) * (-1)^(n - 1) - 
                n^(-1) * choose((dim - 1), (n - 1)) * (-1)^(n - 1) + 
                integral
    }
    return(integral * (-1)^(dim - 1) / factorial(dim))
}


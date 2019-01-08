auxiliaryFunction <- function(x=NA, dim=NA)
# Retrieve the auxiliary function corresponding to the dimension
# See the appendix of the documentation for its derivation
# Input:
#   dim: dimension of the integral
#   a: parameter of this genz function; default = 5
#   u:parameter of this genz function; default = 0.5
{
    if (dim %% 4 == 1){
        funcVal <- sin(x)
    } else if (dim %% 4 == 2){
        funcVal <- -cos(x)
    } else if (dim %% 4 == 3){
        funcVal <- -sin(x)
    } else {
        funcVal <- cos(x)
    }
  return(funcVal)
}

oscillatoryIntegral <- function(dim=NA, u=0.5, a=5)
# Computes the true integral of the d-dimensional oscillatory function
# See the appendix of the documentation for its derivation
# Input:
#   dim: dimension of the integral
#   a: parameter of this genz function; default = 5
#   u:parameter of this genz function; default = 0.5
{
    # Set auxiliary function
    phi <- auxiliaryFunction

    # Compute integral
    integral <- 0
    for (n in dim:0){
      integral <- choose(dim, n) * phi((2*pi*u + n * a), dim) * (-1)^(dim - n) + integral
    }
    return(integral / a^dim)
}


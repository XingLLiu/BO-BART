contIntegral <- function(dim, a=5, u=0.5)
  
  # Description: Calculate the integral of continuous integrand Genz functions in [0,1]^dim
  # Input:
  #        dim: dimension of integrand
  #        a: parameter a in continuous integrand
  #        u: parameter u in continuous integrand, in the domain [0,1]
  # 
  #Output:
  #        integral: integral of the continuous integrand in the domain [0,1]^dim 
  
{
  if ( u < 0 || u >1 )
  {
    stop("u must be in [0,1]")
  }
  
  oneDimensionResult <- 1/a * ( 2 - exp( -a * u ) - exp( a * ( u - 1 ) ) )
  
  result <- oneDimensionResult ^ dim
  
  return(result)
  
}

copeakIntegral <- function(dim, a=5)
# Computes the true integral of the d-dimensional copeak function
# See the appendix of the documentation for its derivation
# Input:
#   dim: dimension of the integral
# Output:
#   Value of integral over [0, 1]^d
{
  # Set a value
  a <- 600/dim^3

  integral <- 0
  for (n in 1:dim){
    integral <- ( (a * n + 1) * (a * (n - 1) + 1) )^(-1) * choose((dim - 1), (n - 1)) * (-1)^(n -1) + integral
  }
  return ( integral / (a^(dim - 1) * factorial(dim)) )
}


discIntegral <- function(dim, a=5, u=0.5)
  
  # Description: Calculate the integral of discontinuous integrand Genz functions in [0,1]^dim
  # Input:
  #        dim: integer, dimension of integrand, must be greater than 1
  #        a: parameter a in continuous integrand
  #        u: parameter u in continuous integrand, in the domain [0,1]
  # 
  #Output:
  #        integral: integral of the Discontinuous integrand in the domain [0,1]^dim
  
{
  if ( u < 0 || u >1 )
  {
    stop("u must be in [0,1]")
  }
  
  if (dim < 2){
    stop("dimension must be greater than 2")
  }
  
  result <- ( ( 1/a ) * ( exp( a * u ) - 1 ) )  ^ 2  *  ( (1/a) * ( exp(a) - 1 ) ) ^ ( dim - 2 ) 
  
  return(result)
  
}

gaussianIntegral <- function(dim, a=5, u=0.5)
  
  # Description: Calculate the integral of Gaussian-Peak integrand Genz functions in [0,1]^dim
  # Input:
  #        dim: dimension of integrand
  #        a: parameter a in continuous integrand
  #        u: parameter u in continuous integrand, in the domain [0,1]
  #
  # errorFunction:
  #        Calculate the value of error function at point x
  # 
  #Output:
  #        integral: integral of the Gaussian Peak integrand in the domain [0,1]^dim 
  
{
  if ( u < 0 || u >1 )
  {
    stop("u must be in [0,1]")
  }
  
  errorFunction <- function (x) pnorm( x, mean = 0, sd = sqrt(0.5) ) - 
                                         pnorm( -x, mean = 0, sd = sqrt(0.5) )
  
  oneDimensionResult <- ( sqrt(pi) / (2 * a) ) * ( errorFunction( a * u ) + errorFunction( a - a* u ) )
  
  result <- oneDimensionResult ^ dim
  
  return(result)
  
}

auxiliaryFunction <- function(x, dim)
# Retrieve the auxiliary function corresponding to the dimension
# This itself is not a Genz function; see the appendix of the documentation for its derivation
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

oscillatoryIntegral <- function(dim, u=0.5, a=5)
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
    for (n in 0:dim){
      integral <- choose(dim, n) * phi((2*pi*u + n * a), dim) * (-1)^(dim - n) + integral
    }
    return(integral / a^dim)
}



productPeakIntegral <- function(dim, a=5, u=0.5)
  
  # Description: Calculate the integral of Product-peak integrand Genz functions in [0,1]^dim
  # Input:
  #        dim: dimension of integrand
  #        a: parameter a in continuous integrand
  #        u: parameter u in continuous integrand, in the domain [0,1]
  # 
  #Output:
  #        integral: integral of the Product Peak integrand in the domain [0,1]^dim 
  
{
  if ( u < 0 || u >1 )
  {
    stop("u must be in [0,1]")
  }
  
  oneDimensionResult <- ( -a ) * ( atan( a * ( u - 1 ) ) - atan( a * u ) )
  
  result <- oneDimensionResult ^ dim
  
  return(result)
}

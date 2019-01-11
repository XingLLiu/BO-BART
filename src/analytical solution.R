cont_integral <- function(dim, a, u)
  
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
    print ("u must be in [0,1]")
    
    break
  }
  
  one_dimension_result <- 1/a * ( 2 - exp( -a * u ) - exp( a * ( u - 1 ) ) )
  
  result <- one_dimension_result ^ dim
  
  return (result)
  
}

disc_integral <- function(dim, a, u)
  
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
    print ("u must be in [0,1]")
    
    break
  }
  
  if (dim < 2){
    
    print ("dimension must be greater than 2")
    
    break
    
  }
  
  result <- ( ( 1/a ) * ( exp( a * u ) - 1 ) )  ^ 2  *  ( (1/a) * ( exp(a) - 1 ) ) ^ ( dim - 2 ) 
  
  return (result)
  
}

gaussian_integral <- function(dim, a, u)
  
  # Description: Calculate the integral of Gaussian-Peak integrand Genz functions in [0,1]^dim
  # Input:
  #        dim: dimension of integrand
  #        a: parameter a in continuous integrand
  #        u: parameter u in continuous integrand, in the domain [0,1]
  #
  # error_function:
  #        Calculate the value of error function at point x
  # 
  #Output:
  #        integral: integral of the Gaussian Peak integrand in the domain [0,1]^dim 
  
{
  if ( u < 0 || u >1 )
  {
    print ("u must be in [0,1]")
    
    break
  }
  
  error_function <- function (x) pnorm( x, mean = 0, sd = sqrt(0.5) ) - 
                                         pnorm( -x, mean = 0, sd = sqrt(0.5) )
  
  one_dimension_result <- ( sqrt(pi) / (2 * a) ) * ( error_function( a * u ) + error_function( a - a* u ) )
  
  result <- one_dimension_result ^ dim
  
  return (result)
  
}

product_peak_integral <- function(dim, a, u)
  
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
    print ("u must be in [0,1]")
    
    break
  }
  
  one_dimension_result <- ( -a ) * ( atan( a * ( u - 1 ) ) - atan( a * u ) )
  
  result <- one_dimension_result ^ dim
  
  return (result)
}
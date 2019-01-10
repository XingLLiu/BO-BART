analytical_sol <- function (FUN, dim){
  
  #Analytical Solution of Genz, for cont, gaussian, prpeak and disc only 
  if (FUN == "disc"){
    
    if (dim == 1){
      
      print ("disc undefined for dimension = 1")
      
      break
      
    }
    
    f <- function (x) exp(5 * x)
    
    result <- (( adaptIntegrate(f,lowerLimit = rep(0, dim),upperLimit = rep(0.5, dim))[[1]] ) ^ (2)) *
                ( adaptIntegrate(f,lowerLimit = rep(0, dim),upperLimit = rep(1, dim))[[1]] ) ^ (dim - 2)
    
    
  }else{
    
    result <- (adaptIntegrate(FUN,lowerLimit = rep(0,3),upperLimit = rep(1,3))[[1]]) ^ (dim)
    
  }
  
  return (result)
  
  
}

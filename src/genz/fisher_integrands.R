h_func <- function(x) 
{
  return (1/(1 + exp(-80*x)))
}

gf_func <- function(x, f) 
{
  result <- rep(0, nrow(x))
  result[abs(x)<=1] <- exp(-1/(1 - abs(x[abs(x)<=1])^2) +
                             cos(f*pi*abs(x[abs(x)<=1])))
  return (result)
}

create_fisher_function <- function(C, R, H, F, P, dim) {
  fisher_function <- function(x) {
    integrand <- 1
    for (i in 1:dim) {
      integrand <- integrand*(
        H[i]*gf_func((matrix(x[,i]) - C[i])/R[i], F[i]) + (-1)^P[i]*(
          0.5 - h_func(matrix(x[,i]) - C[i])
        )
      )
    }
    return (integrand)
  }
  return(fisher_function)
}




library(matrixStats)  
library(MASS)
cont <- function(xx, u=rep(0.5, 1, length(xx)), a=rep(5, 1, length(xx)))
{
  ##########################################################################
  #
  # CONTINUOUS INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u = c(u1, u2, ..., ud) (optional), with default value
  #     c(0.5, 0.5, ..., 0.5)
  # a = c(a1, a2, ..., ad) (optional), with default value c(5, 5, ..., 5)
  #
  ##########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    
    xx <- matrix(xx) 
    
  }
  
  u=rep( 0.5, 1, ncol(xx) )
  
  a=rep(5, ncol(xx)) 
  
  sum <- abs(xx-u) %*% a
  
  y <- exp(-sum)
  
  
  return(y)
}


copeak <- function(xx, u=rep(0.5, 1, ncol(xx)), a=rep(5, ncol(xx)) )
{
  ##########################################################################
  #
  # CORNER PEAK INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u = c(u1, u2, ..., ud) (optional), with default value
  #     c(0.5, 0.5, ..., 0.5)
  # a = c(a1, a2, ..., ad) (optional), with default value c(5, 5, ..., 5)
  #
  #########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx) 
  }
  
  u=rep( 0.5, 1, ncol(xx) )
  a=rep(5, ncol(xx)) 
  
  
  d <- ncol(xx) # change to ncol
  
  sum <- xx %*% a
  
  y <- (1 + sum)^(-d-1)
  
  return(y)
}


disc <- function(xx, u=rep(0.5, 1, length(xx)), a=rep(5, 1, length(xx)))
{
  ##########################################################################
  #
  # DISCONTINUOUS INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u = c(u1, u2, ..., ud) (optional), with default value
  #     c(0.5, 0.5, ..., 0.5)
  # a = c(a1, a2, ..., ad) (optional), with default value c(5, 5, ..., 5)
  #
  #########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    
    xx <- matrix(xx) 
    
  }
  
  # Function only defined for dimension >= 2
  if (ncol(xx) < 2) stop("incorrect dimension. Discrete Genz function only defined for dimension >= 2") 

  x1 <- xx[ ,1]
  x2 <- xx[ ,2]
  u1 <- u[1]
  u2 <- u[2]
  
  u=rep( 0.5, 1, ncol(xx) )
  
  a=rep(5, ncol(xx))
  
  xx [which ( x1 > u1 | x2 > u2), ] <- 0

  sum <- xx %*% a
  
  y <- exp(sum)
  
  y[which (y == 1)] <- 0

  
  return(y)
  
}


gaussian <- function(xx, u=rep(0.5, 1, length(xx)), a=rep(5, 1, length(xx)))
{
  ##########################################################################
  #
  # GAUSSIAN PEAK INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u  = c(u1, u2, ..., ud) (optional), with default value
  #      c(0.5, 0.5, ..., 0.5)
  # a  = c(a1, a2, ..., ad) (optional), with default value c(5, 5, ..., 5)
  #
  ##########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    
    xx <- matrix(xx) 
    
  }
  
  u=rep( 0.5, 1, ncol(xx) )
  
  a=rep(5, ncol(xx)) 
  
  sum <- (xx - u)^2 %*% a^2
  
  y <- exp(-sum)
  
  return(y)
  
}

oscil <- function(xx, u=rep(0.5, 1, length(xx)), a=rep(5, 1, length(xx)))
{
  ##########################################################################
  #
  # OSCILLATORY INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u  = c(u1, u2, ..., ud) (optional), with default value
  #      c(0.5, 0.5, ..., 0.5)
  # a  = c(a1, a2, ..., ad) (optional), with default value c(5, 5, ..., 5)
  #
  ##########################################################################
  
  if (is.matrix(xx) == FALSE) { 
    
    xx <- matrix(xx) 
    
  }
  
  u=rep( 0.5, 1, ncol(xx) )
  
  a=rep(5, ncol(xx)) 
  
  sum <- xx %*% a
  
  y <- cos(2 * pi * u[1] + sum)
  
  return(y)
  
}



prpeak <- function(xx, u=rep(0.5, 1, length(xx)), a=rep(5, 1, length(xx)))
{
  ##########################################################################
  #
  # PRODUC PEAK INTEGRAND FAMILY
  #
  # Authors: Sonja Surjanovic, Simon Fraser University
  #          Derek Bingham, Simon Fraser University
  # Questions/Comments: Please email Derek Bingham at dbingham@stat.sfu.ca.
  #
  # Copyright 2013. Derek Bingham, Simon Fraser University.
  #
  # THERE IS NO WARRANTY, EXPRESS OR IMPLIED. WE DO NOT ASSUME ANY LIABILITY
  # FOR THE USE OF THIS SOFTWARE.  If software is modified to produce
  # derivative works, such modified software should be clearly marked.
  # Additionally, this program is free software; you can redistribute it 
  # and/or modify it under the terms of the GNU General Public License as 
  # published by the Free Software Foundation; version 2.0 of the License. 
  # Accordingly, this program is distributed in the hope that it will be 
  # useful, but WITHOUT ANY WARRANTY; without even the implied warranty 
  # of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
  # General Public License for more details.
  #
  # For function details and reference information, see:
  # http://www.sfu.ca/~ssurjano/
  #
  ##########################################################################
  #
  # INPUTS:
  #
  # xx = c(x1, x2, ..., xd)
  # u  = c(u1, u2, ..., ud) (optional), with default value
  #      c(0.5, 0.5, ..., 0.5)
  # a  = c(a1, a2, ..., ad) (optional), with default value c(5, 5, ..., 5)
  #
  ##########################################################################

  if (is.matrix(xx) == FALSE) { 
    
    xx <- matrix(xx) 
    
  }
  
  u=rep( 0.5, 1, ncol(xx) )
  
  a=rep(5, ncol(xx)) 
  
  sum <- a^( -2 ) + ( xx - u ) ^ 2
  
  y <- rowProds(1/sum)
  
  return(y)
}




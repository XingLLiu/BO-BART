mixtureGenz <- function(xc){

    gen0 <- copeak
    gen1 <- disc
    if (is.matrix(xc) == FALSE) { 
    xc <- matrix(xc, nrow = 1) 
    }

    # Set Genz functions
    gen0 <- prpeak
    gen1 <- gaussian
    dim <- ncol(xc)

    c <- xc[ , dim]
    xx <- xc[ , 1:(dim - 1)]

    if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1)  
    }

    # Compute Genz function values
    indices_0 <- which(c == 0)
    indices_1 <- which(c == 1)

    y0 <- as.numeric(gen0( xx[indices_0, ] ))
    y1 <- as.numeric(gen1( xx[indices_1, ] ))


    y <- rep(0, nrow(xx))
    y[indices_0] <- y0
    y[indices_1] <- y1

    return(y)

}

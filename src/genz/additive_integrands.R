additive_gaussian <- function(xx, a=NA){

    if (is.matrix(xx) == FALSE) { 
        xx <- matrix(xx, ncol = 1) 
    }

    dim <- ncol(xx)
    u <- 0.5
    if ( any(is.na(a)) ){
        a <- 100/dim^2
    }

    if (length(a) == 1){
        a <- rep(a, dim)
    }

    sum <- 0
    for (j in 1:dim){
        sum <- sum + exp(-a[j]^2 * (xx[, j] - u)^2)
    }

    return(sum)
}
# ================================================================
# Rejection method to sample from Genz functions
# ================================================================
# # compute M
# M <- (optimize(phi, c(0, 10), maximum = TRUE, x0 = x0, gamma = gamma))$objective




rejection_samp <- function(N = NA, lower_lim = NA, upper_lim = NA, M = NA, f_func = NA){
    ##########################################################################
    # Rejection sampling.

    # Input :
    #         N = [int] number of realizations to sample.
    #         lower_lim = [double] lower limit of the support of f_func.
    #         upper_lim = [double] upper limit of the support of f_func.
    #         M = [float] sup(f^*) over the support.
    #         f_func = [function] Genz function to sample from.

    # Output:
    #         x_sample = [vector] sampled realizations from f_func.
    ##########################################################################

    # initialize
    x_sample <- matrix(NA, nrow = N, ncol = 1)
    n <- 0   # looping index
    n_total <- 0   # no. of total trials
    while (n < N){
        # generate realizations from Unif(0, 1)
        u <- runif(2, 0, 1)
        # generate realization from g
        y <- (upper_lim - lower_lim) * u[1] + lower_lim
        # accept with probability f/(M * g^*)
        phi_val <- f_func(y)
        if (u[2] * M < phi_val){
            n <- n + 1
            x_sample[n] <- y
        }
        # update no. of total trials
        n_total <- n_total + 1
    }
    cat('Total no. of trials =', n_total, '\n')
    cat('Acceptance rate =', round(N/n_total, 4), '\n')
    return(x_sample)
}





# ================================================================
# Example (1d product peak)
# ================================================================
library(matrixStats)
library(lhs)
# Function to be sampled from
f_func <- function(xx){

    if (is.matrix(xx) == FALSE) { 
    xx <- matrix(xx, nrow = 1)  
    }

    dim <- ncol(xx)
    u <- rep(0.5, dim)
    a <- rep(600/dim^3, dim)

    sum <- a^(-2) + (xx - u)^2
    y <- rowProds(1/sum)

    if (xx <= 1 & xx >= 0){
        return(y)
    }
    else{
        return(0)
    }
}

# Rejection sampling
# N <- 10^4
# M <- 360000
# sample_rej <- rejection_samp(N, lower_lim = 0, upper_lim = 1, M = M, prpeak)

# # Random LHS sampling
# xx <- randomLHS(N, 1)
# sample_rand <- prpeak(matrix(sort(xx), length(xx), 1))

# x_vec <- 0:N/N
# hist(sample_rej, probability = TRUE, breaks = seq(0, 1, l = 500), col = "slategrey",
#     main = "Different sampling methods")
# lines(x_vec, 1/3600 * prpeak(matrix(x_vec, length(x_vec), 1)), col = "red")
# points(sort(xx), 1/3600 * sample_rand, col = "blue", cex = 0.2)
# legend("topright", legend = c("rejection sampling", "true function", "random sampling"),
#        col=c("black", "red", "blue"), cex = 0.8, lty = c(1, 1, 1))

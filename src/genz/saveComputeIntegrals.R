# We compute the exact integrals for all the genz functions for dimensions
# 1,2,3,5,10,20 and then store them as CSV format
# we use the notation from https://www.sfu.ca/~ssurjano/disc.html
# thus alpha corresponds to the a's, and beta corresponds to the b's
# we fix a = [5,5,5,5,5...]
# and    u = [0.5,....]
# Hence, it is sufficient to use scalars a = 5 and u = 0.5.
# we integrate from 0 to 1

source("./genz/analyticalIntegrals.R")

numGenz <- 6
dimensions <- c(1, 2, 3, 5, 10, 20)
integrals <- matrix(rep(NA, length(dimensions) * numGenz), nrow = numGenz)

for (k in 1:length(dimensions)){

    # Compute integrals
    dim <- dimensions[k]
    integrals[1, k] <- contIntegral(dim)
    integrals[2, k] <- copeakIntegral(dim)
    integrals[4, k] <- gaussianIntegral(dim)
    integrals[5, k] <- oscillatoryIntegral(dim)
    integrals[6, k] <- productPeakIntegral(dim)

    # Discontinuous integrand only defined for dim >= 2
    if (dim > 1){
        integrals[3, k] <- discIntegral(dim)
    }
}

# Write results to file
write.table(integrals, file = "./genz/integrals.csv", sep=",", row.names=FALSE, col.names=FALSE)

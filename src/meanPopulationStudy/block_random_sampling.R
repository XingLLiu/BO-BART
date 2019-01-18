computeMean <- function(trainX, trainY, candidateX, candidateY, num_iterations)
{
    # stratify the population
    maleCandidateY <- candidateY[candidateX$SEX == 1]
    femaleCandidateY <- candidateY[candidateX$SEX == 2]

    maleRatio <- sum(trainX$SEX == 1) / nrow(trainX)

    numMaleCandidate <- floor(500 * maleRatio)
    numFemaleCandidate <- 500 - numMaleCandidate

    BRmean <- c()
    BRstandardDeviation <- c()

    for (i in 1:numMaleCandidate) {
        BRmean[i] <- mean(c(trainY, maleCandidateY[1:i]))
        BRstandardDeviation[i] <- sqrt( var(c(trainY, maleCandidateY[1:i])) )
    }

    for (i in 1:numFemaleCandidate) {
        BRmean[i+numMaleCandidate] <- mean(c(trainY, maleCandidateY[1:numMaleCandidate], femaleCandidateY[1:i]))
        BRstandardDeviation[i+numMaleCandidate] <- sqrt( var(c(trainY, maleCandidateY[1:numMaleCandidate], femaleCandidateY[1:i])) )
    }
    
    return(list("BRmean" = BRmean, "BRstandardDeviation" = BRstandardDeviation))

}
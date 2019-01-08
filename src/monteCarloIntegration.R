monteCarloIntegrationUniform <- function(FUN, numSamples, dim)
# Monte Carlo Integration
# input:
# 		FUN: function f in the integral
#		numSamples: number of samples to take
# output:
#		prediction: list containing meanValue2 and standardDeviation of the integral
{
	meanValue2 <- rep(0, numSamples)
	standardDeviation2 <- rep(0, numSamples)
	
	for (i in 1:numSamples) {
	  CandidateSet <- matrix(runif(i*dim), ncol = dim)
	  integration <- mean(FUN(CandidateSet))
	  meanValue2[i] <- integration
	  standardDeviation2[i] <- sqrt(var(meanValue2))
	}
	return(list("meanValueMonteCarlo" = meanValue2, "standardDeviationMonteCarlo" = standardDeviation2))
}

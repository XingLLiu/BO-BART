monteCarloIntegrationUniform <- function(FUN, numSamples, dim)
# Monte Carlo Integration
# input:
# 		FUN: function f in the integral
#		numSamples: number of samples to take
# output:
#		prediction: list containing meanValue2 and standardDeviation of the integral
{
	meanValueMonteCarlo <- rep(0, numSamples)
	standardDeviationMonteCarlo <- rep(0, numSamples)
	
	for (i in 1:numSamples) {
	  CandidateSet <- matrix(runif(i*dim), ncol = dim)
	  ######## Mixed Genz ########
	  CandidateSet <- cbind(randomLHS(i, (dim - 1)), sample(c(0,1), i, replace = TRUE))
	  ############################	

	  functionSamples <- FUN(CandidateSet)
	  meanValueMonteCarlo[i] <- mean(functionSamples)
	  standardDeviationMonteCarlo[i] <- sqrt(var(functionSamples))
	}
	return(list("meanValueMonteCarlo" = meanValueMonteCarlo, "standardDeviationMonteCarlo" = standardDeviationMonteCarlo))
}

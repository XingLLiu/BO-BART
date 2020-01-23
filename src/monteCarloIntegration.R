monteCarloIntegrationUniform <- function(FUN, numSamples, dim, measure)
  #'Crude Monte Carlo Approximation
  #' 
  #'@description The function approximates the integral of interest using curde monte carlo
  #'There is no sequential sampling/adaptive Bayesian quadrature as there are no posterior samples
  #' 
  #'@param FUN Function; The function to be integrated 
  #'@param numSamples Integer; The number of samples used in calculate mean
  #'@param dim Integer; The dimension of the input X
  #'
  #'@return List; A list containing meanValue (appximation) and the variance of crude monte Carlo
{
	meanValueMonteCarlo <- rep(0, numSamples)
	standardDeviationMonteCarlo <- rep(0, numSamples)
	
	for (i in 1:numSamples) {
	  
	  if (measure == "uniform") {
	    CandidateSet <- as.matrix(replicate(dim, runif(i)))
	  } else if (measure == "gaussian") {
	    CandidateSet <- as.matrix(replicate(dim, rtnorm(i, mean = 0.5, lower=0, upper=1)))
	  }

	  functionSamples <- FUN(CandidateSet)
	  meanValueMonteCarlo[i] <- mean(functionSamples)
	  standardDeviationMonteCarlo[i] <- sqrt(var(functionSamples))
	}
	return(list("meanValueMonteCarlo" = meanValueMonteCarlo, "standardDeviationMonteCarlo" = standardDeviationMonteCarlo))
}

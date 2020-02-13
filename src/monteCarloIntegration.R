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
	
	for (i in 50:(50+numSamples)) {
	  
	  if (measure == "uniform") {
	    CandidateSet <- replicate(dim, runif(i, 0, 1))
	  } else if (measure == "gaussian") {
	    CandidateSet <- replicate(dim, rtnorm(i, mean=0.5, lower=0, upper=1))
	  }

	  functionSamples <- FUN(CandidateSet)
	  meanValueMonteCarlo[i-49] <- mean(functionSamples)
	  standardDeviationMonteCarlo[i-49] <- sqrt(var(functionSamples))
	}
	return(list("meanValueMonteCarlo" = meanValueMonteCarlo, "standardDeviationMonteCarlo" = standardDeviationMonteCarlo))
}

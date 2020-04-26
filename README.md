# BO-BART

[![Documentation Status](https://readthedocs.org/projects/bart-bq/badge/?version=latest)](https://bart-bq.readthedocs.io/en/latest/?badge=latest)


Bayesian Optimization using Bayesian Additive Regression Trees (under review)

- Zhu H, Liu X, Briol F-X, Kang R, Shen Z, and Flaxman S

Project Documentation: https://bart-bq.readthedocs.io/en/latest/

## Code directory ##

    .
    |
    └── README.md
    ├───├data: Data used for survey design experiment
    ├───├Figures: Convergence plots for the benchmark testing
    ├───├python: Python codes for hyper-parameter tuning for GP
    ├───├results: Results for Genz testing. This is where we store the results generated by our integral approximation functions, as well as the analytical integrals of the benchmark testing functions.
    ├───├src
    	├── genz: Source files used to compute and store the analytic integrals
    	├── meanPopulationStudy: Source files used to conduct survey design experiment
        ├── BARTBQ.R: Implementation of BART-Int
        ├── GPBQ.R: Implementation of Bayesian Quadrature with Gaussian processes interface
        ├── GPBQRunTime.R: Main class of Bayesian Quadrature with Gaussian processes, with heuristic used to compute median bandwidth
        ├── mixedGenzMain.R: Main class to do BART-Int, GPBQ and Monte Carlo integrations with a mixture of two Genz families (**Archived**)
        ├── monteCarloIntegration.R: Main class of crude Monte Carlo integration
        ├── optimise_gp.R: Source file used to optimise the length scale using Pytorch
        ├── setup.R: Source file used to download all dependencies (**Archived**)
    ├── compute_CV.R: Main class for multiple runs of the entire pipeline
    ├── draw_BART_posterior.R: Main class to draw the posterior densities of BART
    ├── draw_GP_posterior.R: Main class to draw the posterior densities of GP
    ├── gp_test.R: Temporary (**Archived**) 
    ├── integrationMain.R: Main class to do BART-Int, GPBQ and Monte Carlo integrations. Tweak your genz functions and parameters here
    ├── plot_posterior_example.R: Main class to plot a toy example of the posteriors of BART-Int and GP
    ├── plot_step_function_posterior.R: Main class to plot the posteriors of BART-Int and GP with the step function
    ├── test_gaussian_prior.R: Main class to do BART-Int, GPBQ and Monte Carlo integrations with a truncated Gaussian prior
    ├── test_step_function.R: Main class to run the integration pipeline with a step function (**Archived**)
	    



## Dependencies

```r
    yaml
    MASS
    cubature
    lhs
    data.tree
    dbarts
    matrixStats
    mvtnorm
    doParallel
    kernlab
    msm
    dbarts_0.9-8
    caret
```

## To run the genz integrals approximations:

1) Install all the necessary packages

```r
install.packages(c("yaml", "MASS", "cubature", "lhs", "data.tree", "matrixStats", "mvtnorm", "doParallel", "kernlab", "msm", "caret"))

# an old version of dbarts
packageurl <- "https://cran.r-project.org/src/contrib/Archive/dbarts/dbarts_0.9-8.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
```

2) To reproduce the benchmark tests, run `integrationMain.R` with customized inputs. There are 7 arguments in total, of which the last two are optional. The last argument should only be specified when the step function is used (`genz_function_number = 7`), and is set to `1` if not specified. For example:
```
Rscript integrationMain.R dimension num_iterations genz_function_number kernel_name sequential_flag (measure) (number_of_jumps_for_step_function)

```
where `genz_function_number` follows the indexing in this [documentation](https://www.sfu.ca/~ssurjano/integration.html) for the Genz families. The results will be stored in `results`, where you can find the `.csv` and `.RData` files containing the numerical values and the automatically generated graphs.

3) If you want to tune the GP integral, you can also run

```
Rscript GPRunTime.R dimension num_iterations genz_function_number initial_training_set_size
```

Results will also be stored in `results` and `Figures`.

## To test the design with real-life data with provided

1) Install the dependencies in `R`. Make sure you are using **R 3.5.0** or higher.

2) Run

```
Rscript src/meanPopulationStudy/poptMean.R num_new_surveys

```

This will generate and store the results in `results/populationStudy` and `Figures/populationStudy`, where you can find the `.csv` and `.RData` files containing the numerical values and the automatically generated graphs.


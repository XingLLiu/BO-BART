# BO-BART

[![Documentation Status](https://readthedocs.org/projects/bart-bq/badge/?version=latest)](https://bart-bq.readthedocs.io/en/latest/?badge=latest)


Bayesian Optimization using Bayesian Additive Regression Trees (under review)

- Kang R, Liu X, Shen Z, Zhu H and Flaxman S

Project Documentation: https://bart-bq.readthedocs.io/en/latest/

## Code directory ##

    .
    |
    └── README.md
    ├───├data: contains the data used for survey design experiment
    ├───├notes: notes
    ├───├report: report tex files
    ├───├results: results for Genz testing. This is where we store the results generated by our integral approximation functions.
    ├───├src
    	├── genz: stores the source files used to compute the analytic integrals and store them
    	├── meanPopulationStudy: source files used to conduct survey design experiment
    	├── packages: source files to load dependencies (libraries)
        ├── BARTBQ.R: Implementation of Bayesian Quadrature with BART
        ├── GPBQ.R: Implementation of Bayesian Quadrature with Gaussian processes interface
        ├── GPBQRunTime.R: Main class of Bayesian Quadrature with Gaussian processes, with heuristic used to compute median bandwidth
        ├── monteCarloIntegration.R: Main class of crude Monte Carlo integration
        ├── integrationMain.R: Main class to do BART, GP and Monte Carlo integrations. Tweak your genz functions and parameters here
	    



## Dependencies

```r
    MASS
    cubature
    lhs
    data.tree
    dbarts
    matrixStats
    mvtnorm
    msm
```

## To run the genz integrals approximations:

1) Install all the necessary packages

```
install.packages(c("yaml", "MASS", "cubature", "lhs", "data.tree", "matrixStats", "mvtnorm", "doParallel", "kernlab", "msm"))
```
and also 
```
packageurl <- "https://cran.r-project.org/src/contrib/Archive/dbarts/dbarts_0.9-8.tar.gz"
install.packages(packageurl, repos=NULL, type="source")
```

2) To simulate, run the test scripts with customized inputs. There are 6 inputs in total; the last input is optional and only works for the step function (`genz_function_number` = 7). For example:
```
Rscript integrationMain.R dimension num_iterations genz_function_number kernel_name sequential_flag (number_of_jumps_for_step_function)

```
where `genz_function_number` is following the indexing in https://www.sfu.ca/~ssurjano/integration.html. 

This will generate the results in `../results`, where you can take a look at the `.csv` files or automatically generated graphs.

4) If you want to tune the GP integral, you can also run

```
Rscript GPRunTime.R dimension num_iterations genz_function_number initial_training_set_size
```

Results will similarly be in `../results`.

## To test the design in real-life data with provided dataset

1) Install the dependencies in `R`. Make sure you are using **R 3.5.0** or higher.

2) In the terminal, `cd` to `meanPopulationStudy`

3) Run

```
Rscript poptMean.R num_iterations

```

This will generate the results in `meanPopulationStudy` where you can take a look at the `.csv` files or automatically generated graphs.


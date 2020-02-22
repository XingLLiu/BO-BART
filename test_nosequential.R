source("src/genz_pipeline_nosequential.R")

# arguments to be passed to the testing pipeline.
# similar to the arguments in integrationMain.R.
# epochs = number of design points
# num_cv_total = number of runs
args <- list(dim = 20, num_iterations = 1, whichGenz = 9, whichKernel = "matern32", sequential = 0,
             measure = "uniform", epochs = seq(22, 23, by = 10), num_runs = 1)

genz_pipeline_nosequential(args)

plot_args <- list(dims_list = c(20), genz_list = c(9), sequential_list = c("NoSequential"),
                    measure = "uniform", num_cv = 1)
plot_results(plot_args)    
# the current setting in BOBART is correct. See the results for genz1, which seems to be more
# reasonable      

source("src/genz_pipeline_nosequential.R")

# arguments to be passed to the testing pipeline.
# similar to the arguments in integrationMain.R.
# epochs = number of design points
# num_cv_total = number of runs
args <- list(dim = 20, num_iterations = 1, whichGenz = 9, whichKernel = "matern32", sequential = 0,
             measure = "uniform", epochs = seq(25, 100, by = 5), num_runs = 20)

genz_pipeline_nosequential(args)
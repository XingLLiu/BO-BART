library(here)
library(yaml)

config <- read_yaml(file.path(here("config"), "config_local.yml"))
root_dir <- config$root_dir
libraries <- read.csv(file.path(root_dir, "src/packages/requirements.txt"))

for (lib in libraries)
{
    install.packages(lib, file.path(root_dir, "lib"))
}

packageurl <- "https://cran.r-project.org/src/contrib/Archive/dbarts/dbarts_0.9-8.tar.gz"
install.packages(packageurl, repos=NULL, type="source")

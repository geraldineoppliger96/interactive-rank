source("sim_dat.R")
source("tests.R")
source("experiment_wrap.R")

packages = c("magrittr", "reshape2", "ggplot2", "randomForest", "quantregForest",
             "ash", "caTools", "tibble", "dplyr", "tidyr", "doParallel", "foreach",
             "energy", "MASS")
for (one_pack in packages){
  if (!require(one_pack, character.only = TRUE))
  {
    install.packages(one_pack, dep=TRUE, repos = "http://cran.us.r-project.org")
    if(!require(one_pack,character.only = TRUE)) stop("Package not found")
  }
  library(one_pack, character.only = TRUE)
}


cl <- makeCluster(detectCores())
registerDoParallel(cl)


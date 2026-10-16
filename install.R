pk <- c("yaml", "ggplot2", "lme4", "pROC", "MASS", "testthat")
miss <- setdiff(pk, rownames(installed.packages()))
if (length(miss)) install.packages(miss, repos = "https://cloud.r-project.org")

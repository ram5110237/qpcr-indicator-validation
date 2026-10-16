library(testthat)
test_dir(file.path(dirname(normalizePath(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)))), "testthat"),
         reporter = "summary", stop_on_failure = TRUE)

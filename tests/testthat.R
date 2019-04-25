library(testthat)
require("oro.nifti")
library(pTFCE)

#source("tests/test.help.R")

test_check("pTFCE", load_helpers=T, reporter = JunitReporter$new(file = "junit_result.xml"))

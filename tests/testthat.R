library(testthat)
require("oro.nifti")
library(pTFCE)

source("tests/test.help.R")

test_check("pTFCE", reporter = "full", load_helpers=T)

#test cases for pTFCE
context("FWER")
#test_file("../report.log")

test_that("FWER p to Z", {

  expect_equal(fwe.p2z(100, 0.05), 4.046402, tolerance=1e-6)
  expect_equal(fwe.p2z(10, 0.01), 3.843033, tolerance=1e-6)
  expect_equal(fwe.p2z(100000, 0.00001), 7.037125, tolerance=1e-6)
})

test_that("FWER Z to p", {

  expect_equal(fwe.z2p(100, 2.3), 1, tolerance=1e-6)
  expect_equal(fwe.z2p(100, 5), 0.001045908, tolerance=1e-6)
  expect_equal(fwe.z2p(10000, 6), 0.0006233478, tolerance=1e-6)
})

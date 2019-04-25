#test cases for pTFCE
require("oro.nifti")
context("random data")
#test_file("../report.log")

test_that("Volume with random data", {

  expected=readNIfTI(system.file("extdata", "test01A_out.nii.gz", package="pTFCE"))

  ptfce=ptfce(img = readNIfTI(system.file("extdata", "test01A_in.nii.gz", package="pTFCE")),
        V=20*20*20,
        Rd = 20, #magic value
        mask = readNIfTI(system.file("extdata", "test01A_mask.nii.gz", package="pTFCE")),
        Nh=50,
        verbose = F
  )$p

  expect_equal(mean(ptfce), mean(expected))
  expect_equal(sd(ptfce), sd(expected))
  expect_equal(ptfce[10,10,10], expected[10,10,10])
})

test_that("Volume with random data, smoothness estimated", {

  expected=readNIfTI(system.file("extdata", "test01B_out.nii.gz", package="pTFCE"))

  ptfce=ptfce(img = readNIfTI(system.file("extdata", "test01B_in.nii.gz", package="pTFCE")),
              mask = readNIfTI(system.file("extdata", "test01B_mask.nii.gz", package="pTFCE")),
              Nh=50,
              verbose = F
  )$p

  #writeNIfTI(ptfce,"../../data/test01B_out")

  expect_equal(mean(ptfce), mean(expected))
  expect_equal(sd(ptfce), sd(expected))
  expect_equal(ptfce[10,10,10], expected[10,10,10])
})




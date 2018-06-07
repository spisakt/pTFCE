#test cases for pTFCE
require("oro.nifti")
context("random data")
#test_file("../report.log")

test_that("Volume with random data", {

  expected=readNIfTI("../../data/test01A_out.nii.gz")

  ptfce=ptfce(img = readNIfTI("../../data/test01A_in.nii.gz"),
        V=20*20*20,
        Rd = 20, #magic value
        mask = readNIfTI("../../data/test01A_mask.nii.gz", verbose = F)
  )$pTFCE

  expect_equal(mean(ptfce), mean(expected))
  expect_equal(sd(ptfce), sd(expected))
  expect_equal(ptfce[10,10,10], expected[10,10,10])
})

test_that("Volume with random data, smoothness estimated", {

  expected=readNIfTI("../../data/test01B_out.nii.gz")

  ptfce=ptfce(img = readNIfTI("../../data/test01B_in.nii.gz"),
              mask = readNIfTI("../../data/test01B_mask.nii.gz"),
              verbose = F
  )$pTFCE

  #writeNIfTI(ptfce,"../../data/test01B_out.nii.gz")

  expect_equal(mean(ptfce), mean(expected))
  expect_equal(sd(ptfce), sd(expected))
  expect_equal(ptfce[10,10,10], expected[10,10,10])
})




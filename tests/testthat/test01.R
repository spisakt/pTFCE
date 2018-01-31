#test cases for pTFCE
require("oro.nifti")
context("value")
#test_file("../report.log")

test_that("Volume with random data", {

  expected=readNIfTI("../../data/test01A_out.nii.gz")

  ptfce=ptfce(img = readNIfTI("../../data/test01A_in.nii.gz"),
        V=20*20*20,
        Rd = 20, #magic value
        mask = readNIfTI("../../data/test01A_mask.nii.gz")
  )$pTFCE

  expect_equal(mean(ptfce), mean(expected))
  expect_equal(sd(ptfce), sd(expected))
  expect_equal(sd(ptfce[10,10,10]), sd(expected[10,10,10]))
})



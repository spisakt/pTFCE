# helper functions for testing
require("oro.nifti")
simulateVol=function(dim=20, data=rnorm(dim*dim*dim))
{
  vol=array(data, dim=c(20,20,20))
  return(vol)
}

createTestVol=function(num=1, tag="A")
{
  set.seed(1986)
  if (num==1)
  {
    dim=20
    vol=simulateVol(dim)
    mask=simulateVol(dim, 1)
    ret=ptfce(img = vol, V=20*20*20, Rd = 20, mask = mask, length.out = 50)
    writeNIfTI(vol, "data/test01A_in")
    writeNIfTI(ret$pTFCE, "data/test01A_out")
    writeNIfTI(mask, "data/test01A_mask")

    write(paste(date(), "Test", num),file="data/generated.test.data.log",append=TRUE)
  }
  else
  {
    warning("No such test case!")
  }

}


loadTestVol=function(num, tag="A", inout="in")
{
  readNIfTI(paste("data/test", sprintf("%.2d", num), tag, "_", inout, ".nii.gz" , sep=""))
}

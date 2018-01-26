#main function and helper functions for pTFCE R package

# imports:
devtools::use_package("oro.nifti")
devtools::use_package("mmand")

#main function

#' pTFCE main function.
#'
#' The threshold-free cluster enhancement (TFCE) approach integrates cluster information into voxel-wise statistical inference to enhance detectability of neuroimaging signal. Despite the significantly increased sensitivity, the application of TFCE is limited by several factors: (i) generalization to data structures, like brain network connectivity data is not trivial, (ii) TFCE values are in an arbitrary unit, therefore, P-values can only be obtained by a computationally demanding permutation-test.
#' Here, we introduce a probabilistic approach for TFCE (pTFCE), that gives a simple general framework for topology-based belief boosting.
#' The core of pTFCE is a conditional probability, calculated based on Bayes' rule, from the probability of voxel intensity and the threshold-wise likelihood function of the measured cluster size. We provide an estimation of these distributions based on Gaussian Random Field (GRF) theory. The conditional probabilities are then aggregated across cluster-forming thresholds by a novel incremental aggregation method. Our approach is validated on simulated and real fMRI data.
#' pTFCE is shown to be more robust to various ground truth shapes and provides a stricter control over cluster "leaking" than TFCE and, in the most realistic cases, further improves its sensitivity. Correction for multiple comparison can be trivially performed on the enhanced P-values, without the need for permutation testing, thus pTFCE is well-suitable for the improvement of statistical inference in any neuroimaging workflow.
#'
#' @param img Nifti image to enhance ("nifti" class from "oro.nifti" package)
#' @param Rd Resel count (as output by FSL smoothest)
#' @param V Number of voixels in mask
#' @param mask Mask
#' @param length.out Number of thresholds
#' @param logpmin min threshold
#' @param logpmax max threshold
#'
#' @return TFCE object
#' @export
#'
#' @examples Z=readNIfTI("Zmap.nii.gz");
#' MASK=readNIfTI("mask.nii.gz");
#' smooth=read.table("smoothness.txt");
#' ptfce(Z, V=smooth[2,2], Rd = smooth[1,2]*smooth[2,2], mask = MASK, length.out = 100);
ptfce=function(img, Rd, V, mask, length.out=50,   logpmin=0, logpmax=-log(pnorm(max(img), lower.tail = F)) )
{
  logp.thres=seq(logpmin, logpmax, length.out=length.out)
  dh=logp.thres[2]-logp.thres[1]
  p.thres=exp(-logp.thres)
  threshs=qnorm(p.thres, lower.tail = F)
  threshs[1]=-9999999999 # exp(-Inf):=0 #TODO: machine minimum
  ndh = length(threshs)
  CLUST=array(NA, dim=c(dim(img), length(threshs)))
  PVC=array(1, dim=c(dim(img), length(threshs)))
  for (hi in 1:length(threshs))
  {
    h=threshs[hi]
    thr=img
    thr[img>h]=1
    thr[img<=h]=0
    ccc=mmand::components(thr, kernel = mmand::shapeKernel(3, dim = 3, type = "diamond"))
    t=table(ccc)
    sizes=t[ccc]
    sizes=sizes*as.numeric(mask)

    CLUST[,,,hi]=sizes

    for (size in sort(unique(sizes)) )
    {
      if (size==0)
        next;
      PVC[,,,hi][sizes==size]=pvox.clust(V, Rd, size, h)
    }
    # plot if you want
    # orthographic(nifti(CLUST[,,, hi], datatype=16) )
  }
  # calculate pTFCE
  pTFCE=array(apply(PVC, c(1,2,3), function(x){exp( -aggregate.logpvals(-log(x), dh) )}), dim=dim(img))
  pTFCE[pTFCE==0]=4.940656e-324 #underflow #TODO: handle this better
  #TODO: return object
  return(list(pTFCE=pTFCE, CLUST=CLUST, PVC=PVC, thr=threshs, V=V, Rd=Rd, Z=img))
}

# helper functions

#' Aggregate logpvals by the "equidistant incremental logarithmic probability aggregation".
#'
#' @param logpvals vector of enhnaced -logP values correspoinding tho various thresholds
#' @param d delta -logP (giving the equidistant distribution in the log-space)
#'
#' @return aggregated -logP probability
#' @export
#'
#' @examples
aggregate.logpvals=function(logpvals, d) #
{
  logpvals[is.infinite(logpvals)]=745
  s=sum(logpvals)
  return( 0.5*(sqrt(d*(8*s+d))-d) )
}

#' expected value of cluster size at threshold h in img of size V with RESEL count Rd
#' mainly ported from the source code of FSL (https://fsl.fmrib.ox.ac.uk/fsl/fslwiki)
#'
#' @param h image height, that is, Z score threshold
#' @param V Number of voxels
#' @param Rd Rd
#'
#' @return Expected value of cluster size
#'
#' @examples
Es=function(h, V, Rd)
{
  h2=h*h
  ret=  ( log(V)+pnorm(h, log.p = T, lower.tail = F) )
  h2=h2[h>=1.1]
  ret[h>=1.1] = ret[h>=1.1] - ( log(Rd)+log(h2-1)-h2/2+-2*log(2*pi) )
  return( exp(ret) )

}

dvox=function(h) # PDF of Z threshold/voxel value
{
  dnorm(h )
}

pvox=function(h) # p-value for voxel value
{
  pnorm(h, lower.tail = F)
}

dclust=function(h, V, Rd, c) # PDF of cluster extent, given h thershold
{
  if (is.na(c))
    return(dvox(h))
  dcl=function(h, V, Rd, c)
  {
    lambda=(Es(h, V, Rd)/gamma(2.5) )^(-2/3)
    dclust=lambda*exp(-lambda*c^(2/3))
    #dclust[dclust==0]=.Machine$double.xmin # underflow happened
    dclust[h<1.3]=0
    dclust
  }
  dcl(h, V, Rd, c)/integrate(function(x){dcl(x, V, Rd, c)}, -Inf, Inf)$value
}

dvox.clust=function(h, V, Rd, c) # PDF of Z threshold value given cluster size
{
  if (is.na(c))
    #return(rep(0, length(h)))
    return(dvox(h))
  if (is.nan(dclust(h, V, Rd, c)[1])) # underflow
    return(dvox(h)*dclust(h, V, Rd, c))

  dvox(h)*dclust(h, V, Rd, c)/integrate(function(x){dvox(x)*dclust(x, V, Rd, c)}, -Inf, Inf)$value
}

pvox.clust=function(V, Rd, c, actH) # p-value for Z threshold value given cluster size
{
  if (actH<=1.3)
    return(pvox(actH)) # GRF theory might not apply at low thresholds

  if (is.nan(dvox.clust(actH, V, Rd, c)[1])) # underflow
    return(exp(-745)) #TODO: make machnie independent

  integrate(function(x){dvox.clust(x, V, Rd, c)}, actH, Inf)$value
}


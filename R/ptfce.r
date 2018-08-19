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
#' @references
#' T. Spisák, Z. Spisák, M. Zunhammer, U. Bingel, S. Smith, T. Nichols, T. Kincses, Probabilistic TFCE: a generalized combination of cluster size and voxel intensity to increase statistical power, under review.
#'
#' https://github.com/spisakt/pTFCE
#'
#' @param img Nifti Z-score image to enhance ("nifti" class from "oro.nifti" package)
#' @param mask Mask
#' @param Nh Number of thresholds
#' @param Rd Resel count (as output by FSL smoothest)
#' @param V Number of voxels in mask (as output by FSL smoothest)
#' @param residual 4D residual data for a better estimate of image smoothness (optional)
#' @param logpmin min threshold
#' @param logpmax max threshold
#' @param ZestThr Cluster-forming Z threshold below which P(h|c) is estimated as P(h), due to limitation of GRF theory. (default: 1.3)
#' @param verbose boolean: print progress bar and diagnostic messages if true (default)
#'
#' @details The function takes a Z-score image and a mask image (both "nifti" object of the oro.nifti package) as obligatory inputs.
#' Mask can be either binary or continous, in the latter case it will be thresholded at 0.5.
#' Smoothness information Rd (voxels per RESEL) and number of voxels V are optional and are to be interpreted as in FSL "smoothest".
#' If not specified, these values are estimated from the data, internally via the smoothest() function of this package, which is a direct port of the corresponding FSL function.
#' If Rd and/or V is not specified, and residual is specified, image smoothness will be determined basedd on teh 4D residual data (more accurate, see ?smoothest()).
#' The default value of the parameter Nh should work with images having usual ranges of Z-scores. At low values, although the processing becomes faster, the estimation of the enhanced values might become inaccurate and the enhanced p-value distribution becomes increasingly non-uniform. It is not recommended to set it lower than 30.
#' The parameters logpmin and logpmax define the range of values the incremental thresholding procedure covers. By default, these are based on the input data.
#'
#' @return An object of class "ptfce" is a list containing at least the following components:
#' \itemize{
#' \item{\strong{p} (uncorrected) pTFCE enhanced p-values}
#' \item{\strong{logp}  negative logarithm of pTFCE enhanced p-values}
#' \item{\strong{Z} pTFCE enhanced p-values converted to Z-scores}
#' \item{\strong{number_of_resels} number of resels (input for e.g. ptoz) }
#' \item{\strong{fwer0.05.Z} Z-score threshold corresponding to the corrected p<0.05 threshold, controlled for FWER}
#' }
#'
#' @export
#'
#' @examples Z=readNIfTI("Zmap.nii.gz");
#'
#' MASK=readNIfTI("mask.nii.gz");
#'
#' pTFCE=ptfce(Z, MASK);
#'
#' smooth=read.table("smoothness.txt"); # e.g. output by FSL smoothest
#' pTFCE=ptfce(Z, V=smooth[2,2], Rd = smooth[1,2]*smooth[2,2], mask = MASK, Nh = 100);
#' # plot if you want
#' orthographic(pTFCE$pTFCE, datatype=16) )
ptfce=function(img, mask, Nh=100, Rd=NA, V=NA, resels=NA, residual=NA,  logpmin=0, logpmax=-log(pnorm(max(img), lower.tail = F)), ZestThr=1.3, verbose=T )
{
  if (length(img[is.na(img)])>0)
  {
    warning("NAs detected and replaced with zero!")
    img[is.na(img)] = 0
  }
  if (length(img[!is.finite(img)])>0)
  {
    warning("Infinite values detected and replaced with zero!")
    img[!is.finite(img)] = 0
  }
  autosmooth=F
  if (is.na(Rd) || is.na(V))
  {
    autosmooth=T
    if (!is.na(residual))
      warning("Smoothness estimatmion based on 4D residual data is not yet implemented!")
    smooth=smoothest(img, mask, verbose = verbose)
    V=smooth$volume
    Rd=smooth$dLh*V
  }
  logp.thres=seq(logpmin, logpmax, length.out=Nh)
  dh=logp.thres[2]-logp.thres[1]
  p.thres=exp(-logp.thres)
  threshs=qnorm(p.thres, lower.tail = F)
  threshs[1]=-9999999999 # exp(-Inf):=0 #TODO: machine minimum
  ndh = length(threshs)
  CLUST=array(NA, dim=c(dim(img), length(threshs)))
  PVC=array(1, dim=c(dim(img), length(threshs)))

  if (verbose) cat("* Performing pTFCE...\n")
  if (verbose) pb <- txtProgressBar(style = 3)
  for (hi in 1:length(threshs))
  {
    if (verbose) setTxtProgressBar(pb, hi/length(threshs))
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
      PVC[,,,hi][sizes==size]=pvox.clust(V, Rd, size, h, ZestThr=ZestThr)
    }
    # plot if you want
    # orthographic(nifti(CLUST[,,, hi], datatype=16) )
  }
  # calculate pTFCE
  pTFCE=array(apply(PVC, c(1,2,3), function(x){exp( -aggregate.logpvals(-log(x), dh) )}), dim=dim(img))

  # copy nifit header information
  snames = methods::slotNames(img)
  snames = snames[ !snames %in% c(".Data", "dim_") ]
  pTFCE = nifti(img = pTFCE, dim = dim(pTFCE))
  class(pTFCE) = class(img)
  for (islot in snames) {
    methods::slot(pTFCE, islot) = methods::slot(img, islot)
  }
  # done copying header

  pTFCE[pTFCE==0]=.Machine$double.xmin #underflow
  pTFCE[pTFCE==1]=1-.Machine$double.neg.eps #underflow
  if (verbose) close(pb)
  if (autosmooth)
  {
    number_of_resels=smooth$volume / smooth$resels
    fwer0.05.Z=fwe.p2z(number_of_resels, 0.05)
  }
  else
  {
    if (is.na(resels))
    {
      warning("For GRF-based FWER correction, please specify resels, or use smoothness estimation based on the data, by not specifying Rd and V!")
      number_of_resels=NA
      fwer0.05.Z=NA
    }
    else
    {
      number_of_resels=V / resels
      fwer0.05.Z=fwe.p2z(number_of_resels, 0.05)
    }

  }
  return(list(p=pTFCE,
              logp=-log(pTFCE),
              Z=qnorm(pTFCE, lower.tail = F),
              number_of_resels=number_of_resels,
              fwer0.05.Z=fwer0.05.Z
              )
         )
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
#' @param Rd Rd (dLh*V)
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

pvox.clust=function(V, Rd, c, actH, ZestThr=1.3) # p-value for Z threshold value given cluster size
{
  if (actH<=ZestThr)  # ZestThr: GRF theory might not apply at low thresholds
    return(pvox(actH))

  if (is.nan(dvox.clust(actH, V, Rd, c)[1])) # underflow
    return(exp(-745)) #TODO: make machnie independent

  integrate(function(x){dvox.clust(x, V, Rd, c)}, actH, Inf)$value
}

#' Estimate global image smoothness.
#'
#' Ported form the C++ implementation of FSL.
#'
#' For mathematical background, see:
#' https://www.fmrib.ox.ac.uk/datasets/techrep/tr00df1/tr00df1/index.html
#'
#' @param img Z-score image or a 4D residual image, ("nifti" class from
#'   "oro.nifti" package)
#' @param mask image mask ("nifti" class from "oro.nifti" package)
#' @param dof degrees of freedom, obligatory if img is a 4D residual image
#' @param verbose boolean: print progress bar and diagnostic messages if true
#'   (default)
#'
#' @return An object of class "smoothness" is a list containing at least the
#'   following components:
#'   \itemize{
#'   \item \strong{volume}   volume of the mask used for estimating smoothness, in voxels
#'
#'   \item \strong{sigmasq}  sigma squared values in the x, y and z directions
#'
#'   \item \strong{FWHM}     full width at half maximum values of smoothness in the x, y and z
#'   direction (voxels)
#'
#'   \item \strong{dLh}      determininant of Lambda to the half (voxels^-3)
#'
#'   \item \strong{resels}   resel size (voxels per resel) }
#' @export
#'
#' @details The function takes two images (both "nifti" object of the oro.nifti
#'   package): (i) either a Z-score image or a 4D residual image together with
#'   the degrees of freedom, and (ii) a mask image as obligatory inputs. Mask
#'   can be either binary or continous, in the latter case it will be
#'   thresholded at 0.5.
#'
#'   For a Gaussian random field the smoothness is defined as  \deqn{W
#'   =|\Lambda|^{-1/2D} } where D is the dimensionality of the field and
#'   \eqn{\Lambda} the covariance matrix of it's first partial derivatives.
#'
#'   Using Z-score image is less optimal because:
#' \itemize{
#' \item Smoothness estimates need spatial derivatives, which are very noisy
#'   quantities and, for a single 3D map, can be computed just once on each
#'   direction.
#' \item The z-map may contain effects, and these affect smoothness.
#' }
#'
#'
smoothest=function(img, mask, dof=NA, verbose=T)
{
  if (verbose) cat("* Estimating smoothness based on the data...\n")
  if (verbose) pb <- txtProgressBar(style = 3)

  stand = standardise(mask, img);
  mask_volume=stand$count
  #mask=stand$mask
  #img=stand$R
  N=0

  usez=T
  if (dim(img)[3] <= 1)
  {
    usez = F;
    warning("using 2d image mode, but it is not tested yet")
  }

  # Estimate the smoothness of the normalised residual field
  # see TR00DF1 for mathematical description of the algorithm.
  X = 1
  Y = 2
  Z = 3
  SSminus = c(0, 0, 0)
  S2 = c(0, 0, 0)

  # adjust zstart for 2d mode
  zstart=2
  if (!usez) zstart=1
  for (z in 2:dim(img)[3])
  {
    if (verbose) setTxtProgressBar(pb, z/dim(img)[3])
    for (y in 2:dim(img)[2])
    {
      for (x in 2:dim(img)[1])
      {
        # Sum over N
        if( (mask[x, y, z]>0.5) &&
            (mask[x-1, y, z]>0.5) &&
            (mask[x, y-1, z]>0.5) &&
            ( (!usez) || (mask[x, y, z-1]>0.5) ) )
        {
          N=N+1

          if (!is.na(dim(img)[4])) # using res4D
          {
            for ( t in 1:dim(img)[4])
            {
              # Sum over M
              SSminus[X] =  SSminus[X] + img[x, y, z, t] * img[x-1, y, z, t]
              SSminus[Y] = SSminus[Y] + img[x, y, z, t] * img[x, y-1, z, t]
              if (usez) SSminus[Z] = SSminus[Z] + img[x, y, z, t] * img[x, y, z-1, t]

              S2[X] = S2[X] + 0.5 * (img[x, y, z, t]*img[x, y, z, t] + img[x-1, y, z, t]*img[x-1, y, z, t])
              S2[Y] = S2[Y] + 0.5 * (img[x, y, z, t]*img[x, y, z, t] + img[x, y-1, z, t]*img[x, y-1, z, t])
              if (usez) S2[Z] = S2[Z] + 0.5 * (img[x, y, z, t]*img[x, y, z, t] + img[x, y, z-1, t]*img[x, y, z-1, t])
            }
          }
          else # using Z-score map
          {
            SSminus[X] =  SSminus[X] + img[x, y, z] * img[x-1, y, z]
            SSminus[Y] = SSminus[Y] + img[x, y, z] * img[x, y-1, z]
            if (usez) SSminus[Z] = SSminus[Z] + img[x, y, z] * img[x, y, z-1]

            S2[X] = S2[X] + 0.5 * (img[x, y, z]*img[x, y, z] + img[x-1, y, z]*img[x-1, y, z])
            S2[Y] = S2[Y] + 0.5 * (img[x, y, z]*img[x, y, z] + img[x, y-1, z]*img[x, y-1, z])
            if (usez) S2[Z] = S2[Z] + 0.5 * (img[x, y, z]*img[x, y, z] + img[x, y, z-1]*img[x, y, z-1])

          }
        } # endif mask
      } #endz
    } #endy
  } #endx
  norm = 1.0/N
  v = dof # v - degrees of freedom (nu)
  if (!is.na(dim(img)[4])) # using res4D
  {
    print(paste("Non-edge voxels = ", N))
    print(paste("(v - 2)/(v - 1) = ", (v - 2)/(v - 1)))

    norm = (v - 2) / ((v - 1) * N * dim(img)[4]);
  }
  #print(paste("SSminus[X]", SSminus[X], "SSminus[Y]", SSminus[Y],"SSminus[Z]", SSminus[Z],"S2[X]", S2[X],"S2[Y]", S2[Y],"S2[Z]", S2[Z]))

  # for extreme smoothness
  if (SSminus[X]>=0.99999999*S2[X])
  {
    SSminus[X]=0.99999*S2[X]
    warning("Extreme smoothness detected in X - possibly biased global estimate.")
  }
  if (SSminus[Y]>=0.99999999*S2[Y])
  {
    SSminus[Y]=0.99999*S2[Y]
    warning("Extreme smoothness detected in Y - possibly biased global estimate.")
  }
  if (usez)
  {
    if (SSminus[Z]>=0.99999999*S2[Z])
    {
      SSminus[Z]=0.99999*S2[Z]
      warning("Extreme smoothness detected in Z - possibly biased global estimate.")
    }
  }

  # Convert to sigma squared
  sigmasq=rep(NA,3);

  sigmasq[X] = -1.0 / (4 * log(abs(SSminus[X]/S2[X])))
  sigmasq[Y] = -1.0 / (4 * log(abs(SSminus[Y]/S2[Y])))
  if (usez)
  {
    sigmasq[Z] = -1.0 / (4 * log(abs(SSminus[Z]/S2[Z])))
  }
  else
  {
    sigmasq[Z]=0
  }
  names(sigmasq)=c("x", "y", "z")

  # the following is determininant of Lambda to the half
  # i.e. dLh = | Lambda |^(1/2)
  # Furthermore, W_i = 1/(2.lambda_i) = sigma_i^2 =>
  #   det(Lambda) = det( lambda_i ) = det ( (2 W_i)^-1 ) = (2^D det(W))^-1
  #   where D = number of dimensions (2 or 3)
  if (usez)
  {
    dLh=((sigmasq[X]*sigmasq[Y]*sigmasq[Z])^-0.5)*(8^-0.5)
  }
  else
  {
    dLh=((sigmasq[X]*sigmasq[Y])^-0.5)*(4^-0.5)
  }
  #correcting for temporal DOF!!!!!
  if(!is.na(dim(img)[4]))
    dLh = dLh * interpolate(v);

  # Convert to full width half maximum
  FWHM=rep(NA,3)
  FWHM[X] = sqrt(8 * log(2) * sigmasq[X])
  FWHM[Y] = sqrt(8 * log(2) * sigmasq[Y])
  if (usez)
  {
    FWHM[Z] = sqrt(8 * log(2) * sigmasq[Z])
  }
  else
  {
    FWHM[Z]=0
  }
  names(FWHM)=c("x", "y", "z")

  resels = FWHM[X] * FWHM[Y]
  if (usez)
    resels = resels * FWHM[Z]

  if (verbose) close(pb)

  return(list(volume=mask_volume,
              sigmasq=sigmasq,
              FWHM=FWHM,
              dLh=dLh,
              resels=resels))
}

standardise=function(mask, R)
{
  count=0
  M=dim(R)[4]
  count=sum(mask[mask>0.5])
  if( !is.na(M) )
  {
    # For each voxel
    #    standardise data in 4th direction
    R=aperm(apply(R, c(1,2,3), scale), c(2,3,4,1))
  }#endif 4d
  #TODO: updtae mask if it differs from the actual data
  return(list(count=count, maks=mask, R=R))
}

interpolate=function(v)
{
  lut=rep(NA,500)
    lut[5]   = 1.5423138; lut[6]   = 1.3757105; lut[7]   = 1.2842680;
    lut[8]   = 1.2272151; lut[9]   = 1.1885232; lut[10]  = 1.1606988;
    lut[11]  = 1.1398000; lut[12]  = 1.1235677; lut[13]  = 1.1106196;
    lut[14]  = 1.1000651; lut[15]  = 1.0913060; lut[16]  = 1.0839261;
    lut[17]  = 1.0776276; lut[18]  = 1.0721920; lut[19]  = 1.0674553;
    lut[20]  = 1.0632924; lut[25]  = 1.0483053; lut[30]  = 1.0390117;
    lut[40]  = 1.0281339; lut[50]  = 1.0219834; lut[60]  = 1.0180339;
    lut[70]  = 1.0152850; lut[80]  = 1.0132621; lut[90]  = 1.0117115;
    lut[100] = 1.0104851; lut[150] = 1.0068808; lut[200] = 1.0051200;
    lut[300] = 1.0033865; lut[500] = 1.0020191;

  if (v<6) return(1.1) # ?? no idea - steve ??

  if (v>500)
  {
    retval=1.0321/v + 1
  }
  else
  {
    lut.low=lut[1:floor(v)]
    lut.high=lut[ceiling(v):length(lut)]

    j.first=which(lut==min(lut.low, na.rm = T))
    j.second=lut[j.first]

    i.first=which(lut==max(lut.high, na.rm = T))
    i.second=lut[i.first]

    retval = (j.second - i.second)/(j.first - i.first)*(v - i.first) + j.second
  }
  return(retval^0.5)

}

#' Convert Z-score value to FWER corrected p-value
#'
#' https://github.com/spisakt/pTFCE
#'
#' @param resel_count resel count
#' @param Z z-score value
#'
#' @details The Z-score value is converted to to p-value, corrected for multiple comparisons by controlling for famili-wise error rate
#' The parameter resel_count is the number of resels (resolution elements) in the image, and can be obtained e.g by smoothest() (see examples).
#'
#' @return FWER-corrected p-value
#'
#' @export
#'
#' @examples
#' s=smoothest(zmap, mask)
#' fwe.z2p(resel_count=s$volume/s$resels, Z=2.3)
#'
#' pTFCE=ptfe(zmap, mask)
#' fwe.z2p(pTFCE$number_of_resels, Z=3.1)
fwe.z2p=function(resel_count, Z)
{
  # based on FSL
  #the probability of a family wise error is approximately equiv- alent to the expected Euler Characteristic
  if (Z<2)
    p=1 # Below z of 2 E(EC) becomes non-monotonic
  else
    p = resel_count * 0.11694 * exp(-0.5*Z*Z)*(Z*Z-1)  #/* 0.11694 = (4ln2)^1.5 / (2pi)^2 */
  return(min(p,1.0))
}

#' Get Z-score threshold for an FWER corrected p-value
#'
#' https://github.com/spisakt/pTFCE
#'
#' @param resel_count resel count
#' @param FWEP FWER corrected p-value
#'
#' @details The p-value, corrected for multiple comparisons by controlling for famili-wise error rate, is converted to a Z-score threshold.
#' The parameter resel_count is the number of resels (resolution elements) in the image, and can be obtained e.g by smoothest() (see examples).
#'
#' @return Z-score threshold
#'
#' @export
#'
#' @examples
#' s=smoothest(zmap, mask)
#' fwe.p2z(resel_count=s$volume/s$resels, FWEP=0.05)
#'
#' pTFCE=ptfe(zmap, mask)
#' fwe.p2z(pTFCE$number_of_resels, FWEP=0.05)
fwe.p2z=function(resel_count, FWEP=0.05)
{
  # based on FSL ptoz implementation
  p=FWEP/(resel_count*0.11694) #/* (4ln2)^1.5 / (2pi)^2 */

  #if (p<MINP) p=MINP;

  if (p>0.5)
    return(0)
  else
  {
      l=2
      u=100
      z=0
      while (u-l>0.001)
      {
        z=(u+l)/2
        pp=exp(-0.5*z*z)*(z*z-1)
        if (pp<p)
          u=z
        else
          l=z
      }

  }
  return(z)
}


Es=function(h, V, Rd) # expected value of cluster size at threshold h in img of size V with RESEL count Rd
{
  #h[h<=1.1]=1 # TODO: how to substitute h<=1 ????
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

ptfce=function(img, Rd, V, mask, length.out=50,   logpmin=0, logpmax=-log(pnorm(max(img), lower.tail = F)) )
{
  #print(paste("start: ", logpmin, logpmax, dh))
  logp.thres=seq(logpmin, logpmax, length.out=length.out)
  dh=logp.thres[2]-logp.thres[1]
  p.thres=exp(-logp.thres)
  threshs=qnorm(p.thres, lower.tail = F)
  threshs[1]=-9999999999 # exp(-Inf):=0
  ndh = length(threshs)
  
  CLUST=array(NA, dim=c(dim(img), length(threshs)))
  PVC=array(1, dim=c(dim(img), length(threshs)))
  
  #V=sum(img>0)

  for (hi in 1:length(threshs))
  {
    h=threshs[hi]
    #print(h)
    thr=img
    thr[img>h]=1
    thr[img<=h]=0
    ccc=components(thr, kernel = shapeKernel(3, dim = 3, type = "diamond"))
    t=table(ccc)
    sizes=t[ccc]
    sizes=sizes*as.numeric(mask)
    
    CLUST[,,,hi]=sizes
    
    for (size in sort(unique(sizes)) )
    {
      if (size==0)
        next;
      #print(paste(" * ", size))
      PVC[,,,hi][sizes==size]=pvox.clust(V, Rd, size, h)
    }
    # plot if you want
    #  orthographic(nifti(CLUST[,,, hi], datatype=16) )     
  }
  # calculate pTFCE
  pTFCE=array(apply(PVC, c(1,2,3), function(x){exp( -integrate.logpvals(-log(x), dh) )}), dim=dim(img))
  pTFCE[pTFCE==0]=4.940656e-324
  return(list(pTFCE=pTFCE, CLUST=CLUST, PVC=PVC, thr=threshs, V=V, Rd=Rd, Z=img))
}

analyzeAtThreshold=function(TFCE, x,y,z,thr)
{
  par(mfrow=c(3,1));
  h=TFCE$thr
  hx=h[-1]
  actH=hx[hx>thr][1]
  i=which(hx==actH)
  print(paste("Actual threshold: ", actH))
  #############################
  # TEST CASES
  #############################
  #x=55; y=22; z=17 #Z=3.7 FP
  #x=39; y=33; z=17 #Z=3.7 TP
  #############################
  #x=34; y=22; z=17 # Z=2.6, FP
  #x=29; y=35; z=17 # Z=2.6 TP
  #############################
  #x=27;y=31; z=20 # 2.38, FP leaker
  #x=56;y=24;z=19 #Z=2.39, FP
  #############################
  cx=TFCE$CLUST[x,y,z,];
  c=cx[i]
  
  plot(hx, dvox(hx), type="l",
       main="Prior Probability Density Function",
       xlab="Z threshold",
       ylab="PDF(threshold)");
  text(actH+0.3, 0.8*max(dvox(hx)), paste("p=", sprintf("%0.5f",integrate(dvox, actH, Inf)$value)))
  abline(v=actH, lty=2)
  
  
  
  plot(hx, dclust(hx, V, Rd, c), type="l",
       main=paste("Likelihood Function for cluster size:", c),
       xlab="Z threshold",
       ylab="Likelihood(cluster|threshold)");
  text(actH+0.3, 0.8*max(dclust(hx, V, Rd, c)), paste("p=", sprintf("%0.5f",integrate(function(x){dclust(x, V, Rd, c)}, actH, Inf)$value))) 
  abline(v=actH, lty=2)
  
  plot(hx, dvox.clust(hx, V, Rd, c), type="l",
       main="Posterior Probability Density Function",
       xlab="Z threshold",
       ylab="PDF(threshold|cluster)");
  text(actH+0.3, 0.8*max(dvox.clust(hx, V, Rd, c)), paste("p=", sprintf("%0.5f",integrate(function(x){dvox.clust(x, V, Rd, c)}, actH, Inf)$value))) 
  abline(v=actH, lty=2)
  print(integrate(function(x){dvox.clust(x, V, Rd, c)}, actH, Inf)$value)
  
  par(mfrow=c(1,1));
}

integrate.logpvals=function(logpvals, d) #given that d is equidistant in logpvals
{
  logpvals[is.infinite(logpvals)]=745
  s=sum(logpvals)
  return( 0.5*(sqrt(d*(8*s+d))-d) )
}


analyzeAllThresholds=function(TFCE, x, y, z)
{
  hx=TFCE$thr
  #print(hx)
  #############################
  # TEST CASES
  #############################
  #x=55; y=22; z=17 #Z=3.7 FP
  #x=39; y=33; z=17 #Z=3.7 TP
  #############################
  #x=34; y=22; z=17 # Z=2.6, FP
  #x=29; y=35; z=17 # Z=2.6 TP
  #############################
  #x=27;y=31; z=20 # 2.38, FP leaker
  #x=56;y=24;z=19 #Z=2.39, FP
  #############################
  cx=TFCE$CLUST[x,y,z,];
  cx[is.na(cx)]=0
  
  logpvc=-log(TFCE$PVC[x,y,z,])
  logpvc[cx==0]=0
  
  logpv=-log(pvox(hx[1:length(hx)]))
  logpv[cx==0]=0
  d=-log(pnorm(hx, lower.tail = F))[2]--log(pnorm(hx, lower.tail = F))[1]
  
  plot(-log(pnorm(hx, lower.tail = F)), logpv, type="p", col="red", 
       xlab="-logP(h)",
       ylab="-logP(voxel|c)",
       ylim=c(0, max(logpvc, logpv)))
  lines(-log(pnorm(hx, lower.tail = F)), logpvc, type="l" )
  text(max(-log(pnorm(hx, lower.tail = F)))*0.8,max(logpv )*0.2, paste("P(thresholds)=", sprintf("%0.5f", exp( -integrate.logpvals( logpv, d ) ) ) ) )
  text(max(-log(pnorm(hx, lower.tail = F)))*0.8,max(logpv )*0.5, paste("P(thresholds|clusters)=", sprintf("%0.5f", exp( -integrate.logpvals( logpvc, d ) ) ) ) )
  

  print(TFCE$Z[x,y,z])
  print((pnorm(TFCE$Z[x,y,z], lower.tail=F)))
  #print(paste( "P(thresholds) =",exp(-2*mean(logpv))) )
  print(paste( "P(thresholds) =",exp( -integrate.logpvals( logpv, d ) ) ) )
  
  print(paste( "P(thresholds) =",exp( -integrate.logpvals( logpvc, d ) ) ) )
}


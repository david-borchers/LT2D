#' @title Determine minimum forward distance for simulated data that will yield a given \code{(g0)<1}
#' 
#' @description This is a helper function for simulating data with \code{g(0)<1}.  The function
#'  calculates a minimum forward distance, \code{y_{min}}, to truncate simulated data, i.e. simulated sightings 
#'  with \code{y_i<y_{min}} are deleted.  
#'  
#'  @details The idea behind this function is that simulated data with \code{y_i<y_{min}} would got have been
#'  detected at a given \code{g(0)<1}.  The \code{y_{min}} for a given \code{g(0)<1} is solved using 
#'  \code{optimise()}.  The resulting \code{y_{min}} is substracted from the remaining simulated observations.
#'  @param gy desired detectability at \code{y=0}
#'  @param ymax largest forward distance
#'  @param b vector of detection hazard function parameters
#'  @param hfun detection hazard function
#'  @param tolerance passed into \code{optimize()}
#' @export  
#' @examples 
#' \dontrun{
#' ###No bump.
#'minY=yDistForg(gy=0.8,ymax=dfit.hn$ystart,b=dfit.hn$b,hfun=dfit.hn$hr)
#'minY
#'1-Sy(y=minY,x=0,ymax=dfit.hn$ystart,b=dfit.hn$b,hfun=dfit.hn$hr)
#'}

yDistForg=function(gy,ymax,b,hfun,tolerance=0.0001,...){
  #start by determining g(0)
  g0chk=1-Sy(y=0,x=0,ymax=ymax,b=b,hfun=hfun)
  if(gy>g0chk)
    stop('Model fit estimated g(0)=',round(g0chk,2),' select gy ARG less than estimated g(0)')
  gOPT=function(y,ymax,b,hfun,gy) {
    g=1-Sy(y=y,x=0,ymax=ymax,b=b,hfun=hfun)
    return(abs(g-gy))}
  opt=optimise(f = gOPT,interval=c(0,ystart),ymax=ymax,b=b,hfun=hfun,gy=gy,tol=tolerance)
  return(opt$minimum)
}

#'Data simulator 
#'
#'Simulate a fixed number of observations potentially when \code{g(0)<1}.  This is essentially a wrapper function for \code{\link{simXY}}
#'@param n=NULL if n is a positive integer, simulator will run until number of observations equals n
#'@param ymin=NULL if n is a positive real number, simulator will discard simulated positions <ymin
#'@param ... other arguments to be passed into \link{\code{simXY}}
#'@details  If n=NULL and ymin=NULL then use \link{\code{simXY}}.  When n!=NULL this function iteratively calls \link{\code{simXY}} until n simulated observations are reached.
#'When ymin!=NULL this function discards simulated data with Y-coordinates <ymin.  It is intended that \code{ymin}
#' is calculated by a call of \code{yDistForg} to determine the y distance at which g(y)= a required detection probability
#' @examples 
#' \dontrun{
#'minY=yDistForg(gy=0.8,ymax=dfit.hn$ystart,b=dfit.hn$b,hfun=dfit.hn$hr)
#'minY
#'simDatg0dot8=NULL
#'simDatg0dot8=simXY2(pi.x=dfit.hn$pi.x,logphi=dfit.hn$logphi,
#'                    hr=dfit.hn$hr,b=dfit.hn$b,w=dfit.hn$w,ystart=dfit.hn$ystart,
#'                    n=2000,ymin=minY)
#'min(simDatg0dot8$y)
#' }
simXY2=function(pi.x,logphi,hr,b,w,ystart,xSampL=1e3,discardNotSeen=TRUE,n=NULL,ymin=NULL,...)
{
  if(is.null(n) & is.null(ystart))
    stop('n and ymin are NULL call simXY() instead')
  hr=match.fun(hr);pi.x=match.fun(pi.x)
  Nrows=0
  while(Nrows<n){
    obs=simXY(N=1,pi.x=pi.x,logphi=logphi,hr=hr,b=b,w=w,ystart=ystart,discardNotSeen=discardNotSeen)
    obs=obs$locs
    if(nrow(obs)==0)
      next else{
        if(!is.null(ymin)) {
          obs=subset(obs,y>ymin)
          if(nrow(obs)==0) next
        }
      }
    if(Nrows==0)
      out=obs else out=rbind.data.frame(out,obs)
      Nrows=Nrows+1
  } #end while 
  if(!is.null(ymin)) out$y=out$y-ymin
return(out)
}

#' @title  Simulate observations and fit a model 
#'
#' @note Simulate observations from a hazard rate and perpendicular distribution and fit a model to the simulated observations.
#'Currently the function uses the same hazard rate, prependicular density distribuion and parameters for the simulated data and initial values for fitting. 
#'@param n=NULL if n is a positive integer, simulator will run until number of observations equals n
#'@param ymin=NULL if n is a positive real number, simulator will discard simulated positions <ymin
#'@param optimx use the \code{\link{optimx::optimx()}} function for fitting, Default is FALSE.
#'@param debug=FALSE additional output for debugging. See details
#'@param ... other arguments to be passed into \link{\code{simXY}}
#'@details  If n=NULL and ymin=NULL then use \link{\code{simXY}}.  When n!=NULL this function iteratively calls \link{\code{simXY}} until n simulated observations are reached.
#'When ymin!=NULL this function discards simulated data with Y-coordinates <ymin.  It is intended that \code{ymin}
#' is calculated by a call of \code{yDistForg} to determine the y distance at which g(y)= a required detection probability
#' If \code{debug=TRUE} a list object including the model fit is returned. 
#' @export
#' @author Martin J. Cox
#'@examples
#'system.time({out=simulator(n=60,fit=fit.n.ip1,optimx=FALSE)})
simulator=function(N=NULL,fit,n=NULL,ymin=NULL,optimx=FALSE,debug=FALSE)  {
  pi.x=match.fun(fit$pi.x);logphi=fit$logphi
  hr=match.fun(fit$hr);b=fit$b;w=fit$w;ystart=fit$ystart
  
  if(is.null(n) & is.null(ymin)){
    sim=simXY(N=N,pi.x=pi.x,logphi=logphi,
              hr=hr,b=b,w=w,
              ystart=ystart,xSampL=1e3,discardNotSeen=TRUE)
    x=sim$locs$x;y=sim$locs$y
  }else{
    sim=simXY2(pi.x=pi.x,logphi=logphi,
               hr=hr,b=b,w=w,
               ystart=ystart,xSampL=1e3,
               discardNotSeen=TRUE,n=n,ymin=ymin)
    x=sim$x;y=sim$y}
  
  #hr=match.fun(sim$settings$hr)
  #pi.x=match.fun(sim$settings$pi.x)
  #w=sim$settings$w
  if(optimx) {
    require(optimx)
    require(numDeriv)
    mod=try(fityxOptimx(y=y,x=x,b=b,hr=hr,ystart=ystart,
                        pi.x=pi.x,logphi=logphi,w=w,
                        itnmax=5000,hessian=TRUE))} else {
                          mod=try(fityx(y,x,
                                        b=b,
                                        hr=hr,
                                        ystart=ystart,
                                        pi.x=pi.x,
                                        logphi=logphi,
                                        w=w,
                                        hessian=TRUE,
                                        control=list(maxit=5000)))}
  
  if(length(mod)==1){
    warning('Model failed to converge')
    if(debug) return(list(simVal=rep('convFail',25),fit=mod))
    return(rep('convFail',25))
  }
  if(any(is.na(mod$CVpar))) {warning('Hessian issue')
    if(debug) return(list(simVal=rep('H',25),fit=mod))
    return(rep('H',25))
  }
  
  
  res=phatInterval(mod)
  
  
  res$w=w
  res$Varphat=(res$CV.phat*res$phat)**2
  res$effectiveStrip=res$phat*w
  res$VareffectiveStrip=res$Varphat*w**2
  res$CVeffectiveStrip=sqrt(res$VareffectiveStrip)/res$effectiveStrip
  
  res$N=N
  res$n=length(y)
  #1-Sy(x=0,y=0,ymax=ystart,b=tfit$b,hfun=hr)
  res$g0hat=1-Sy(0,0,ymax=ystart,b=mod$b,hfun=hr)
  #res$Nhat2DLT=n/res$phat
  #res$relBias2DLT=(res$Nhat2DLT-res$N)/res$N
  #res$NhatLow2DLT=n/res$upper.bound
  #res$NhatHigh2DLT=n/res$lower.bound
  #res$covered2DLT=TRUE
  #res$covered2DLT[res$NhatLow2DLT>res$N]=FALSE
  #res$covered2DLT[res$NhatHigh2DLT<res$N]=FALSE
  #fwdGoF=GoFy(fit=mod)
  #res$LT2DFWDp.cvm=fwdGoF$p.cvm
  #res$LT2DFWDp.kol=fwdGoF$p.kolomogarov
  #perpGoF=GoFy(fit=mod)
  #res$LT2DPERPp.cvm=perpGoF$p.cvm
  #res$LT2DPERPp.kol=perpGoF$p.kolomogarov
  # Fit using CDS:
  # Must have column names "Region.Label", "Area", "Sample.Label", "Effort", "distance", "size"
  n.all=length(x)#sim$locs$x)
  ones=rep(1,n.all)
  ddat=data.frame(Region.Label=ones, Area=ones, Sample.Label=ones, Effort=ones, 
                  distance=x, size=ones)
  fit.hr<-ds(ddat[ddat$distance<=w,],key="hr",adjustment=NULL)
  
  summary(fit.hr)
  #gofCDS=ddf.gof(fit.hr$ddf)$dsgof
  
  fit.hr.sum=summary(fit.hr$ddf)
  res$phatCDS=fit.hr.sum$average.p
  res$phatCVCDS=fit.hr.sum$average.p.se/fit.hr.sum$average.p
  #res$NhatCDS=fit.hr$ddf$Nhat # abundance in covered region
  #res$NCVCDS=as.vector(fit.hr.sum$Nhat.se/fit.hr$ddf$Nhat)
  #res$NhatLowCDS=as.vector(res$NhatCDS-1.96*fit.hr.sum$Nhat.se)
  #res$NhatHighCDS=as.vector(res$NhatCDS+1.96*fit.hr.sum$Nhat.se)
  #res$coveredCDS=TRUE
  #res$coveredCDS[res$NhatLowCDS<res$N]=FALSE
  #res$coveredCDS[res$NhatHighCDS>res$N]=FALSE
  #res$relBiasCDS=(res$NhatCDS-res$N)/res$N
  
  #res$gofkspCDS=gofCDS$ks$p
  #res$gofCvMpCDS=gofCDS$CvM$p
  
  par=mod$par
  names(par)=paste('parEst',1:length(par),sep='')
  parCV=mod$CVpar
  names(parCV)=paste('parEstCV',1:length(par),sep='')
  res=cbind.data.frame(res,t(par),t(parCV))
  if(debug)
    return(list(simVal=res,fit=mod))
  return(res)
}
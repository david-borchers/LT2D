#'@title Maximum likelihood estimation for unknown hazard and perpendicular distance distribution
#'
#'@description Uses \code{\link{optim}} to obtain a MLE for the hazard function and animal
#'perpendicular distance distribution.  Functional forms for the hazard and perpendicular distance 
#'distribution must be specified.
#'
#'@param y forward distance observations
#'@param x perpendicular distance observations
#'@param b two to four element vector of hazard rate parameters, some of which may be logged  
#'@param hr hazard rate function 
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param pi.x perpendicular distance density distribution
#'@param logphi parameters for pi.x (some maybe logged)
#'@param w perpendicular truncation distance.
#'@param method optimisation method to be used by \code{\link{optimx}}. Default is \code{"Nelder-Mead"}
#'@param lower, upper Bounds for parameters for use with methods such as \code{"L-BFGS-B"}. See \code{\link{optimx}}.
#'@param control see \code{\link{optimx}} control
#'@param itnmax maximum number of iterations for the \code{Nelder-Mead} method.  NB this is not passed in via the control argument
#'@param hessian return hessian.  See also \code{\link{optimx}}.
#'@param corrFlag=0.7 Absolute parameter correlation value above which a warning is issued.
#'@param ... arguments to be passed into \code{\link{optimx}}
#'@return 
#'\code{\link{optim}} fit object and \cr
#'\code{$hr} = hazard rate function used.\cr
#'\code{$pi.x} = perpendicular distance function used.\cr
#'\code{$ystart} = ystart max forward distance detection used.\cr
#'\code{$w} = perpendicular truncation distance used.\cr
#'\code{$b} = estimated hazard parameters\cr
#'\code{$dat} = data frame with data (\code{$x} and \code{$y})\cr
#'\code{$logphi} \cr
#'\code{AIC} AIC value\cr
#'And if \code{hessian=TRUE}:\cr
#'\code{vcov} variance covariance matrix.  Will warn if there is a problem inverting 
#'the hessian.\cr
#'\code{CVpar} Coefficient of variation for each paramter estimate. \cr
#'\code{error} Boolean, \code{TRUE} if convergence!=0 or problem inverting the hessian, 
#'or parameter correlation is exceeded.\cr
#'@details Must to ensure the hazard function has decayed to (very close to) zero by 
#'\code{ystart}.
#'@export
#'
#'@examples
#'\dontrun{
#'#Data preparation:
#'w=0.03;ystart=0.05
#'library(xlsx)
#'dat=read.xlsx("~/Dropbox/packages/LT2D/data/Jantho Primate Line.xlsx",1)
#'x=dat$PP.Distance
#'y=dat$Forward.Distance
#'nas=which(is.na(y))
#'x=x[-nas]
#'y=y[-nas]
#'gtw=which(x>w)
#'x=x[-gtw]
#'y=y[-gtw]
#'#Example model fits:
#'
#'b=c(-7.3287948, 0.9945317)
#'logphi=c(.01646734, -4.67131261)
#'fit.n.optx=NULL
#'fit.n.optx=fityxOptimx(y,x,b=b,
#'            hr=h1,ystart=ystart,
#'            control=list(trace=5),hessian=TRUE,
#'            pi.x=pi.norm,logphi=logphi,w=w)
#'
#'b=c(-7.3329141,0.9948721)
#'logphi=c(-0.05,-4.7)
#'fit.chn.optx=NULL
#'fit.chn.optx=fityxOptimx(y,x,b=b,hr=h1,ystart=ystart,
#'              pi.x=pi.chnorm,logphi=logphi,w=w,itnmax=5000,
#'              hessian=TRUE,control=list(trace=5))
#'
# Normal bump with ip1 hazard function:
#'b=c(5.2919208, -0.2205593, 8.4701307)
#'logphi=c(0.01784102, -4.42209067)
#'fit.n.ip1.optx=NULL
#'fit.n.ip1.optx=fityxOptimx(y,x,b=b,hr=ip1,ystart=ystart,
#'                pi.x=pi.norm,logphi=logphi,w=w,
#'                hessian=TRUE,control=list(trace=6))
#'}
#'@seealso \code{\link{negloglik.yx}} \code{\link{fityx}} \code{\link{optimx}}
fityxOptimx=function(y,x,b,hr,ystart,pi.x,logphi,w,method="Nelder-Mead",
                     lower=-Inf, upper=Inf,control=list(),
                     itnmax=NULL,hessian=FALSE,corrFlag=0.7,...)
{
  require(optimx)
  require(numDeriv)
  
  if(is.function(hr)){
    eval(parse(text=fNameFinder(hr)))
    hrname=fName  
  }else{
    hrname=hr
    hr=match.fun(hr)
  }
  
  if(is.function(pi.x)){
    eval(parse(text=fNameFinder(pi.x)))
    piname=fName  
  }else{
    piname=pi.x
    pi.x=match.fun(pi.x)
  }
  if(piname=="pi.const") pars=b else pars=c(b,logphi)
  
  length.b=length(b)
 
  fitx=optimx(par=pars, fn=negloglik.yx, y=y,x=x,hr=hr,ystart=ystart,pi.x=pi.x,w=w,
                   length.b=length.b, 
                   method=method, lower=lower, upper=upper,itnmax=itnmax,
                   control=control,hessian=hessian,...)
  fit=list(par=as.vector(unlist(fitx[grep('p',names(fitx))])))
  
  fit$error=FALSE
  if(fitx$convcode!=0){
    warning('Convergence issue (code = ', fitx$convcode,') . Check optimx::optimx() help.')
    fit$error=TRUE
  }
  fit$hr=hrname
  fit$pi.x=piname
  fit$ystart=ystart
  fit$w=w
  fit$b=fit$par[1:length.b]  # ***
  if(length.b!=length(pars)){
    fit$logphi=fit$par[(1+length.b):length(pars)]
  }else{
    fit$logphi=NA
  }    # ***
  fit$AIC=2*fitx$value+2*length(fit$par)
  fit$dat=data.frame(x=x,y=y) # attach data to fitted object

  if(hessian){
    
    negloglik.yxH=function(x1,...) negloglik.yx(x=x1,...)
    fit$hessian=numDeriv::hessian(func=negloglik.yxH, 
                                  x=fit$par,y=y,x1=x,hr=hr,
                                  ystart=ystart,pi.x=pi.x,w=w,
                      length.b=length.b)
    
    
    mNames=paste('b',1:length.b,sep='')
    if(!all(is.na(fit$logphi))) mNames=c(mNames,paste('logphi',1:length(fit$logphi),sep=''))
    fit$vcov=solve(fit$hessian)
    if(any(diag(fit$vcov)<=0)){
      warning('Failed to invert hessian.  Model covergance problem in fityxOptimx?')
      fit$error=TRUE
      fit$CVpar=rep(NA,length(fit$par))} else {
        fit$CVpar=sqrt(diag(solve(fit$hessian)))/abs(fit$par)}
    fit$corr=cov2cor(fit$vcov)
    row.names(fit$corr)=mNames
    colnames(fit$corr)=mNames
    corr=fit$corr
    corr[upper.tri(corr,diag=TRUE)]=NA
    corrIND=which(abs(corr)>corrFlag,arr.ind=T)
    if(nrow(corrIND)){
      warning('absolute correlation exceeds ',corrFlag,' in parameter estimates: ',
              paste(paste(mNames[corrIND[,1]],mNames[corrIND[,2]],sep=' to '),collapse='; '))
      fit$error=TRUE}
  }
  return(fit)
}

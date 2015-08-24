#' @title Kolmogarov-Smirnov goodness-of-fit p-value.
#'
#' @description
#' Kolmogarov-Smirnov goodness-of-fit p-value calculation.
#
#' @param x value of Kolmogarov-distributed random variable at which to evaluate.
#' @param inf approximation to infinity (a large number).
#' @param dp approximation convergence criterion.
#' @return 1 - p-value for Kolmogarov distribution
#' @details
#' Calculates p-value for Kolmogarov distribution at x, approximating infite sum
#' \code{sqrt(2*pi)/x*sum{i=1}^infty exp(-(2i-1)^2*pi^2/(8x^2)))}
#' by a finite sum to inf (default 1000) if sum to inf and inf+1 differ by less
#' than dp (default 1e-4), else sum until difference is less than dp.
p.kolomogarov=function(x,inf=1000,dp=1e-4)
{
  infsum=rep(0,inf)
  i=1:inf
  K=sqrt(2*pi)/x
  p=p1=K*sum(exp(-(2*i-1)^2*pi^2/(8*x^2)))
  dp=1
  while(dp>1e-4) {
    inf=inf+1
    p=p1+K*exp(-(2*inf-1)^2*pi^2/(8*x^2))
    dp=abs(p-p1)
  }
  return(p=p)
}

#'@title Goodness-of-Fit in the forward direction (y)
#'
#' @description Calculates the Goodness-of-Fit in the forward direction (y)
#'
#'@description Plot f(y) and forward distance distribution
#'resulting from a call of \code{\link{fityx}}.
#'
#'@param fit object resulting from a call of \link{fityx}.
#'@param plot=FALSE boolean TRUE - display Q-Q plot.
#'@details Calculates goodness-of-fit for forward distances and calculates p-values 
#'using kolomogarov and Cramer-von Mises.
#'@seealso \code{\link{fityx}}, \code{\link{}}
#'
#'@return
#'data frame with these elements
#'\code{$p.cvm} = Cramer-von Mises p-value.
#'\code{$D.kolomogarov} = x value of Kolmogarov-distributed random variable
#'\code{$p.kolomogarov} = kolomogarov p-value.
#'@examples
#'\dontrun{
#'ystart=0.05;w=0.03
#'logphi=c(0.0180552, -4.4215995)
#'b=c(-23.725809,  -3.136638  , 2.122910)
#'N=200 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.norm,logphi=logphi,
#'hr=ip1,b=b,w=w,ystart=ystart)
#'
#'x=simDat$locs$x; y=simDat$locs$y
#'est.yx=fityx(y,x,b,hr=ip1,ystart,pi.x=pi.norm,logphi,w)
#'plotfit.y(y,x,est.yx,nclass=10)
#'GoFy(fit=est.yx,plot=TRUE)
#'}
#'@seealso \code{\link{fityx}} \code{\link{p.kolomogarov}}
#'\code{\link{GoFx}}
#'@export
GoFy=function(fit,plot=FALSE,dotitle=FALSE){
  ystart=fit$ystart
  w=fit$w
  hr=match.fun(fit$hr)
  if(fit$hr=="h1") analytic.F0=TRUE else analytic.F0=FALSE
  ystart=fit$ystart
  pi.x=match.fun(fit$pi.x)
  b=fit$b
  logphi=fit$logphi
  
  x=fit$dat$x;y=fit$dat$y
  n=length(x)
  if(analytic.F0) {
    zeros=(y==0)
    Fy=rep(NA,length(x))
    if(length(zeros)>0){
      Fy[zeros]=HBhr(x[zeros],h1.to.HB(b))
      Fy[!zeros]=(1-Sy(x=x[!zeros],y=y[!zeros],ymax=ystart,b=b,hfun=hr))
    }
    F0=HBhr(x,h1.to.HB(b))
#    F0=(1-Sy(x=x,y=rep(0,length(y)),ymax=ystart,b=b,hfun=hr))
  } else {
    Fy=(1-Sy(x=x,y=y,ymax=ystart,b=b,hfun=hr))
    F0=(1-Sy(x=x,y=rep(0,length(y)),ymax=ystart,b=b,hfun=hr))
  }
    
  Fy0=Fy/F0
  Fy0.order=order(Fy0)
  yy=y[Fy0.order]
  xx=x[Fy0.order]
  cdf=Fy0[Fy0.order]
  e.cdf=order(cdf)/n
  
  # K-S statistic
  dF=cdf-e.cdf
  worst=which(abs(dF)==max(abs(dF))) # mark point on which Kolmogarov test hinges
  Dn=max(abs(dF))*sqrt(n)
  
  p.cvm=goftest::cvm.test(Fy0)$p.value
  p.kolomogarov=1-p.kolomogarov(Dn)
  
  if(plot){
    main=""
    if(dotitle) main="Forward Dist. Q-Q Plot"
    plot(e.cdf,cdf,xlab="Empirical Distribution Function",
         ylab="Cumulative Distribution Function",
         main=main,
         xlim=c(1,0),ylim=c(0,1),pch="+")
    lines(c(0,1),c(0,1))
    if(dotitle) mtext(paste('p-values: Cramer-von Mises=',round(p.cvm,2),' ;
                kolomogarov=',round(p.kolomogarov,2)))
    points(e.cdf[worst],cdf[worst],col="red") # mark point on which Kolmogarov test hinges
  }
  
  pvals=c(p.cvm,p.kolomogarov);names(pvals)=c("Cramer-von Mises","Kolmogarov-Smirnov")
  return(data.frame(pvals=pvals,D.kolomogarov=Dn))
}


#' @title Parameter conversion for Hayes+Buckland hazard rate model.
#'
#' @description
#' Converts parameters from the 2-dimensional a*r^{-b} form of the Hayes and Buckland (1983) 
#' hazard rate model used in \link{\code{h1}}, to parameters of the perpendicular distance form 
#' of the model used in conventional distance sampling.
#
#' @param b parameter vector of the \link{\code{h1}}.
#' @return
#' parameter vector c(a1,b1) for Hayes and Buckland (1983) hazard rate model 
#' g(x)=1-exp(-a1*x^(-(b1-1))).
#' @seealso \code{\link{h1}}
#' @examples
#' h1.to.HB(c(-7, 1))
h1.to.HB=function(b){
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  b1=exp(b[2])
  a=exp(b[1])
  a1=a*gamma((b1-1)/2)*gamma(0.5)/(2*gamma(b1/2))
  return(c(a1,b1))
}

#' @title Hayes+Buckland hazard rate model.
#'
#' @description
#' Evaluates the Hayes and Buckland (1983) hazard rate model (HB model) as parameterised on 
#' page 36 of that paper.
#
#' @param x perpendicular distance(s).
#' @param theta parameter vector with theta[1] being HB model parameter a1 and theta[2] being b.
#' @return
#' g(x)=1-exp(-a1*x^(-(b-1))).
#' @seealso \code{\link{h1.to.HB}}
#' @examples
#' xx=seq(0,0.03,length=100)
#' p=HBhr(xx,h1.to.HB(c(-7, 0.85)))
#' plot(xx,p,type="l",xlab="Perpendicular distance (x)",ylab="p(x)",ylim=c(0,1))
HBhr=function(x,theta) 1-exp(-theta[1]*x^(-(theta[2]-1)))



#' @title Goodness-of-fit in perpendicular dimension.
#'
#' @description
#' Calculates goodness-of-fit in perpendicular dimension, plots fit, and returns p-value and
#' other stuff. Returns two p-values: \code{p.ks} is the Kolmogarov-Smirnov p-value (which is
#' based on only the largest difference between emprical and theoretical cdfs), and 
#' Cramer-von Mises p-value (which is based on all cdf values).
#
#' @param hmltm fitted model, as output by \code{\link{est.hmltm}}
#' @param plot If TRUE, does Q-Q plot. Point corresponding to largest difference between
#' empirical and theoretical cdf (on which the Kolmogarov-Smirnov test is based) is circled in red.
#'@return
#'data frame with these elements
#'\code{$p.cvm} = Cramer-von Mises p-value.
#'\code{$D.kolomogarov} = x value of Kolmogarov-distributed random variable
#'\code{$p.kolomogarov} = kolomogarov p-value (which is
#' based on only the largest difference between emprical and theoretical cdfs).
#' \code{$qq.x} = empirical distribution function values.
#'\code{$qq.y} = cumulative distribution function values.
#'\code{$x} = x values.
#'@seealso \code{\link{fityx}} \code{\link{p.kolomogarov}} \code{\link{GoFy}}
#' @examples
#'\dontrun{
#'ystart=0.05;w=0.03
#'logphi=c(0.0180552, -4.4215995)
#'b=c(-23.725809,  -3.136638  , 2.122910)
#'N=200 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.norm,logphi=logphi,
#'hr=ip1,b=b,w=w,ystart=ystart)
#'
#'x=simDat$locs$x; y=simDat$locs$y
#'est.yx=fityx(y,x,b,hr=ip1,ystart,pi.x=pi.norm,logphi,w)
#' plotfit.x(x=x,est=est.yx)
#' rug(x=est.yx$dat$x)
#' tst=GoFx(fit=est.yx,plot=TRUE)
#' }
GoFx=function(fit,plot=FALSE,nint=100,dotitle=FALSE){
  ystart=fit$ystart
  w=fit$w
  hr=match.fun(fit$hr)
  ystart=fit$ystart
  pi.x=match.fun(fit$pi.x)
  b=fit$b
  logphi=fit$logphi
  
  x=fit$dat$x;y=fit$dat$y
  n=length(x)
  
  edf=(1:n)/n
  
  f0V=vector(mode = 'numeric',length=length(x))
  for(i in 1:length(x))
    f0V[i]=integrate(p.pi.x,0,x[i],b,hr,ystart,pi.x,logphi,w)$value
  
  Af0=integrate(p.pi.x,0,w,b,hr,ystart,pi.x,logphi,w)$value
  
  cdf=f0V/Af0
  
  cdf.order=order(cdf)
  cdf=cdf[cdf.order]
  e.cdf=cdf.order/n
  # K-S statistic
  dF=cdf-edf
  worst=which(abs(dF)==max(abs(dF)))
  Dn=max(abs(dF))*sqrt(n)
  p.ks=1-p.kolomogarov(Dn)
  p.cvm=goftest::cvm.test(cdf)$p.value # Under model, cdf values are from uniform; default for cvm.test is "punif"
  if(plot) {
    main=""
    if(dotitle) main="Perp. Dist. Q-Q Plot"
    plot(edf,cdf,pch="+",xlim=c(0,1),ylim=c(0,1),
         xlab="Empirical Distribution Function",
         ylab="Cumulative Distribution Function",main=main)
    lines(c(0,1),c(0,1))
    if(dotitle) mtext(paste('p-values: Cramer-von Mises=',round(p.cvm,2),' ;
                kolomogarov=',round(p.ks,2)))
    points(edf[worst],cdf[worst],col="red")
  }
  
  pvals=c(p.cvm,p.ks);names(pvals)=c("Cramer-von Mises","Kolmogarov-Smirnov")
  return(list(pvals=pvals,D.kolomogarov=Dn,qq.x=edf,qq.y=cdf,x=x[cdf.order]))
}
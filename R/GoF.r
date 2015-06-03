#' @title Kolmogarov-Smirnov goodness-of-fit p-value.
#'
#' @description
#' Kolmogarov-Smirnov goodness-of-fit p-value calculation.
#
#' @param x value of Kolmogarov-distributed random variable at which to
evaluate.
#' @param inf approximation to infinity (a large number).
#' @param dp approximation convergence criterion.
#' @return 1 - p-value for Kolmogarov distribution
#' @details
#' Calculates p-value for Kolmogarov distribution at x, approximating
infite sum
#' \code{sqrt(2*pi)/x*sum{i=1}^infty exp(-(2i-1)^2*pi^2/(8x^2)))}
#' by a finite sum to inf (default 1000) if sum to inf and inf+1 differ
by less
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
#'@depends \link{p.kolomogarov}
#'@details Calculates goodness-of-fit for forward distances and
calculates p-values using kolomogarov and Cramer-von Mises.
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
\code{\link{GoFx}}
#'@export
GoFy=function(fit,plot=FALSE){
  ystart=fit$ystart
  w=fit$w
  hr=match.fun(fit$hr)
  ystart=fit$ystart
  pi.x=match.fun(fit$pi.x)
  logphi=fit$logphi
  
  x=fit$dat$x;y=fit$dat$y
  n=length(x)
  Fy=(1-Sy(x=x,y=y,ymax=ystart,b=b,hfun=hr))
  F0=(1-Sy(x=x,y=rep(0,length(y)),ymax=ystart,b=b,hfun=hr))
  
  Fy0=Fy/F0
  Fy0.order=order(Fy0)
  yy=y[Fy0.order]
  xx=x[Fy0.order]
  cdf=Fy0[Fy0.order]
  e.cdf=order(cdf)/n
  
  # K-S statistic
  dF=cdf-e.cdf
  worst=which(abs(dF)==max(abs(dF))) # mark point on which Kolmogarov
  test hinges
  Dn=max(abs(dF))*sqrt(n)
  
  p.cvm=nortest::cvm.test(Fy0)$p.value
  p.kolomogarov=1-p.kolomogarov(Dn)
  
  if(plot){
    plot(1-e.cdf,cdf,xlab="Empirical Distribution Function",
         ylab="Cumulative Distribution Function",
         main="Forward Dist. Q-Q Plot",
         xlim=c(0,1),ylim=c(0,1),pch="+")
    lines(c(0,1),c(1,0))
    mtext(paste('p-values: Cramer-von Mises=',round(p.cvm,2),' ;
                kolomogarov=',round(p.kolomogarov,2)))
    points(1-e.cdf[worst],cdf[worst],col="red") # mark point on which
    Kolmogarov test hinges
  }
  
  return(data.frame(p.cvm=p.cvm,D.kolomogarov=Dn,p.kolomogarov=p.kolomogarov))
}

#' @title Goodness-of-fit in perpendicular dimension.
#'
#' @description
#' Calculates goodness-of-fit in perpendicular dimension, plots fit, and
returns p-value and
#' other stuff. Returns two p-values: \code{p.ks} is the
Kolmogarov-Smirnov p-value (which is
                            #' based on only the largest difference between emprical and theoretical
                            cdfs), and Cramer-von Mises
#' p-value (which is based on all cdf values).
#
#' @param hmltm fitted model, as output by \code{\link{est.hmltm}}
#' @param plot If TRUE, does Q-Q plot. Point corresponding to largest
difference between
#' empirical and theoretical cdf (on which the Kolmogarov-Smirnov test
is based) is circled in red.
#'@return
#'data frame with these elements
#'\code{$p.cvm} = Cramer-von Mises p-value.
#'\code{$D.kolomogarov} = x value of Kolmogarov-distributed random variable
#'\code{$p.kolomogarov} = kolomogarov p-value (which is
#' based on only the largest difference between emprical and theoretical
cdfs).
#' \code{$qq.x} = empirical distribution function values.
#'\code{$qq.y} = cumulative distribution function values.
#'\code{$x} = x values.
#'@depends nortest
#'@seealso \code{\link{fityx}} \code{\link{p.kolomogarov}}
\code{\link{GoFy}}
#' @examples
#' #'\dontrun{
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

GoFx=function(fit,plot=FALSE,nint=100){
  ystart=fit$ystart
  w=fit$w
  hr=match.fun(fit$hr)
  ystart=fit$ystart
  pi.x=match.fun(fit$pi.x)
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
  p.cvm=nortest::cvm.test(cdf)$p.value # Under model, cdf values are
  from uniform; default for cvm.test is "punif"
  if(plot) {
    plot(edf,cdf,pch="+",xlim=c(0,1),ylim=c(0,1),
         xlab="Empirical Distribution Function",
         ylab="Cumulative Distribution Function",main="Perp. Dist. Q-Q
         Plot")
    lines(c(0,1),c(0,1))
    mtext(paste('p-values: Cramer-von Mises=',round(p.cvm,2),' ;
                kolomogarov=',round(p.ks,2)))
    points(edf[worst],cdf[worst],col="red")
  }
  
  
  return(list(p.cvm=p.cvm,D.kolomogarov=Dn,p.kolomogarov=p.ks,qq.x=edf,qq.y=cdf,x=x[cdf.order]))
}
}
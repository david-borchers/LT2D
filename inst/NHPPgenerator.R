library(spatstat)

ymax=20
y=seq(0.001,ymax,length=200)
b=log(c(0.01,1))
h=h1(y,0,b)
plot(y,h,type="l")
plot(y,log(h),type="l")

ymin=ymax*1e-8
Sy=function(x,y,ymax,b,hfun) exp(-integrate(match.fun(hfun),y,ymax,x=x,b=b)$value)
fyx=function(y,x,ymin,ymax,hfun,b,pi.x,logphi,W,lscale=1){
  h=match.fun(hfun)
  pix=match.fun(pi.x)
  nx=length(x)
  ax=abs(x)
  if(length(y)!=nx) stop("Lengths of x and y must be the same.")
  px=pix(ax,logphi,W)
  f=p=p0=rep(NA,nx)
  for(i in 1:nx){
    p0[i]=1-Sy(ax[i],ymin,ymax,b,hfun)
    p[i]=1-Sy(ax[i],y[i],ymax,b,hfun)
    f[i]=h(y[i],ax[i],b)*(1-p[i])
  }
#  return(list(f=f,p=p,p0=p0))
  return(f*px*lscale)
}

calc.lpars=function(n,ymin,ymax,W,hfun,b,pi.x,logphi,nx=100,ny=100,inflate=1.05){
  x=seq(-W,W,length=nx)
  y=seq(ymin,ymax,length=ny)
  a=diff(x)[1]*diff(y)[1] # grid cell area
  lambda=outer(y,x,fyx,ymin=ymin,ymax=ymax,hfun=hfun,b=b,pi.x=pi.x,logphi=logphi,W=W)
  En=sum(lambda*a)
  lscale=n/En
  lmax=max(lambda*lscale)*inflate # bigger than observed in case higher between grid cells
  outlist=list(lscale=lscale,lmax=lmax)
  class(outlist)="ppscale"
  return(outlist)
}

sim.n=function(n,ymin,ymax,W,hfun,b,pi.x,logphi,fix.n=TRUE,intscale=NULL,nbuffer=NULL){
  # calculate scaling needed to give E[n]=n
  if(is.null(intscale)) 
    lpars=calc.lpars(n,ymin,ymax,W,hfun,b,pi.x,logphi)
  if(class(intscale)!="ppscale") stop("intscale must be class `ppscale' (output of calc.lpars).")
  window=owin(c(0,ymax),c(-W,W)) # create observation window
  lmax=lpars$lmax # maximum value of intensity function on grid used by calc.lpars
  lscale=lpars$lscale # multiplier required to get intensity with expected sampls size n
  if(is.null(nbuffer)){
    # increase E[n] by 25% to reduce prob that sample size is < n:
    nbuffer=ifelse(fix.n,1.25,1)
  }
  pp=rpoispp(fyx,lmax,window,ymin=ymin,ymax=ymax,b=b,hfun=hfun,pi.x=pi.x,logphi=logphi,W=W,lscale=lscale*nbuffer)
  if(fix.n) {
    while(pp$n<n){ # crude way of generating big enough sample size:
      pp2=rpoispp(fyx,lmax,window,ymin=ymin,ymax=ymax,b=b,hfun=h2,pi.x=pi.norm,logphi=logphi,W=W,lscale=lscale*nbuffer)
      pp$n=pp$n+pp2$n
      pp$x=c(pp$x,pp2$x)
      pp$y=c(pp$y,pp2$y)
    }
    pp$n=n
    pp$x=pp$x[1:n]
    pp$y=pp$y[1:n]
  }
  return(pp)
}

dat=sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi,intscale=NULL,nbuffer=NULL)
quartz(h=6,w=8)
plot(density(pp))
contour(density(pp),add=TRUE)
plot(pp,pch="+",cex=0.75,add=TRUE)
hist(abs(pp$y),xlab="Perpendicular distance",main="")
hist(pp$x,xlab="Forward distance",main="")

system.time(for(i in 1:20) dat<-sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi))
system.time(for(i in 1:20) dat<-sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi,intscale=lpars))

fy=function(y,x,ymin,ymax,hfun,b,intenscale=1){
  h=match.fun(hfun)
  pix=match.fun(pi.x)
  n=length(x)
  ax=abs(x)
  if(length(y)!=n) stop("Lengths of x and y must be the same.")
  f=p=p0=rep(NA,n)
  for(i in 1:n){
    p0[i]=1-Sy(ax[i],ymin,ymax,b,hfun)
    p[i]=1-Sy(ax[i],y[i],ymax,b,hfun)
    f[i]=h(y[i],ax[i],b)*(1-p[i])
  }
  #  return(list(f=f,p=p,p0=p0))
  return(f*intenscale)
}


b=log(c(0.5,0.2))
nf=100
ymax=50
ys=seq(ymin,ymax,length=nf)
f=fy(ys,rep(0,nf),ymin,ymax,b,h1)
plot(ys,f,type="l")
#plot(ys,fp$p,type="l")

ymax=100
ymin=ymax/1e9
ys=seq(ymin,ymax,length=nf)
f=fy(ys,rep(0,nf),ymin,ymax,b=0.1,h.const)
plot(ys,f,type="l")
#plot(ys,fp$p,type="l")

b=log(c(0.75,1))
ymax=5
ymin=0.01
ys=seq(ymin,ymax,length=nf)
plot(ys,h2(ys,rep(0,nf),b),type="l")
f=fy(ys,rep(0,nf),ymin,ymax,b,h2)
plot(ys,f,type="l")
#plot(ys,fp$p,type="l")

ymax=5
W=2
window=owin(c(0,ymax),c(-W,W))
isc=500
f=fy(ys,rep(0,nf),ymin,ymax,h2,b,intenscale=isc)
plot(ys,f,type="l")
lmax=max(f)
pp=rpoispp(fy,lmax,window,ymin=ymin,ymax=ymax,hfun=h2,b=b,intenscale=isc)
quartz(h=6,w=8)
plot(density(pp))
contour(density(pp),add=TRUE)
plot(pp,pch="+",cex=0.75,add=TRUE)
hist(abs(pp$y),nclass=35)
hist(pp$x,nclass=35)

ymax=5
W=2
window=owin(c(0,ymax),c(-W,W))
isc=500
f=fyx(ys,rep(0,nf),ymin=ymin,ymax=ymax,b=b,hfun=h2,pi.x=pi.norm,logphi=c(0.5,log(0.3)),W=W,intenscale=isc)
plot(ys,f,type="l")
lmax=max(f)
pp=rpoispp(fyx,lmax,window,ymin=ymin,ymax=ymax,b=b,hfun=h2,pi.x=pi.norm,logphi=c(0.5,log(0.3)),W=W,intenscale=isc)
quartz(h=6,w=8)
plot(density(pp))
contour(density(pp),add=TRUE)
plot(pp,pch="+",cex=0.75,add=TRUE)
hist(abs(pp$y),nclass=35)
hist(pp$x,nclass=35)

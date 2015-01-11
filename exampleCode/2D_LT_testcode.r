

source("~/dropbox/packages/2D distance sampling with time/R/2D_LT_functions V2a1.r")
#source("/Users/david/Research/2DLT/2D_LT_functions V2a1.r")
#source("/Users/david/Research/2DLT/2D_LT_functions V2a.r")
#source("C:\\DLB\\Research\\2DLT\\2D_LT_functions V2.r")
#source("C:\\DLB\\SoftwareDev\\Misc R functions\\utils.r")

#simdat=read.table("C:\\DLB\\Research\\2DLT\\simdat.csv")
#x=simdat$x
#y=simdat$y
#y[y==0]=1e-10
## set x and y limits:
#w=ceiling(max(x))
#ystart=1000

ystart=3
w=1
# grid for various plots:
gridx=seq(w/100,w,length=100)
gridy=seq(ystart/100,ystart,length=100)

# animal density function specification
#pi.x=pi.norm; logphi=c(10,log(400))
#pi.x=pi.norm; logphi=c(60,log(50))
mu.pi=0.5
sd.pi=0.2
#pi.x=pi.norm; logphi=c(mu.pi,log(sd.pi))
#pi.x=pi.const; logphi=c(1,NA)
pi.x=pi.hr2; logphi=log(c(0.75,1))
# plot it
adbn=pi.x(gridx,logphi,w)
plot(gridx,adbn,type="l",ylim=c(0,max(adbn)),xlab="perp dist (x)",ylab="pi(x)")

# hazard function specification
#hr=h1; b=log(c(0.75,0.2))
#hr=h1; b=log(c(0.001,1))
#hr=h.okamura; b=log(c(50,25))
#hr=h2; b=log(c(0.75,1)) # this produces a sensible h-r p(x) shape
hr=h.const; b=c(1,NA)
# check hazard function:
#hr(100,100,b)
# plot hazard function:
#windows(w=18); par(mfrow=c(1,3))
quartz(h=6,w=10); par(mfrow=c(1,3))
haz=outer(gridx,gridy,FUN=hr,b=b)
persp(gridx,gridy,log(haz),theta=45,phi=35)#,zlim=c(0,max(haz)))
range(haz)

# test fyx()
f=outer(gridy,gridx,FUN=fyx,b=b,hr=hr,ystart=ystart)
persp(gridx,gridy,t(f),theta=45,phi=35,zlim=c(0,max(f)),zlab="f(y|x)")

# test p.x()=1-S() - survival function, forward distance, y, is equivalent to time 
# - when moving at a constant speed
p.x=px(gridx,b,hr,ystart,nint=100)
plot(gridx,p.x,type="l",ylim=c(0,max(p.x)),xlab="x",ylab="p(x)")

# simulate detections and y-data for N animals:
# draw N perp. distances from adbn:
N=50
#u=runif(N,min=pnorm(0,mu.pi,sd.pi),max=pnorm(w,mu.pi,sd.pi))

plot(seq(0,1,length=1000),
     pi.hr2(x=seq(0,1,length=1000),logphi=logphi,w=w),
     type='l',
     xlab='perp. distance (x)',ylab='pi(x)')

hist(sample(x=seq(0,1,length=1000),size=1e5,replace=TRUE,
       prob=pi.x(x=seq(0,1,length=1000),logphi=logphi,w=w)),
     xlab='perp. distance (x)',main='Sampling obs. from pi(x)')
x=sample(x=seq(0,1,length=1000),size=N,replace=TRUE,
         prob=pi.hr2(x=seq(0,1,length=1000),logphi=logphi,w=w))
####

####
#x=qnorm(u,mu.pi,sd.pi)
y=sim.nhPP(x,b,ystart,hr)
keep=which(y>=0)
n=length(keep)
x=x[keep]
y=y[keep]
par(mfrow=c(2,2))
plot(x,y,ylim=c(0,ystart),cex=0.6)
rug(x,ticksize=0.1)
hist(y,freq=FALSE,xlab="forward dist. (y)")
hist(x,freq=FALSE,xlab="perp. dist. (x)")

# fit to y and x data, taking pi.x as known:
pars=b[1]
negloglik.yx.w(y,x,pars,hr,ystart,pi.x,logphi,w)

#library(profr)
est.yx.w=fityx.w(y,x,pars,hr,ystart,pi.x,logphi,w,control=list(trace=5))


# plot fitted functions
par(mfrow=c(2,2))
par(mfrow=c(1,1))
# perp dist dimension:
#  plot fitted lines
#est.yx.w$b=1
plotdat.yx.w=plotfit.x(x,est.yx)#.w)
#  add true lines:
hr=h.const
b=c(1,NA)  
p.xpi=p.pi.x(gridx,b,hr,ystart,pi.x,logphi,w)
mu=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)$value 
f.x=p.xpi/mu
lines(gridx,f.x,col="green")
ptot=integrate(f=px,lower=0,upper=w,b=b,hr=hr,ystart=ystart)$value
p.x.std=p.x/ptot
lines(gridx,p.x.std,col="green",lty=2)
#lines(gridx,p.x*plotdat.yx.w$p.xfit.std,col="green",lty=3)
lines(gridx,adbn,col="green",lty=3)
legend("topleft",title="True",legend=c("f(x)","p(x)","pi(x)"),col=c("green","green","green"),lwd=c(2,1,1),lty=c(1,2,3))

# forward dist dimension:
plotdat.yx.w.y=plotfit.yx(y,x,est.yx.w,nclass=4)
temp=est.yx.w; temp$b=b
junk=plotfit.yx(y,x,temp,lineonly=TRUE,col="green",lty=2)
# point estimate of abundance
Nhat.yx.w=n/plotdat.yx.w$mufit
cat("N=",N,"; n=",n,"; Nhat.yx.w=",signif(Nhat.yx.w,3),";``bias''=",signif((Nhat.yx.w/N-1)*100,3),"%\n",sep="")

# try fitting to x-data only:
#pi.x=pi.norm; logphi=c(mu.pi,log(sd.pi))
pars=b
negloglik.x(x,pars,hr,ystart,pi.x,logphi,w)
est.x=fitx(x,pars,hr,ystart,pi.x,logphi,w,control=list(trace=5))

# plot fitted functions
plotdat.x=plotfit.x(x,est.x)
lines(gridx,f.x,col="green")
lines(gridx,p.x*plotdat.x$p.xfit.std,col="green",lty=3)
lines(gridx,adbn,col="green",lty=3)
# point estimate of abundance
Nhat.x=n/plotdat.x$mufit
cat("N=",N,"; n=",n,"; Nhat.x=",signif(Nhat.x,3),";``bias''=",signif((Nhat.x/N-1)*100,3),"%\n",sep="")







for(i in 1:10) {

N=50
u=runif(N,min=pnorm(0,mu.pi,sd.pi),max=pnorm(w,mu.pi,sd.pi))
x=qnorm(u,mu.pi,sd.pi)
y=sim.nhPP(x,b,ystart,hr)
keep=which(y>=0)
n=length(keep)
x=x[keep]
y=y[keep]

# fit to y and x data, taking pi.x as UNknown:
  pars=c(logphi,logphi)
negloglik.yx(y,x,pars,hr,ystart,pi.x,w)
est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w,control=list(trace=5))
# plot fitted functions
#windows() # uncomment this line on a PC
#dev.new() # uncomment this line on a Mac
quartz()
par(mfrow=c(2,1))

# perp dist dimension:
est.yx
plotdat.yx=plotfit.x(x,est.yx)


p.xpi=p.pi.x(gridx,b,hr,ystart,pi.x,logphi,w)
mu=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)$value 
f.x=p.xpi/mu
lines(gridx,f.x,col="green",lwd=2)
lines(gridx,p.x*plotdat.yx$p.xfit.std,col="green",lty=2,lwd=2)
lines(gridx,adbn,col="green",lty=3,lwd=2)
legend("topleft",title="True",legend=c("f(x)","p(x)","pi(x)"),col=c("green","green","green"),lwd=c(2,1,1),lty=c(1,2,3))

# forward dist dimension:
plotdat.yx.y=plotfit.yx(y,x,est.yx,nclass=4)
temp=est.yx; temp$b=b
junk=plotfit.yx(y,x,temp,lineonly=TRUE,col="green",lty=2)
# point estimate of abundance
Nhat.yx=n/plotdat.yx$mufit
cat("N=",N,"; n=",n,"; Nhat.yx=",signif(Nhat.yx,3),";``bias''=",signif((Nhat.yx/N-1)*100,3),"%\n",sep="")

}
fn=function(x,pars) exp(-(x-pars[1])/pars[2])
lines(seq(0.05,1,0.05), fn(seq(0.05,1,0.05),c(0.4,0.1)),col='green')
# Progress: 
# 6/8/10: 
# -------
# 1. Use expected lifetime CONDITIONAL ON LIFETIME BEING LESS THAN (ymax-0) to get GoF in y-dimension.
#    - Q-Q plot

######
#try estimating pi(x) and p(x) parameters using known forms of p(x) uniform and 
#pi(x) - hr2 with starting pars

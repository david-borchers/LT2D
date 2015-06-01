#22 May
library(xlsx)
library(Distance)
dat=read.xlsx("/Users/dlb/Research/2DLT/Jantho Primate Line.xlsx",1)
x=dat$PP.Distance
y=dat$Forward.Distance
quartz(h=3,w=12)
par(mfrow=c(1,3))
pdlab="Perpendicular distance"
fdlab="Distance along transect"
sp=unique(dat$species)
plot(y,x,pch="+",ylab=pdlab,xlab=fdlab,main="",col=sp)
legend("topright",legend=sp,col=sp,pch="+",title="species",cex=0.75)
hist(y,breaks=seq(0,max(na.omit(y)),length=16),xlab=fdlab,main="")
hist(x,breaks=seq(0,max(na.omit(x)),length=12),xlab=pdlab,main="")

plot(jitter(y,1,0),jitter(x),pch="+",ylab=pdlab,xlab=fdlab,main="",col=sp)

nas=which(is.na(y))
x=x[-nas]
y=y[-nas]
w=0.03;ystart=0.05
# some stuff to get reasonable start values:
logphi=c(0.015,log(0.01))
mu=logphi[1]
sigma=exp(logphi[2])
f=dnorm(x,mean=mu,sd=sigma)/(pnorm(w,mean=mu,sd=sigma)-pnorm(0,mean=mu,sd=sigma))
hist(x[x<=w],freq=FALSE)
points(x[x<=w],f[x<=w])
b=log(c(0.02,1.75))
hr1.to.p(seq(0,0.03,length=20),b=b,w)
h=hr1.to.p(x,b,w)
plot(x[x<=w],h[x<=w],pch="+",ylim=c(0,1))

# Try fitting a few models:
fit.n=fityx(y[x<=w],x[x<=w],b=b,hr=h1,ystart=ystart,
            pi.x=pi.norm,logphi=logphi,w=w,hessian=FALSE,control=list(trace=5))
plotfit.x(x,fit.n)
fit.unif=fityx(y[x<=w],x[x<=w],b=b,hr=h1,ystart=ystart,
               pi.x=pi.const,logphi=logphi,w=w,hessian=FALSE,control=list(trace=5))
plotfit.x(x,fit.unif)
aics=c(fit.n$AIC,fit.unif$AIC)
names(aics)=c("fit.n$AIC","fit.unif$AIC")
daics=aics-min(aics)
sort(daics)

#27 May
#————————————— Code —————————————
library(xlsx)
library(boot) # just for logit and inv.logit functions
source('~/Documents/packages/LT2D/R/2DLTfunctions.r')

dat=read.xlsx("~/Dropbox/packages/LT2D/data/Jantho Primate Line.xlsx",1)
x=dat$PP.Distance
y=dat$Forward.Distance
nas=which(is.na(y))
x=x[-nas]
y=y[-nas]

quartz(h=3,w=12)
par(mfrow=c(1,3))
pdlab="Perpendicular distance"
fdlab="Distance along transect"
sp=unique(dat$species)
plot(y,x,pch="+",ylab=pdlab,xlab=fdlab,main="",col=sp)
legend("topright",legend=sp,col=sp,pch="+",title="species",cex=0.75)
hist(y,breaks=seq(0,max(na.omit(y)),length=16),xlab=fdlab,main="")
hist(x,breaks=seq(0,max(na.omit(x)),length=12),xlab=pdlab,main="")

plot(jitter(y,1,0),jitter(x),pch="+",ylab=pdlab,xlab=fdlab,main="",col=sp)

w=0.03;ystart=0.05
#y=y[x<=w]
#x=x[x<=w]
# some stuff to get reasonable start values:
logphi=c(-3,350)
f=pi.invlogit(x,logphi,w)
plot(x[x<=w],f[x<=w])
mu=logphi[1]
sigma=exp(logphi[2])
f=dnorm(x,mean=mu,sd=sigma)/(pnorm(w,mean=mu,sd=sigma)-pnorm(0,mean=mu,sd=sigma))
hist(x[x<=w],freq=FALSE)
points(x[x<=w],f[x<=w])
b=log(c(0.02,1.75))
hr1.to.p(seq(0,0.03,length=20),b=b,w)
h=hr1.to.p(x,b,w)
plot(x[x<=w],h[x<=w],pch="+",ylim=c(0,1))

# Try fitting a few models:
w=0.03;ystart=0.05
# Normal bump:
b=c(-7.3287948, 0.9945317)
logphi=c(.01646734, -4.67131261)
fit.n=fityx(y[x<=w],x[x<=w],b=b,hr=h1,ystart=ystart,
            pi.x=pi.norm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5))
plotfit.x(x,fit.n)
plotfit.y(y,x,fit.n,nclass=13)
# Uniform
b=c(-7.0744154,0.9876447)
logpi=c(1,1)
fit.unif=fityx(y[x<=w],x[x<=w],b=b,hr=h1,ystart=ystart,
               pi.x=pi.const,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5))
plotfit.x(x,fit.unif)
plotfit.y(y,x,fit.unif,nclass=13)
# Normal dip:
b=c(-7.3329141,0.9948721)
logphi=c(-0.05,-4.7)
fit.chn=fityx(y[x<=w],x[x<=w],b=b,hr=h1,ystart=ystart,
              pi.x=pi.chnorm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5))
plotfit.x(x,fit.chn)
plotfit.y(y,x,fit.chn,nclass=13)

# Normal bump with IP1 hazard function:
b=c(logit(0.5),log(c(0.0005,1)))
logphi=c(.01646734, -4.67131261)
# ests from initial fit s
b=c(-23.694670, -3.143091, 2.120717)
logphi=c(0.01808988, -4.41129077)
fit.ng0=fityx(y[x<=w],x[x<=w],b=b,hr=ip1,ystart=ystart,
              pi.x=pi.norm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5,maxit=1000))
plotfit.x(x,fit.ng0)
plotfit.y(y,x,fit.ng0,nclass=13)
# g0:
1-Sy(0,0,ystart,fit.ng0$b,ip1)
plotfit.y(x=rep(0,length(y)),est=fit.ng0,lineonly=TRUE) # f(y|x=0)
yy=seq(0,ystart,length=100)
plot(yy,1-Sy(rep(0,length(yy)),yy,ystart,fit.ng0$b,ip1),type="l",
     xlab="Forward distance (y)",ylab="p(detect by y)") # p(y|x=0)

# Look at estimates of mean p:
phat(fit.n);phat(fit.unif);phat(fit.chn);phat(fit.ng0)
phatInterval(fit.n)
phatInterval(fit.unif)
phatInterval(fit.chn)
phatInterval(fit.ng0)

# Look at AICs
aics=c(fit.n$AIC,fit.unif$AIC,fit.chn$AIC,fit.ng0$AIC)
names(aics)=c("fit.n$AIC","fit.unif$AIC","fit.chn$AIC","fit.ng0$AIC")
daics=aics-min(aics)
sort(daics)

#simulator
simoutFn='~/dropbox/packages/LT2D/simulations/20150530/avoidanceN200.csv'
#voidanceN200.csv -
N=200
nSim=1000
for(i in 1:nSim)
{ #start simulation loop
sim=simXY(N=N,pi.x=match.fun(fit.ng0$pi.x),logphi=fit.ng0$logphi,
      hr=fit.ng0$hr,b=fit.ng0$b,w=fit.ng0$w,
      ystart=fit.ng0$ystart,xSampL=1e3,discardNotSeen=TRUE)
hr=match.fun(sim$settings$hr)
w=sim$settings$w

fit=fityx(sim$locs$y,sim$locs$x,
          b=sim$settings$b,
          hr=hr,
          ystart=sim$settings$ystart,
              pi.x=pi.norm,
          logphi=sim$settings$logphi,
          w=sim$settings$w,
          hessian=TRUE,
          control=list(trace=5,maxit=1000))
res=phatInterval(fit)
res$N=N
res$n=length(sim$locs$y)
res$Nhat2DLT=n/res$phat
res$relBias2DLT=(res$Nhat2DLT-res$N)/res$N
res$NhatLow2DLT=n/res$upper.bound
res$NhatHigh2DLT=n/res$lower.bound
res$covered2DLT=TRUE
res$covered2DLT[res$NhatLow2DLT>res$N]=FALSE
res$covered2DLT[res$NhatHigh2DLT<res$N]=FALSE

# Fit using CDS:
# Must have column names "Region.Label", "Area", "Sample.Label", "Effort", "distance", "size"
n.all=length(sim$locs$x)
ones=rep(1,n.all)
ddat=data.frame(Region.Label=ones, Area=ones, Sample.Label=ones, Effort=ones, 
                distance=sim$locs$x, size=ones)
fit.hr<-ds(ddat[ddat$distance<=w,],key="hr",adjustment=NULL)

summary(fit.hr)
gofCDS=ddf.gof(fit.hr$ddf)$dsgof

fit.hr.sum=summary(fit.hr$ddf)
res$phatCDS=fit.hr.sum$average.p
res$phatCVCDS=fit.hr.sum$average.p.se/fit.hr.sum$average.p
res$NhatCDS=fit.hr$ddf$Nhat # abundance in covered region
res$NCVCDS=as.vector(fit.hr.sum$Nhat.se/fit.hr$ddf$Nhat)
res$NhatLowCDS=as.vector(res$NhatCDS-1.96*fit.hr.sum$Nhat.se)
res$NhatHighCDS=as.vector(res$NhatCDS+1.96*fit.hr.sum$Nhat.se)
res$coveredCDS=TRUE
res$coveredCDS[res$NhatLowCDS<res$N]=FALSE
res$coveredCDS[res$NhatHighCDS>res$N]=FALSE
res$relBiasCDS=(res$NhatCDS-res$N)/res$N
res$gofkspCDS=gofCDS$ks$p
res$gofCvMpCDS=gofCDS$CvM$p

par=fit$par
names(par)=paste('parEst',1:length(par),sep='')
parCV=fit$CVpar
names(parCV)=paste('parEstCV',1:length(par),sep='')
res=cbind.data.frame(res,t(par),t(parCV))
if(i==1)
  out=res else
  out=rbind.data.frame(out,res)  
  write.csv(out,simoutFn)
}#end of simulation loop.
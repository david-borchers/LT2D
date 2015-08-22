library(LT2D)
library(xlsx)
library(goftest)

# Primate data:
# =============
dat=read.xlsx("./data/Jantho Primate Line.xlsx",1)
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
#plot(jitter(y,1,0),jitter(x),pch="+",ylab=pdlab,xlab=fdlab,main="",col=sp)
plot(jitter(y,1,0),jitter(x),pch="+",ylab=pdlab,xlab=fdlab,main="")
#legend("topright",legend=sp,col=sp,pch="+",title="species",cex=0.75)
hist(y,breaks=seq(0,max(na.omit(y)),length=16),xlab=fdlab,main="")
hist(x,breaks=seq(0,max(na.omit(x)),length=12),xlab=pdlab,main="")

# Try fitting a few models:
w=0.03;ystart=0.05

# Normal bump with hazard h1:
b=c(-7.3287948, 0.9945317)
logphi=c(.01646734, -4.67131261)
fit.n=fityx(y[x<=w],x[x<=w],b=b,hr=h1,ystart=ystart,
            pi.x=pi.norm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5))
plotfit.x(x[x<=w],fit.n,nclass=20);rug(x[x<=w])
GoFx(fit.n,plot=TRUE)$pvals
#plotfit.y(y[x<=w],x,fit.n,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fit.n,nclass=32);rug(x=y[x<=w])
GoFy(fit.n,plot=TRUE)$pvals
#EHSW:
phatInterval(fit.n)
phatInterval(fit.n)*w
# p(0):
p0.n=1-Sy(0,0,ystart,fit.n$b,h1);p0.n
plotfit.smoothfy(fit.n,xmax=0.004)

# Uniform with h1
b=c(-7.0744154,0.9876447)
logpi=c(1,1)
fit.unif=fityx(y[x<=w],x[x<=w],b=b,hr=h1,ystart=ystart,
               pi.x=pi.const,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5))
plotfit.x(x[x<=w],fit.unif,nclass=20);rug(x[x<=w])
GoFx(fit.unif,plot=TRUE)$pvals
#plotfit.y(y[x<=w],x,fit.unif,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fit.unif,nclass=32);rug(x=y[x<=w])
GoFy(fit.unif,plot=TRUE)$pvals
#EHSW:
phatInterval(fit.unif)
phatInterval(fit.unif)*w
# p(0):
p0.unif=1-Sy(0,0,ystart,fit.unif$b,h1);p0.unif
plotfit.smoothfy(fit.unif,xmax=0.004)

# Normal dip:
b=c(-7.3329141,0.9948721)
logphi=c(-0.05,-4.7)
fit.chn=fityx(y[x<=w],x[x<=w],b=b,hr=h1,ystart=ystart,
              pi.x=pi.chnorm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5))
plotfit.x(x[x<=w],fit.chn,nclass=20);rug(x[x<=w])
GoFx(fit.chn,plot=TRUE)$pvals
#plotfit.y(y[x<=w],x,fit.chn,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fit.chn,nclass=20);rug(x=y[x<=w])
GoFy(fit.chn,plot=TRUE)$pvals
#EHSW:
phatInterval(fit.chn)
phatInterval(fit.chn)*w
# p(0):
p0.chn=1-Sy(0,0,ystart,fit.chn$b,h1);p0.chn
plotfit.smoothfy(fit.chn,xmax=0.004)


# ip1 with normal bump:
b=c(5.2919208, -0.2205593, 8.4701307)
logphi=c(0.01784102, -4.42209067)
fit.n.ip1=fityx(y[x<=w],x[x<=w],b=b,hr=ip1,ystart=ystart,
                pi.x=pi.norm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5,maxit=1000))
plotfit.x(x[x<=w],fit.n.ip1,nclass=20);rug(x[x<=w])
GoFx(fit.n.ip1,plot=TRUE)$pvals
#plotfit.y(y[x<=w],x,fit.n.ip1,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fit.n.ip1,nclass=32);rug(x=y[x<=w])
GoFy(fit.n.ip1,plot=TRUE)$pvals
#EHSW:
phatInterval(fit.n.ip1)
phatInterval(fit.n.ip1)*w
# p(0):
p0.ip1=1-Sy(0,0,ystart,fit.n.ip1$b,ip1);p0.ip1
plotfit.smoothfy(fit.n.ip1,xmax=0.004)


# Normal bump with ip0 hazard function (SELECTED MODEL):
b=c(5.2919208, 8.4701307)
logphi=c(0.01784102, -4.42209067)
fit.n.ip0=fityx(y[x<=w],x[x<=w],b=b,hr=ip0,ystart=ystart,
                pi.x=pi.norm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5,maxit=1000))
quartz();par(mfrow=c(2,2))
plotfit.x(x[x<=w],fit.n.ip0,nclass=20);rug(x[x<=w])
GoFx(fit.n.ip0,plot=TRUE)$pvals
#plotfit.y(y[x<=w],x,fit.n.ip0,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fit.n.ip0,nclass=32);rug(x=y[x<=w])
GoFy(fit.n.ip0,plot=TRUE)$pvals
#EHSW:
phatInterval(fit.n.ip0)
phatInterval(fit.n.ip0)*w
# p(0):
p0.ip0=1-Sy(0,0,ystart,fit.n.ip0$b,ip0);p0.ip0
plotfit.smoothfy(fit.n.ip0,xmax=0.004)


# Normal bump with EP1 hazard function:
b=c(1.555219, -3.698782, 1.330832)
logphi=c(0.01575465, -4.45259285)
fit.n.ep1=fityx(y[x<=w],x[x<=w],b=b,hr=ep1,ystart=ystart,
                pi.x=pi.norm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5,maxit=1000))
plotfit.x(x[x<=w],fit.n.ep1,nclass=20);rug(x[x<=w])
GoFx(fit.n.ep1,plot=TRUE)$pvals
#plotfit.y(y[x<=w],x,fit.n.ep1,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fit.n.ep1,nclass=32);rug(x=y[x<=w])
GoFy(fit.n.ep1,plot=TRUE)$pvals
#EHSW:
phatInterval(fit.n.ep1)
phatInterval(fit.n.ep1)*w
# p(0):
p0.ep1=1-Sy(0,0,ystart,fit.n.ep1$b,ep1);p0.ep1
plotfit.smoothfy(fit.n.ep1,xmax=0.004)


# Genralized hazard rate with half-normal bump:
b2=c(fit.n$b,-10) # exact same model but in ghy terms
b2=c(-19.42374920,2.05969961,-3.36703594)
logphi=fit.n$logphi
logphi=c(,0.01230915,-4.71841683)
fit.ghy=fityx(y[x<=w],x[x<=w],b=b2,hr=ghy,ystart=ystart,
              pi.x=pi.norm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5,maxit=2000))
plotfit.x(x[x<=w],fit.ghy,nclass=20);rug(x[x<=w])
GoFx(fit.ghy,plot=TRUE)$pvals
#plotfit.y(y[x<=w],x,fit.ghy,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fit.ghy,nclass=32);rug(x=y[x<=w])
GoFy(fit.ghy,plot=TRUE)$pvals
#EHSW:
phatInterval(fit.ghy)
phatInterval(fit.ghy)*w
# p(0):
p0.ghy=1-Sy(0,0,ystart,fit.ghy$b,ghy);p0.ghy
plotfit.smoothfy(fit.ghy,xmax=0.004)


# Genralized hazard rate with 2 scale parameters and half-normal bump:
b3=c(fit.ghy$b,fit.ghy$b[1]) # exact same model but in ghy terms
pars=c(b3,fit.ghy$logphi)
b=c(-36.904485, 2.513848 ,-2.661960,-25.971073)
logphi=c(0.01919996,-4.26128280)
fit.ghy2=fityx(y[x<=w],x[x<=w],b=b3,hr=ghy2,ystart=ystart,
               pi.x=pi.norm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5,maxit=1000))
plotfit.x(x[x<=w],fit.ghy2,nclass=20);rug(x[x<=w])
GoFx(fit.ghy2,plot=TRUE)$pvals
#plotfit.y(y[x<=w],x,fit.ghy2,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fit.ghy2,nclass=32);rug(x=y[x<=w])
GoFy(fit.ghy2,plot=TRUE)$pvals
#EHSW:
phatInterval(fit.ghy2)
phatInterval(fit.ghy2)*w
# p(0):
p0.ghy2=1-Sy(0,0,ystart,fit.ghy2$b,ghy2);p0.ghy2
plotfit.smoothfy(fit.ghy2,xmax=0.004)


# Look at AICs 
aics=c(fit.n$AIC,fit.unif$AIC,fit.chn$AIC,fit.n.ip1$AIC,fit.n.ep1$AIC,fit.n.ip0$AIC,fit.ghy$AIC,fit.ghy2$AIC)
names(aics)=c("fit.n","fit.unif","fit.chn","fit.n.ip1","fit.n.ep1","fit.n.ip0","fit.ghy","fit.ghy2")
daics=aics-min(aics)
sort(daics)

# Fit using CDS:
library(Distance)
# Must have column names "Region.Label", "Area", "Sample.Label", "Effort", "distance", "size"
n.all=length(x)
ones=rep(1,n.all)
ddat=data.frame(Region.Label=ones, Area=ones, Sample.Label=ones, Effort=ones, distance=x, size=ones)
fit.hr<-ds(ddat[ddat$distance<=w,],key="hr",adjustment=NULL)
fit.hn<-ds(ddat[ddat$distance<=w,],key="hn",adjustment=NULL)
summary(fit.hr)
# Use HR for its shoulder:
plot(fit.hr,nc=20,showpoints=FALSE,pl.den=0)
ddf.gof(fit.hr$ddf)
# extract ESHW:
fit.hr.summ=summary(fit.hr)
eshw=fit.hr.summ$ds$average.p*w
cv.eshw=fit.hr.summ$ds$average.p.se/fit.hr.summ$ds$average.p
eshw;cv.eshw

fit.hr$ddf$Nhat # abundance in covered region
n=length(ddat$distance[ddat$distance<=w])
# Compare with same model fit, but using forward distances:
fit.hr$ddf$Nhat
n/phat(fit.unif)
(n/phat(fit.unif))/fit.hr$ddf$Nhat
# Compare with best model fit using forward distances:
n/phat(fit.n.ip1)
(n/phat(fit.n.ip1))/fit.hr$ddf$Nhat
# same, but if believe that p(0)=1:
(p0.n.ip1*n/phat(fit.n.ip1))/fit.hr$ddf$Nhat
# Compare CVs (0.058 for CDS)
phatInterval(fit.n.ip1)



# Dolphin data:
# ============
adat=read.xlsx("./data/AnasData.xlsx",1)
rdist=adat$Radial.distance
ang=adat$Angle
xd=abs(rdist*sin(ang*pi/180))
yd=abs(rdist*cos(ang*pi/180))
size=adat$Cluster.size
wd=0.15;ystartd=0.55
quartz(h=3,w=12)
par(mfrow=c(1,3))
pdlab="Perpendicular distance"
fdlab="Forward distance"
#plot(jitter(yd[xd<=wd],1,0),jitter(xd[xd<=wd],1,0),pch="+",ylab=pdlab,xlab=fdlab,main="",cex=round(log(size[xd<=wd])))
plot(jitter(yd[xd<=wd],1,0),jitter(xd[xd<=wd],1,0),pch="+",ylab=pdlab,xlab=fdlab,main="")
hist(yd[xd<=wd],breaks=seq(0,max(yd[xd<=wd]),length=26),xlab=fdlab,main="")
hist(xd[xd<=wd],breaks=seq(0,max(xd[xd<=wd]),length=19),xlab=pdlab,main="")


# Hazard-rate with half-normal bump: (SELECTED MODEL)
b=c(-7.3287948, 0.9945317)
logphi=-0.4811025
dfit.hn=fityx(yd[xd<=wd],xd[xd<=wd],b=b,hr=h1,ystart=ystartd,
              pi.x=pi.hnorm,logphi=logphi,w=wd,hessian=TRUE,control=list(trace=5))
quartz();par(mfrow=c(2,2))
plotfit.x(xd[xd<=wd],dfit.hn,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.hn,plot=TRUE)$pvals
#plotfit.y(yd[xd<=wd],xd,dfit.hn,nclass=20);rug(x=yd[xd<=wd])
plotfit.smoothfy(dfit.hn,nclass=32);rug(x=yd[xd<=wd])
GoFy(dfit.hn,plot=TRUE)$pvals
#EHSW:
phatInterval(dfit.hn)
phatInterval(dfit.hn)*wd
# p(0):
p0.hn=1-Sy(0,0,ystartd,dfit.hn$b,h1);p0.hn
plotfit.smoothfy(dfit.hn,xmax=0.01)
# Density estimate:
n=length(dfit.hn$dat$x)
L=1672.77 # from Canadas et al. (in nm)
Dhat=(n/phatInterval(dfit.hn)$phat)/(2*wd*L)
Dhat

# Generalized hazard-rate (ghy) with half-normal bump:
b=c(-7.3287948, 0.9945317, log(0.025))
logphi=-0.4811025
dfit.hn.ghy=fityx(yd[xd<=wd],xd[xd<=wd],b=b,hr=ghy,ystart=ystartd,
                  pi.x=pi.hnorm,logphi=logphi,w=wd,hessian=TRUE,control=list(trace=5))
plotfit.x(xd[xd<=wd],dfit.hn.ghy,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.hn.ghy,plot=TRUE)$pvals
plotfit.smoothfy(dfit.hn.ghy,nclass=32);rug(x=yd[xd<=wd])
#plotfit.y(yd[xd<=wd],xd,dfit.hn.ghy,nclass=20);rug(x=yd[xd<=wd])
GoFy(dfit.hn.ghy,plot=TRUE)$pvals
#EHSW:
phatInterval(dfit.hn.ghy)
phatInterval(dfit.hn.ghy)*wd
# p(0):
p0.hn.ghy=1-Sy(0,0,ystartd,dfit.hn.ghy$b,ghy);p0.hn.ghy
plotfit.smoothfy(dfit.hn.ghy,xmax=0.01)


# Hazard rate with uniform
b=c(-7.0744154,0.9876447)
logphi=NULL
dfit.unif=fityx(yd[xd<=wd],xd[xd<=wd],b=b,hr=h1,ystart=ystartd,
                pi.x=pi.const,logphi=logphi,w=wd,hessian=TRUE,control=list(trace=5))
plotfit.x(xd[xd<=wd],dfit.unif,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.unif,plot=TRUE)$pvals
#plotfit.y(yd[xd<=wd],xd,dfit.unif,nclass=20);rug(x=yd[xd<=wd])
plotfit.smoothfy(dfit.unif,nclass=32);rug(x=yd[xd<=wd])
GoFy(dfit.unif,plot=TRUE)$pvals
#EHSW:
phatInterval(dfit.unif)
phatInterval(dfit.unif)*wd
# p(0):
p0.unif=1-Sy(0,0,ystartd,dfit.unif$b,h1);p0.unif
plotfit.smoothfy(dfit.unif,xmax=0.01)

# ip1 bump with uniform :
b=c(4.0203610, -3.4051984, 0.6552985)
logphi=NULL
dfit.ip1.u=fityx(yd[xd<=wd],xd[xd<=wd],b=b,hr=ip1,ystart=ystartd,
                 pi.x=pi.const,logphi=logphi,w=wd,hessian=TRUE,control=list(trace=5,maxit=1000))
plotfit.x(xd[xd<=wd],dfit.ip1.u,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.ip1.u,plot=TRUE)$pvals
#plotfit.y(yd[xd<=wd],xd,dfit.ip1.u,nclass=20);rug(x=yd[xd<=wd])
plotfit.smoothfy(dfit.ip1.u,nclass=32);rug(x=yd[xd<=wd])
GoFy(dfit.ip1.u,plot=TRUE)$pvals
#EHSW:
phatInterval(dfit.ip1.u)
phatInterval(dfit.ip1.u)*wd
# p(0):
p0.ip1.u=1-Sy(0,0,ystartd,dfit.ip1.u$b,ip1);p0.ip1.u
plotfit.smoothfy(dfit.ip1.u,xmax=0.01)


# Half-normal bump with ip1:
b=c(4.0203610, -3.4051984, 0.6552985)
logphi=-2.279662
dfit.n.ip1=fityx(yd[xd<=wd],xd[xd<=wd],b=b,hr=ip1,ystart=ystartd,
                 pi.x=pi.hnorm,logphi=logphi,w=wd,hessian=TRUE,control=list(trace=5,maxit=1000))
plotfit.x(xd[xd<=wd],dfit.n.ip1,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.n.ip1,plot=TRUE)$pvals
#plotfit.y(yd[xd<=wd],xd,dfit.n.ip1,nclass=20);rug(x=yd[xd<=wd])
plotfit.smoothfy(dfit.n.ip1,nclass=32);rug(x=yd[xd<=wd])
GoFy(dfit.n.ip1,plot=TRUE)$pvals
#EHSW:
phatInterval(dfit.n.ip1)
phatInterval(dfit.n.ip1)*wd
# p(0):
p0.n.ip1=1-Sy(0,0,ystartd,dfit.n.ip1$b,ip1);p0.n.ip1
plotfit.smoothfy(dfit.n.ip1,xmax=0.01)


# Half-normal bump with ep1:
b=c(10.0078331, -2.9408222, -0.1744689)
logphi=-2.329975
dfit.n.ep1=fityx(yd[xd<=wd],xd[xd<=wd],b=b,hr=ep1,ystart=ystartd,
                 pi.x=pi.hnorm,logphi=logphi,w=wd,hessian=TRUE,control=list(trace=5,maxit=1000))
plotfit.x(xd[xd<=wd],dfit.n.ep1,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.n.ep1,plot=TRUE)$pvals
#plotfit.y(yd[xd<=wd],xd,dfit.n.ep1,nclass=20);rug(x=yd[xd<=wd])
plotfit.smoothfy(dfit.n.ep1,nclass=32);rug(x=yd[xd<=wd])
GoFy(dfit.n.ep1,plot=TRUE)$pvals
#EHSW:
phatInterval(dfit.n.ep1)
phatInterval(dfit.n.ep1)*wd
# p(0):
p0.n.ep1=1-Sy(0,0,ystartd,dfit.n.ep1$b,ep1);p0.n.ep1
plotfit.smoothfy(dfit.n.ep1,xmax=0.01)



# Half-normal bump with ip2:
b2=c(3.6737110, -3.352613, 0.8568386,-3.1155813)
logphi=-2.459988
dfit.n.ip2=fityx(yd[xd<=wd],xd[xd<=wd],b=b2,hr=ip2,ystart=ystartd,
                 pi.x=pi.hnorm,logphi=logphi,w=wd,hessian=TRUE,control=list(trace=5,maxit=1000))
plotfit.x(xd[xd<=wd],dfit.n.ip2,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.n.ip2,plot=TRUE)$pvals
#plotfit.y(yd[xd<=wd],xd,dfit.n.ip2,nclass=20);rug(x=yd[xd<=wd])
plotfit.smoothfy(dfit.n.ip2,nclass=32);rug(x=yd[xd<=wd])
GoFy(dfit.n.ip2,plot=TRUE)$pvals
#EHSW:
phatInterval(dfit.n.ip2)
phatInterval(dfit.n.ip2)*wd
# p(0):
p0.n.ip2=1-Sy(0,0,ystartd,dfit.n.ip2$b,ip2);p0.n.ip2
plotfit.smoothfy(dfit.n.ip2,xmax=0.01)


# Half-normal bump with ep2:
b=c(-2.7688390, -2.8147403, -0.2064417, -2.9908584)
logphi=-2.426298
dfit.n.ep2=fityx(yd[xd<=wd],xd[xd<=wd],b=b,hr=ep2,ystart=ystartd,
                 pi.x=pi.hnorm,logphi=logphi,w=wd,hessian=TRUE,control=list(trace=5,maxit=1000))
plotfit.x(xd[xd<=wd],dfit.n.ep2,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.n.ep2,plot=TRUE)$pvals
#plotfit.y(yd[xd<=wd],xd,dfit.n.ep2,nclass=20);rug(x=yd[xd<=wd])
plotfit.smoothfy(dfit.n.ep2,nclass=32);rug(x=yd[xd<=wd])
GoFy(dfit.n.ep2,plot=TRUE)$pvals
#EHSW:
phatInterval(dfit.n.ep2)
phatInterval(dfit.n.ep2)*wd
# p(0):
p0.n.ep2=1-Sy(0,0,ystartd,dfit.n.ep2$b,ep2);p0.n.ep2
plotfit.smoothfy(dfit.n.ep2,xmax=0.01)


library(Distance)
# Must have column names "Region.Label", "Area", "Sample.Label", "Effort", "distance", "size"
nd.all=length(xd)
ones=rep(1,nd.all)
dddat=data.frame(Region.Label=ones, Area=ones, Sample.Label=ones, Effort=ones, distance=xd, size=ones)
dfit.hn.d<-ds(dddat[dddat$distance<=wd,],key="hn",adjustment=NULL)
dfit.hr.d<-ds(dddat[dddat$distance<=wd,],key="hr",adjustment=NULL) # HR is best by AIC
n=length(dddat[dddat$distance<=wd,1])
summary(dfit.hr.d)
ddf.gof(dfit.hr.d$ddf)
plot(dfit.hr.d,nc=30,showpoints=FALSE,pl.den=0)
dfit.hr.d$ddf$Nhat
# extract ESHW:
dfit.hr.summ=summary(dfit.hr.d)
deshw=dfit.hr.summ$ds$average.p*wd
cv.deshw=dfit.hr.summ$ds$average.p.se/dfit.hr.summ$ds$average.p
deshw;cv.deshw


aics=c(dfit.hn$AIC,dfit.unif$AIC,dfit.n.ip1$AIC,dfit.n.ep1$AIC,dfit.n.ip2$AIC,dfit.n.ep2$AIC,dfit.ip1.u$AIC,dfit.hn.ghy$AIC)
names(aics)=c("dfit.hn","dfit.unif","dfit.n.ip1","dfit.n.ep1","dfit.n.ip2","dfit.n.ep2","dfit.ip1.u","dfit.hn.ghy")
daics=aics-min(aics)
sort(daics)

p0s=c(p0.hn,p0.unif,p0.n.ip1,p0.n.ep1,p0.n.ip2,p0.n.ep2,p0.ip1.u,p0.hn.ghy)
names(p0s)=c("dfit.hn","dfit.unif","dfit.n.ip1","dfit.n.ep1","dfit.n.ip2","dfit.n.ep2","dfit.ip1.u","dfit.hn.ghy")
dp0s=p0s-min(p0s)
sort(p0s)

# Density estimate:
#Nd=nd/(0.112/wd)
#Nd/(2*wd*1672.77)

# Primate and Dolphin plots for paper:
# ===================================

# Scatterplots
# ------------
pdlab="Perpendicular distance"
fdlab="Forward distance"
quartz(h=5,w=12);par(mfrow=c(1,2))
plot(jitter(y,1,0),jitter(x),pch="+",ylab=pdlab,xlab=fdlab,main="")
lines(c(-1,1.1*max(y)),rep(w,2),lty=2,col="gray")
lines(c(-1,1.1*max(y)),rep(xmax,2),lty=3,col="gray")
plot(jitter(yd,1,0),jitter(xd,1,0),pch="+",ylab=pdlab,xlab=fdlab,main="")
lines(c(-1,1.1*max(yd)),rep(wd,2),lty=2,col="gray")
lines(c(-1,1.1*max(yd)),rep(xdmax,2),lty=3,col="gray")

# Primate fits:
# -------------
quartz();par(mfrow=c(2,2))
plotfit.x(x[x<=w],fit.n.ip0,nclass=20);rug(x[x<=w])
GoFx(fit.n.ip0,plot=TRUE)$pvals
#plotfit.y(y[x<=w],x,fit.n.ip0,nclass=20);rug(x=y[x<=w])
plotfit.smoothfy(fit.n.ip0,nclass=32);rug(x=y[x<=w])
GoFy(fit.n.ip0,plot=TRUE)$pvals

# Dolphin fits:
# -------------
quartz();par(mfrow=c(2,2))
plotfit.x(xd[xd<=wd],dfit.hn.ghy,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.hn.ghy,plot=TRUE)$pvals
plotfit.smoothfy(dfit.hn.ghy,nclass=32);rug(x=yd[xd<=wd])
#plotfit.y(yd[xd<=wd],xd,dfit.hn.ghy,nclass=20);rug(x=yd[xd<=wd])
GoFy(dfit.hn.ghy,plot=TRUE)$pvals

# f(y) plots for x close to zero
# ------------------------------
pcin=0.11 # percentage of points to include
xdsort=sort(xd);edfxd=order(xdsort)/length(xd)
xdcut=length(edfxd[edfxd<=pcin])
xsort=sort(x);edfx=order(xsort)/length(x)
xcut=length(edf[edfx<=pcin])
xdmax=xdsort[xdcut]
xmax=xsort[xcut]
quartz(h=3.5);par(mfrow=c(1,2))
plotfit.smoothfy(fit.chn,xmax=xmax)
plotfit.smoothfy(dfit.hn,xmax=xdmax)

# CDS plots
# ---------
quartz(h=3.5);par(mfrow=c(1,2))
plot(fit.hr,nc=20,showpoints=FALSE,pl.den=0)
plot(dfit.hr.d,nc=14,showpoints=FALSE,pl.den=0)






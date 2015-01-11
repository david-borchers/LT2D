library(distance2Dt)

hr=h2; b=log(c(0.75,1))
pi.x=pi.norm; logphi=c(0.5,log(0.2))
N=50 #true number of animals
w=1;ystart=2;ymin=1e-8
#generate some observations
simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,hr=hr,b=b,w=w,ystart=ystart)
x=simDat$locs$x; y=simDat$locs$y
fityx(y,x,b,hr,ystart,pi.x,logphi,w)

quartz(h=5)
n=100000
simdat100K=sim.n(n,ymin,ystart,w,hr,b,pi.x,logphi,fix.n=TRUE,intscale=NULL,nbuffer=NULL)
dens100K=density(simdat100K$pp)
plot(dens100K,main="",ribbon=FALSE,asp=0.5)
n=200
simdat=sim.n(n,ymin,ystart,w,hr,b,pi.x,logphi,fix.n=TRUE,intscale=NULL,nbuffer=NULL)
plot(simdat$pp,pch="+",add=TRUE)
lines(c(0,ystart),rep(0,2))

# plot detections in plan view
plot(simdat$pp$x,simdat$pp$y,pch="+",asp=0.5,ylim=c(-w,w),xlim=c(0,ystart),
     xlab="Forward distance",ylab="Perpendicular distance",main="")
lines(c(0,ystart),rep(0,2))
points(0,0,pch=19)

# perp and forward dist histograms
quartz()
par(mfrow=c(2,1))
hist(simdat$locs$x,main="Perp dist dbn",xlab="Perpendicular distance")
hist(simdat$locs$y,main="Forward dist dbn",xlab="Forward distance")

# plot f(y|x=0)
ny=100
y=seq(0.001,ystart,length=ny)
fy0=fyx(y,rep(0,ny),b,hr,ystart)
quartz(h=4)
plot(y,fy0,type="l",xlab="Forward distance",ylab="f(y|x) at x=0")

# plot 2D hazard
gridx=seq(0,w,length=50); gridy=seq(0,ystart,length=50)
f=outer(gridy,gridx,FUN=fyx,b,hr,ystart)
persp(gridx,gridy,t(f),theta=45,phi=35,zlab="f(y|x)")

# plot f(y|x=0.75)
xx=3
yst=3
ny=100
y=seq(0.001,yst,length=ny)
fyxx=fyx(y,rep(xx,ny),b,hr,ystart=yst)
quartz(h=4)
plot(y,fyxx,type="l",ylim=c(0,max(fyxx)),xlab="Forward distance",ylab="f(y|x) at x=0.25")

# plot p(x)
nx=100
x=seq(0,w,length=nx)
pxs=px(x,b,hr,ystart)
quartz(h=4)
plot(x,pxs,type="l",xlab="Perpendicular distance",ylab="p(x)=1-S(T|x)")

# fit to data
x=simdat$locs$x; y=simdat$locs$y
fit=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
quartz(h=4)
plotfit.x(simdat$locs$x,fit,addTruth=TRUE,true.pi.x=pi.x,true.logphi=logphi,true.hr=hr,true.b=b)



# Try with attraction:

logphi=c(0,log(0.33))

quartz(h=5)
n=10000
bigsimdat=sim.n(n,ymin,ystart,w,hr,b,pi.x,logphi,fix.n=TRUE,intscale=NULL,nbuffer=NULL)
bigdens=density(bigsimdat$pp)
plot(bigdens,main="",ribbon=FALSE,asp=0.5)
n=200
simdat=sim.n(n,ymin,ystart,w,hr,b,pi.x,logphi,fix.n=TRUE,intscale=NULL,nbuffer=NULL)
plot(simdat$pp,pch="+",add=TRUE)
lines(c(0,ystart),rep(0,2))

plot(simdat$pp$x,simdat$pp$y,pch="+",asp=0.5,ylim=c(-w,w),xlim=c(0,ystart),
     xlab="Forward distance",ylab="Perpendicular distance",main="")
lines(c(0,ystart),rep(0,2))
points(0,0,pch=19)

quartz()
par(mfrow=c(2,1))
hist(simdat$locs$x,main="Perp dist dbn",xlab="Perpendicular distance")
hist(simdat$locs$y,main="Forward dist dbn",xlab="Forward distance")

quartz(h=4)
x=simdat$locs$x; y=simdat$locs$y
fit=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
plotfit.x(simdat$locs$x,fit,addTruth=TRUE,true.pi.x=pi.x,true.logphi=logphi,true.hr=hr,true.b=b)

# fit with bad start values
blogphi=c(0.2,log(0.5))
bfit=fityx(y,x,b,hr,ystart,pi.x,blogphi,w)
plotfit.x(simdat$locs$x,bfit,addTruth=TRUE,true.pi.x=pi.x,true.logphi=logphi,true.hr=hr,true.b=b)

## det fns look suspicious; plot separately
#x=seq(0,1,length=100)
#truep=px(x,b,hr,ystart)
#estp=px(x,fit$par[1:2],hr,ystart)
#bestp=px(x,bfit$par[1:2],hr,ystart)
#ylim=range(c(0,estp,bestp,truep))
#plot(x,truep,type="l")
#lines(x,estp,col="blue")
#lines(x,estp,col="red")

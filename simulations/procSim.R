#dat=read.csv('/Users/martincox/Downloads/simRes-g1dotNormN30(10).csv')
dat1=read.csv('/Users/martincox/Dropbox/packages/2D distance sampling with time/simulations/simRes-g1dotNorm0.5,log(0.2).csv')
dat2=read.csv('/Users/martincox/Dropbox/packages/2D distance sampling with time/simulations/simRes-g1dotNorm0.5,log(0.2)-2.csv')
dat=rbind.data.frame(dat1,dat2)
######CHANGE HERE FOR VARIABLE n #######
summary(dat$n)
cutV=seq(0,330,by=30)
dat$ng=dat$n #use this for fixed n
#dat$ng=cut(x=dat$n,breaks=cutV) #use for variable n
any(is.na(dat$ng))
nrow(dat)
table(dat$ng)
quartz('sim',3,5)
par(mar=c(2.5,3,1,1),mgp=c(1.5,0.25,0),tcl=-0.25,mfrow=c(3,1))
##panel 1
boxplot(phat.bias.percent~ng,data=dat,xlab='n',ylab='% bias')
title('% p hat relative bias')
abline(h=0,lty=3)
##panel 2
sapply(split(dat$phat.bias.percent,dat$ng),summary)
mV=with(dat,sapply(split(phat.bias.percent,ng),mean))
######CHANGE HERE FOR VARIABLE n #######
x=as.numeric(names(mV))
#x=cutV[-length(cutV)]+diff(cutV)[1]/2
##
xLim=c(min(x)-diff(x)[1]/2,max(x)+diff(x)[1]/2)
CI95=with(dat,sapply(split(phat.bias.percent,ng),function(x) sort(var(x)/length(x))))
lower=mV-CI95;upper=mV+CI95
yLim=c(min(lower),max(upper))
plot(x,mV,col=grey(0.2),type='o',xlim=xLim,ylim=yLim,
     xlab='n',ylab='Mean % bias')
title('Mean % p hat relative bias')
legend('bottomright',c('mean','95%CI'),lty=c(1,2),pch=c(1,NA))
abline(h=0,lty=3)
lines(x,lower,lty=2);lines(x,upper,lty=2)
##panel 3
cT=with(dat,table(covered,ng))
coverProb=data.frame(cT[1,]/apply(cT,2,sum))
names(coverProb)='coverp'
#plot(as.numeric(row.names(coverProb)),coverProb$coverp,type='o',pch=19,xlim=xLim,
plot(x,coverProb$coverp,type='o',pch=19,xlim=xLim,
  ylim=c(0,1.1*max(coverProb$coverp)),xlab='n',ylab='cover probability')
abline(h=c(0.05),lty=3)
title('Coverage prob. (prop. true p outside 95%CI)')
Qs=with(dat,sapply(split(phat,n),quantile,probs=c(0.025,0.975)))
pM=with(dat,sapply(split(phat,n),mean))
#########################
#this doesn't work for variable n:
require(plotrix)
par(mfrow=c(1,1))
#plotCI(x=as.numeric(colnames(Qs)),y=rep(unique(dat$p),ncol(Qs)),ui=Qs[2,],li=Qs[1,],
#       xlab='n',ylab='p hat')

plotCI(x=as.numeric(names(pM)),y=pM,ui=Qs[2,],li=Qs[1,],
       xlab='n',ylab='p hat')
abline(h=unique(dat$p),lty=3)
title('95%CI from non-parametric bootstrap')
########################

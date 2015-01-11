#source("2D_LT_functions V2a1.r")
source("~/dropbox/packages/2D distance sampling with time/R/2DLTfunctions.r")

# sims=list(g1dot0N60=list(ystart=3,w=1,
#                  pi.x=pi.norm, logphi=c(0.5,log(0.2)),
#                  hr=h2,b=log(c(0.75,1)),
#                  N=60,
#                  miss=FALSE,
#                  nSim=1e3),
#   g1dot0N200=list(ystart=3,w=1,
#                  pi.x=pi.norm, logphi=c(0.5,log(0.2)),
#                  hr=h2,b=log(c(0.75,1)),
#                  N=200,
#                  miss=FALSE,
#                  nSim=1e3))

# sims=list(g1dot0N60missing=list(ystart=3,w=1,
#                  pi.x=pi.norm, logphi=c(0.5,log(0.2)),
#                  hr=h2,b=log(c(0.75,1)),
#                  N=100,
#                  miss=TRUE,
#                  nSim=1e3),
#   g1dot0N200missing=list(ystart=3,w=1,
#                  pi.x=pi.norm, logphi=c(0.5,log(0.2)),
#                  hr=h2,b=log(c(0.75,1)),
#                  N=400,
#                  miss=TRUE,
#                  nSim=1e3),
#   g0dot8N60missing=list(ystart=3,w=1,
#                         pi.x=pi.norm, logphi=c(0.5,log(0.2)),
#                         hr=h21,b=c(log(c(0.75,1)),qlogis(0.8)),
#                         N=100,
#                         miss=TRUE,
#                         nSim=1e3),
#   g0dot8N200missing=list(ystart=3,w=1,
#                          pi.x=pi.norm, logphi=c(0.5,log(0.2)),
#                          hr=h21,b=c(log(c(0.75,1)),qlogis(0.8)),
#                          N=400,
#                          miss=TRUE,
#                          nSim=1e3)) 

sims=list(g1dot0N60missingHN=list(ystart=3,w=1,
                                pi.x=pi.norm, logphi=c(0,log(0.2)),
                                hr=h2,b=log(c(0.75,1)),
                                N=100,
                                miss=TRUE,
                                nSim=1e3),
          g1dot0N200missingHN=list(ystart=3,w=1,
                                 pi.x=pi.norm, logphi=c(0,log(0.2)),
                                 hr=h2,b=log(c(0.75,1)),
                                 N=400,
                                 miss=TRUE,
                                 nSim=1e3),
          g0dot8N60missingHN=list(ystart=3,w=1,
                                pi.x=pi.norm, logphi=c(0,log(0.2)),
                                hr=h21,b=c(log(c(0.75,1)),qlogis(0.8)),
                                N=100,
                                miss=TRUE,
                                nSim=1e3),
          g0dot8N200missingHN=list(ystart=3,w=1,
                                 pi.x=pi.norm, logphi=c(0,log(0.2)),
                                 hr=h21,b=c(log(c(0.75,1)),qlogis(0.8)),
                                 N=400,
                                 miss=TRUE,
                                 nSim=1e3)) 
#############################################################################
NV=seq(30,300,30)
sims=list(g1dotNormN30=list(ystart=3,w=1,
                             pi.x=pi.norm, logphi=c(0.5,log(0.2)),
                             hr=h2,b=log(c(0.75,1)),
                             N=NV[1],
                             nFix=TRUE,
                             nSim=200))
for(i in 2:length(NV)){
  sims[[i]]=sims[[1]]
  sims[[i]]$N=NV[i]
  names(sims)[i]=paste('g1dotNormN',NV[i],sep="_")
}
modFitL=list()


for(j in 1:length(sims))
{
  cs=sims[[j]]
  gridx=seq(cs$w/100,cs$w,length=100)
  gridy=seq(cs$ystart/100,cs$ystart,length=100)
  (pars=c(cs$b,cs$logphi))
  dput(cs,paste('simResSettings-',names(sims)[j],'.robj',sep=''))
  
  for(i in 1:cs$nSim){
    message(Sys.time(),' - ',names(sims)[j],' : starting simulation ',i,'\n')
    simRes=sim.n(n=cs$N,ymin=0.01,ystart=cs$ystart,w=cs$w,hr=cs$hr,
                 b=cs$b,pi.x=cs$pi.x,logphi=cs$logphi,fix.n=cs$nFix)
    x=simRes$locs$x; y=simRes$locs$y
    
    est.yx=NULL
    est.yx=fityx(y=y,x=x,b=cs$b,hr=cs$hr,ystart=cs$ystart,
                 pi.x=cs$pi.x,logphi=cs$logphi,w=cs$w,hessian=TRUE)
    if(i==1 & j==1) {modFitL[[1]]=est.yx} else {modFitL[[length(modFitL)+1]]=esy.yx}
    sim=coveragep(fit=est.yx,true.hr=cs$hr,true.b=cs$b,interval=0.95,
              true.pi.x=cs$pi.x,true.logphi=cs$logphi,verbose=TRUE)
    sim$phat.bias.percent=(sim$phat/sim$p-1)*100
    sim$n=length(x);sim$N=cs$N
    sim$Nhat=sim$n/sim$phat
    sim$Nhat.bias.percent=(sim$Nhat/sim$N-1)*100
    if(i==1 & j==1) {sim.out=sim} else {sim.out=rbind.data.frame(sim.out,sim)}
    write.csv(sim.out,paste('simRes-',names(sims)[1],'.csv',sep=''),row.names=FALSE)
    dput(est.yx,paste('allMods-',names(sims)[1],'.robj',sep=''))
  #boxplot(lapply(split(sim.out$phat.bias.percent,sim.out$n),summary))
  }
  
} #end jth loop.
  
####process sim data with fixed n
#sims=read.csv('~/Downloads/simRes-g0equals08.csv')
sims=read.csv('~/Downloads/simRes-g0equals08.csv')
dim(sims)
names(sims)
ts=table(sims$covered)
ts[1]/ts[2]
par(mfrow=c(1,1))
hist((sims$phat/sims$p-1)*100,xlab='effective strip width % bias',
     main=paste(nrow(sims),'simulations with g(0)=0.8'))
summary((sims$phat/sims$p-1)*100)

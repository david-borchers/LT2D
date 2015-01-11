#20140731: simulation scenarios
#1
#detectability = uniform
#perpendicular density gradient = h2 (density c. 02.; at perpendicular trucnation distance,w=1)
sims=list(g1dot0detectUNIFpixhr2p=list(ystart=3,w=1,
                            pi.x=hr2.to.p, logphi=c(0.5,log(0.2)),
                            hr=h.const,b=NULL,
                            N=300,
                            nFix=TRUE,
                            nSim=200))


cs=sims[[1]]
gridx=seq(cs$w/100,cs$w,length=100)
gridy=seq(cs$ystart/100,cs$ystart,length=100)
(pars=c(cs$b,cs$logphi))
simRes=sim.n(n=cs$N,ymin=0.01,ystart=cs$ystart,w=cs$w,hr=cs$hr,
               b=cs$b,pi.x=cs$pi.x,logphi=cs$logphi,fix.n=cs$nFix)
  x=simRes$locs$x; y=simRes$locs$y
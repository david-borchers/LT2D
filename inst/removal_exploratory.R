p=0.1
h= -log(1-p)
t=1:200
Ft=function(t,h) {1-exp(t*h)}
ft=function(t,h) {(1-Ft(t-1,h))*h}

plot(t,ft(t,h),type="l")


h=function(t,b) exp(b*t)
Ft=function(t,b) {
  Fi=rep(NA,length(t))
  for(i in 1:length(t)){
    tt=1:t[i]
    Fi[i]=1-exp(sum(h(tt,b)))
  }
  return(Fi)
}
ft=function(t,b) {
  fi=rep(NA,length(t))
  for(i in 1:length(t)){
    fi[i]=(1-Ft(t[i]-1,b))*h(t[i],b)
  }
  return(fi)
}

h=function(t,b) 0.5/(1+exp(-b*t))


b=0.01
plot(t,h(t,b),type="l")
plot(t,ft(t,b),type="l")



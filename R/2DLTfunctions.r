#'@title Detection hazard function \code{h1} prob(detect | available at x,y)
#'
#'@description  This hazard function has the form k(r,y)=a*r^(-b) from Hayes and Buckland (1983)
#' p36. Note: This function uses x for perp. dist., they use y.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta), and the function returns
#'theta[1]*(y^2+x^2)^(-theta[2]/2).
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h1(0.5,0.5,b=log(c(0.001,1)))
#'@seealso \code{\link{h2}}, \code{\link{ghy}}, \code{\link{ghy2}}
#'@export
h1=function(y,x,b)
{
  fName='h1'
  ## Test comment.
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  return(theta[1]*(y^2+x^2)^(-theta[2]/2))
}

#'@title Detection hazard function \code{ghy} prob(detect | available at x,y)
#'
#'@description  This hazard function is a generalization of the form k(r,y)=a*r^(-b) from 
#'Hayes and Buckland (1983) p36, the generalization being that a parameter to be estimated
#'is added to y. When this parameter is zero you get the form of Hayes and Buckland (1983) p36.
#'Note: This function uses x for perp. dist., they use y.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta); theta[1] is as per \link{\code{h1}} 
#'parameter theta[1]; theta[2] is as per \link{\code{h1}} theta[2]; theta[3] is parameter that
#'is added to y to shift forward distance origin and allow p(0)<1.
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h1(0.5,0.5,b=log(c(0.001,1)))
#'ghy(0.5,0.5,b=log(c(0.001,1,0.0)))
#'ghy(0.5,0.5,b=log(c(0.001,1,0.01)))
#'@seealso \code{\link{h2}}, \code{\link{ghy2}}
#'@export
ghy=function(y,x,b)
{
  fName='ghy'
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b)
  #theta=c(exp(b[1:2]),b[3])
  theta1=theta[1]^(1/theta[2])
  return(((x/theta1)^2+((y+theta[3])/theta1)^2)^(-theta[2]/2))
}

#'@title Detection hazard function \code{ghy2} prob(detect | available at x,y)
#'
#'@description  This hazard function is a generalization of the form k(r,y)=a*r^(-b) from 
#'Hayes and Buckland (1983) p36, the generalization being that (1) a parameter to be estimated
#'is added to y, and (2) x and y have separate scale parameters. It is a generalization of
#'\link{\code{ghy}} to allow x and y to have separate scale parameters. 
#'Note: This function uses x for perp. dist., they use y.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta); theta[1] is as per \link{\code{h1}} 
#'parameter theta[1]; theta[2] is as per \link{\code{h1}} theta[2]; theta[3] is parameter that
#'is added to y to shift forward distance origin and allow p(0)<1; theta[4] is the equivalent of 
#'theta[1], but specific to y, whereas theta[1] is specific to x in this function.
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h1(0.5,0.5,b=log(c(0.001,1)))
#'ghy(0.5,0.5,b=log(c(0.001,1,0.0)))
#'ghy(0.5,0.5,b=log(c(0.001,1,0.01)))
#'ghy2(0.5,0.5,b=log(c(0.001,1,0.01,0.001)))
#'ghy2(0.5,0.5,b=log(c(0.001,1,0.01,0.005)))
#'@seealso \code{\link{h2}}, \code{\link{ghy2}}
#'@export
ghy2=function(y,x,b)
{
  
  fName='ghy2'
  if(length(b)!=4) {
    cat(b,"\n")
    stop("b must be vector of length 4.")
  }
  theta=exp(b)
  thetax=theta[1]^(1/theta[2])
  thetay=theta[4]^(1/theta[2])
  return(((x/thetax)^2+((y+theta[3])/thetay)^2)^(-theta[2]/2))
}

#'@title Detection hazard function \code{h2} prob(detect | available at x,y)
#'
#'@description  This hazard function has the form k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2 from Hayes and Buckland (1983) p37.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta)
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h2(0.5,0.5,b=log(c(0.75,1)))
#'@export
#'@seealso \code{\link{h1}}
h2=function(y,x,b)
  #-------------------------------------------------------------------------------
# Detection hazard function prob(detect | available at x,y),
# Corresponding to Hayes and Buckland (1983) k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2
# on p37.
# Note: I use x for perp. dist., they use y.
# Inputs:
#  b: log(theta), where theta is vector of hazard rate parameters
#-------------------------------------------------------------------------------
{
  fName='h2'
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  return(theta[1]*y*(y^2+x^2)^(-(theta[2]+3)/2))
}

#' @title Hazard detection function of form \code{h2} with g(0)<1
#' 
#'@description  This hazard function has the form \eqn{k(r,y)=c*[a*\sqrt(r^2-y^2)/r^{(b+1)}]; b>2} modified 
#'from Hayes and Buckland (1983) p37, where c is imperfect detectability i.e. \eqn{g(0)<1}.
#'
#'@references Hayes, R. J., and S. T. Buckland. "Radial-distance models for the line-transect method." Biometrics (1983): 29-42.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[2]} is log(theta) and \code{b[3]} is qlogis(b[3]) detectability at the observer i.e. \eqn{g(0)=c}
#'@return probability of detection given that an animal is availabe at location x,y
#'#'@examples
#'h21(0.5,0.5,b=c(log(c(0.75,1)),qlogis(0.9)))
#' @export
h21=function(y,x,b)
  #-------------------------------------------------------------------------------
# Detection hazard function prob(detect | available at x,y),
# Corresponding to Hayes and Buckland (1983) k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2
# on p37.
# Note: I use x for perp. dist., they use y.
# Inputs:
#  b: log(theta), where theta is vector of hazard rate parameters
#-------------------------------------------------------------------------------
{
  fName='h21'
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b[1:2])
  dF=function(y,theta) theta[1]*y*(y^2+x^2)^(-(theta[2]+3)/2)
  g0=plogis(b[3])
  return(dF(y,theta)*g0)
}

#'@title Three-parameter inverse power hazard detection function
#' 
#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*[theta[2]/(sqrt{theta[2]^2+x^2+y^2})]^(theta[3]+1).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[1]} is plogis(theta[1])  \code{b[2]} is 
#'log(theta[2]) and \code{b[3]} is log(b[3]). 
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'b=c(-23.725809, -3.136638,2.122910)
#'ip1(0.5,0.5,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ip1(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ip1=function(y,x,b)
{
  fName='ip1'
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
#  theta=exp(b[2:3])
#  dF=function(y,x,theta) (theta[1]/(theta[1]^2+x^2+y^2))^(theta[2]+1)
#  g0=plogis(b[1])
  theta=exp(b)
  p=theta[1]*(1/sqrt(1+(x/theta[2])^2+(y/theta[2])^2))^(theta[3]+1)
#  p=theta[1]*(theta[2]/sqrt(theta[2]^2+x^2+y^2))^(theta[3]+1)
  return(p)
}


#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*(1/sqrt(1+(x)^2+(y)^2))^(theta[2]+1).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b 2-parameter vector, where \code{b[1]} is log(theta[1]) and  \code{b[2]} is 
#'log(theta[2]). 
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'b=c(5.2919208, 8.4701307)
#'ip0(0.05,0.05,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ip0(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ip0=function(y,x,b)
{
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  p=theta[1]*(1/sqrt(1+(x)^2+(y)^2))^(theta[2]+1)
  return(p)
}


#'@title Three-parameter exponential power hazard detection function 
#' 
#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*exp(-(x^theta[3]+y^theta[3])/(theta[2]^theta[3])).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[1]} is plogis(theta[1])  \code{b[2]} is 
#'log(theta[2]) and \code{b[3]} is log(b[3]). 
#'@return probability of detection given that an animal is availabe at location x,y
#'#'@examples
#'b=c(1, -4, 1)
#'ep1(0.5,0.5,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ep1(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ep1=function(y,x,b)
{
  fName='ep1'
  if(length(b)!=3) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b[2:3])
  dF=function(y,x,theta) exp(-(x^theta[2]+y^theta[2])/(theta[1]^theta[2]))
  g0=plogis(b[1])
  return(dF(y,x,theta)*g0)
}


#'@title Four-parameter inverse power hazard detection function
#' 
#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*(1/(1+(x/theta[2])^2+(y/theta[4])^2))^(theta[3]+1).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[1]} is plogis(theta[1]) \code{b[2]} is
#'log(theta[2]), where \code{theta[2]} is the scale parameter for x, \code{b[3]} is log(theta[3])
#'and \code{b[2]} is log(theta[4]), where \code{theta[4]} is the scale parameter for y. 
#'@return probability of detection given that an animal is availabe at location x,y
#'#'@examples
#'b=c(-23.725809,-3.136638,2.122910,-3.136638)
#'ip2(0.5,0.5,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ip2(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ip2=function(y,x,b)
{
  fName='ip2'
  if(length(b)!=4) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
#  theta=exp(b[2:4])
#  dF=function(y,x,theta) (1/(1+(x/theta[1])^2+(y/theta[3])^2))^(theta[2]+1)
#  g0=plogis(b[1])
#  return(dF(y,x,theta)*g0)
  theta=exp(b)
  p=theta[1]*(1/sqrt(1+(x/theta[2])^2+(y/theta[4])^2))^(theta[3]+1)
  return(p)
}



#'@title Four-parameter exponential power hazard detection function
#' 
#'@description  Inverse power hazard function, as per Borchers and Langrock (in press):
#'Has form h(y,x)=theta[1]*exp(-(x^theta[3]+y^theta[3])/(theta[2]^theta[3])).
#'
#'@references Borchers, D.L and Langrock, R."Double-observer line transect surveys with Markov-
#'modulated Poisson process models for animal availability" Biometrics (in press).
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, where \code{b[1]} is plogis(theta[1]) \code{b[2]} is
#'log(theta[2]), where \code{theta[2]} is the scale parameter for x, \code{b[3]} is log(theta[3])
#'and \code{b[2]} is log(theta[4]), where \code{theta[4]} is the scale parameter for y. 
#'@return probability of detection given that an animal is availabe at location x,y
#'#'@examples
#'b=c(1, -4, 1)
#'ep1(0.5,0.5,b=b)
#'yy=seq(0,0.03,length=100);xx=rep(0,100)
#'hh=ep1(yy,xx,b=b)
#'plot(yy,hh,type="l")
#' @export
ep2=function(y,x,b)
{
  fName='ep2'
  if(length(b)!=4) {
    cat(b,"\n")
    stop("b must be vector of length 3.")
  }
  theta=exp(b[2:4])
  dF=function(y,x,theta) exp(-((x/theta[1])^theta[2]+(y/theta[3])^theta[2]))
  g0=plogis(b[1])
  return(dF(y,x,theta)*g0)
}


#'@title Detection hazard function \code{h.exp2} prob(detect | available at x,y)
#'
#'@description  2-paramter Exponential power hazard model of Skaug & Schweder 1999. The gamma parameter is fixed at 2.

#'@references Skaug, Hans J., and Tore Schweder. "Hazard models for line transect surveys with independent observers." Biometrics 55.1 (1999): 29-36.
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, may be logged parameter values
#'@return probability of detection given that an animal is availabe at location x,y
#'@examples
#'h.exp2(0.5,0.5,b=log(c(0.75,0.9)))
#'@seealso \code{\link{h1}} \code{\link{h2}}
#'@export
h.exp2=function(y,x,b=c(0,0))
  #----------------------------------------------------------
# Detection hazard function prob(detect | available at x,y),
# 2-paramter Exponential power hazard model of Skaug & Schweder 1999.
# (gama fixed equal to 2).
#----------------------------------------------------------
{
  fName='h.exp2'
  mu=exp(b[1])
  sigma=exp(b[2])
  gama=2
  hr=mu*exp(-(x^gama + y^gama)/sigma^gama)
  return(hr)
}

#'@title Detection hazard function of Hiroshi et al. (2003)
#'#'
#'@description  2-parameter hazard model of Hiroshi, et al. 2003.
#'
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b parameter vector, may be logged parameter values
#'@references Okamura, Hiroshi et al. (2003) Abundance Estimation of Diving Animals by the Double-Platoform Line Transect Method, Biometrics 59(3):512-520.
#'@return probability of detection given that an animal is available at location x,y
#'@examples
#'h.okamura(0.5,0.5,b=log(c(0.5,0.9)))
#'@seealso \code{\link{h1}} \code{\link{h2}} \code{\link{h.exp2}}
#'@export
h.okamura=function(y,x,b=c(0,0))
  #----------------------------------------------------------
# Detection hazard function prob(detect | available at x,y).
# From Okamura's paper
#----------------------------------------------------------
{
  fName='h.okamura'
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  sigma.x=exp(b[1])
  sigma.y=exp(b[2])
  return(exp(-(x/sigma.x + y/sigma.y)))
}


#'@title Detection hazard function \code{h.const} prob(detect | available at x,y)
#'
#'@description  This function is for constant detectability throughout 
#'
#'@param y Forward distance
#'@param x perpendicular distance
#'@param b single value parameter vector (giving probability of detection).
#'@return probability of detection given that an animal is available at location x,y
#'@examples
#'h.const(0.5,0.5,b=1)
#'@export
#'@seealso \code{\link{h1}} \code{\link{h2}} \code{\link{h.exp2}} \code{\link{h.okamura}}
h.const=function(y,x,b=1) {fName='h.const'
  return(rep(b[1],length(y)))}

#'@title Half-normal form for perpendicular animal density function
#'
#'@description Half-normal distribution of perpendicular animal density.  
#'Truncation occurs at x=w (the perpendicular truncation distance)
#'
#'@param x prependicular trackline distance
#'@param logphi numeric vector; parameters, some of which may be logged
#'@param w perpendicular truncation distance
#'@return \eqn{\pi(x)} animal density at distance x
#'@examples
#'plot(seq(0,1,length=100),pi.hnorm(x=seq(0,1,length=100),logphi=0.5,w=1),
#'type='l',xlab='Perp. distance, x',ylab=expression(pi(x)))
#'@export
pi.hnorm=function(x,logphi,w){
  fName='pi.hnorm'
  hnF=function(x,logphi) exp(-x^2/(2*exp(logphi[1])^2))
  return(hnF(x,logphi)/integrate(hnF,0,w,logphi)$value)
}

#'@title Complementary half-normal form for perpendicular animal density function
#'
#'@description Complementary half-normal distribution of perpendicular animal density.  
#'Truncation occurs at x=w (the perpendicular truncation distance)
#'
#'@param x prependicular trackline distance
#'@param logphi numeric vector; parameters, some of which may be logged
#'@param w perpendicular truncation distance
#'@return \eqn{\pi(x)} animal density at distance x
#'@examples
#'plot(seq(0,1,length=100),pi.hnorm(x=seq(0,1,length=100),logphi=0.5,w=1),
#'type='l',xlab='Perp. distance, x',ylab=expression(pi(x)))
#'@export
pi.chnorm=function(x,logphi,w){
  fName='pi.chnorm'
  chnF=function(x,logphi) 1-exp(-(x-logphi[1])^2/(2*exp(logphi[2])^2))
  return(chnF(x,logphi)/integrate(chnF,0,w,logphi)$value)
}


#'@title Truncated normal form for perpendicular animal density function
#'
#'@description Truncated normal distribution of perpendicular animal density.  
#'Truncation occurs at x=0 (on the track line) and x=w (perpendicular truncation distance)
#'
#'@param x prependicular trackline distance
#'@param logphi numeric vector; parameters, some of which may be logged
#'@param w perpendicular truncation distance
#'@return \eqn{\pi(x)} animal density at distance x
#'@examples
#'plot(seq(0,1,length=100),pi.norm(x=seq(0,1,length=100),logphi=c(0.5,log(0.3)),w=1))
#'@export
pi.norm=function(x,logphi,w)
  #-------------------------------------------------------------------------------
# Animal density function (with respect to perp dist, x).
# Inputs:
#  logtphi: theta is vector of parameters, some of which may be logged
#           (see below).
#  w      : perp. truncation dist.
#-------------------------------------------------------------------------------
{
  fName='pi.norm'
  if(length(logphi)!=2) {
    cat(logphi,"\n")
    stop("logphi must be vector of length 2.")
  }
  if(any(x>w)) stop("x can't be greater than w")
  mu=logphi[1]
  sigma=exp(logphi[2])
  f=dnorm(x,mean=mu,sd=sigma)
  denom=(pnorm(w,mean=mu,sd=sigma)-pnorm(0,mean=mu,sd=sigma))
  if(denom>0) f=f/denom else f=0
  return(f)
}

#'@title Uniform perpendicular animal density function
#'
#'@description Uniform distribution for perpendicular animal density.  
#'
#'@param x prependicular trackline distance
#'@param logphi Not used
#'@param w perpendicular truncation distance
#'@return \eqn{\pi(x)} constant animal density of 1/w 
#'@examples
#'plot(seq(0,1,length=100),pi.const(x=seq(0,1,length=100),w=1))
#'@export
pi.const=function(x,logphi=NULL,w){
  fName='pi.const'
  return(rep(1/w,length(x)))}



#' @title Calculate HR perp dist function from hazard 
#' 
#' @description Implements the g(y) on page 36 or that on page 37 of Hayes and Buckland (1983) 
#' from the hazard function given on page 37 or 38, respectively, of that paper (and implemented 
#' in function \code{\link{h2}}).
#'   
#' @param x perpendicular distance.
#' @param b vector of parameters of \code{\link{h2}} (on log scale).
#' @param hr hazard rate function name (must be character): only "h1" and "h2" valid
#' @return value of hazard rate detection function at x
#' @examples
#' b=log(c(0.75,0.1))
#' x=seq(0,1,length=50)
#' plot(x,hr.to.p(x,b,hr='h2'),type="l",ylim=c(0,1),xlab="Perpendicular distance",ylab="P(detect)")
#'@export
hr.to.p=function(x,b,hr){
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  if(theta[2]<=1) {
    warning("exp(b[2])<=1 so setting it equal to 1+1e-10")
    theta[2]=1+1e-10 # gamma can't deal with zero, so get close to it
  }
  a1=(theta[1]*gamma((theta[2]-1)/2)*gamma(0.5))/(2*gamma(theta[2]/2))
  return(1-exp(-a1*x^(-(theta[2]-1))))
}


#' @title Perpendicular animal density function calulated from hazard rate \code{h1}
#' 
#' @description Calculates perpendicular animal density \eqn{\pi(x)} using the prependicular 
#' distance function of the hazard rate function \eqn{k(x,y)=a*r^{-b}} on page 36  of Hayes and Buckland (1983) 
#' from the hazard function given on page 37, of that paper (and implemented 
#' in function \code{\link{hr1.to.p}}).
#'   
#' @param x perpendicular distance.
#' @param logphi vector of parameters of \code{\link{h1}} (on log scale).
#' @param w perpendicular truncation distance
#' @return Animal density at x calculated from hazard rate \code{\link{h1}}
#' @examples
#' logphi=log(c(0.01,1.01))
#' x=seq(0,1,length=50)
#' plot(x,pi.hr1(x,logphi,w=1),type="l",xlab="Perpendicular distance",
#' ylab=expression(pi(x)))
#' @seealso \code{\link{h1}} \code{\link{hr1.to.p}}
#' @export
pi.hr1=function(x,logphi,w)
  #-------------------------------------------------------------------------------
# Animal density function (with respect to perp dist, x).
# Inputs:
#  logtphi: theta is vector of parameters, some of which may be logged
#           (see below).
#  w      : perp. truncation dist.
#-------------------------------------------------------------------------------
{
  fName='pi.hr1'
  if(length(logphi)!=2) {
    cat(logphi,"\n")
    stop("logphi must be vector of length 2.")
  }
  if(any(x>w)) stop("x can't be greater than w")
  f=hr1.to.p(x,b=logphi)/integrate(hr1.to.p,lower=0,upper=w,b=logphi)$value
  return(f)
}

#' @title Calculate harzard rate perpendicular distance function from hazard rate \code{h1} 
#' 
#' @description Implements the a prependicular distance function of the hazard rate function
#' k(x,y)=a*r^(-b) on page 36  of Hayes and Buckland (1983) 
#' from the hazard function given on page 37 of that paper (and implemented 
#' in function \code{\link{h1}}).
#'   
#' @param x perpendicular distance.
#' @param b vector of parameters of \code{\link{h1}} (on log scale).
#' @param w perpendicular truncation distance
#' @return value of hazard rate detection function at x
#' @examples
#' b=log(c(0.01,1.01))
#' x=seq(0,1,length=50)
#' plot(x,hr1.to.p(x,b),type="l",ylim=c(0,1),xlab="Perpendicular distance",ylab="P(detect)")
#' @seealso \code{\link{h1}} \code{\link{pi.hr1}}
#' @export
hr1.to.p=function(x,b,w){
  fName='hr1.to.p'
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  if(theta[2]<=1) {
    warning("exp(b[2])<=1 so setting it equal to 1+1e-10")
    theta[2]=1+1e-10 # gamma can't deal with zero, so get close to it
  }
  a1=(theta[1]*gamma((theta[2]-1)/2)*gamma(0.5))/(2*gamma(theta[2]/2))
  return(1-exp(-a1*x^(-(theta[2]-1))))
}


#' @title Perpendicular animal density function calulated from hazard rate \code{hr2}
#'  
#' @description Calculates perpendicular animal density \eqn{\pi(x)} using the prependicular 
#' distance function of the hazard rate function \eqn{k(r,y)=a \sqrt{(r^2-y^2)}/r^{(b+1)}; b>2}
#' on page 36  of Hayes and Buckland (1983) #' from the hazard function given on page 37, 
#' of that paper (and implemented in function \code{\link{hr1.to.p}}).
#'   
#' @param x perpendicular distance.
#' @param logphi vector of parameters of \code{\link{h1}} (on log scale).
#' @param w perpendicular truncation distance
#' @return Animal density at x calculated from hazard rate \code{\link{h1}}

#' @examples
#' x=seq(0,1,length=50)
#'plot(x,pi.hr2(x,logphi=c(-0.2876821, -2.3025851),w=1),
#'xlab='x',ylab=expression(pi(x)),type='l')
#' @seealso \code{\link{h1}} \code{\link{h2}} \code{\link{pi.hr1}}
#' @export
pi.hr2=function(x,logphi,w)
  #-------------------------------------------------------------------------------
# Animal density function (with respect to perp dist, x).
# Inputs:
#  logtphi: theta is vector of parameters, some of which may be logged
#           (see below).
#  w      : perp. truncation dist.
#-------------------------------------------------------------------------------
{
  fName='pi.hr2'
  if(length(logphi)!=2) {
    cat(logphi,"\n")
    stop("logphi must be vector of length 2.")
  }
  if(any(x>w)) stop("x can't be greater than w")
  f=hr2.to.p(x,b=logphi)/integrate(hr2.to.p,lower=0,upper=w,b=logphi)$value
  return(f)
}

#' @title Calculate harzard rate perpendicular distance function from hazard rate \code{h2} 
#' 
#' @description Implements the a prependicular distance function of the hazard rate function
#' \eqn{k(r,y)=a \sqrt(r^2-y^2)/r^{(b+1)}; b>2} on page 37  of Hayes and Buckland (1983) 
#' from the hazard function given on page 38 of that paper (and implemented 
#' in function \code{\link{h2}}).
#'   
#' @param x perpendicular distance.
#' @param b vector of parameters of \code{\link{h2}} (on log scale).
#' @param w perpendicular truncation distance
#' @return value of hazard rate detection function at x
#' @examples
#' x=seq(0,1,length=50)
#' plot(x,hr2.to.p(x,b=c(-0.2876821, -2.3025851),w=1),ylim=c(0,1),
#' xlab='Perp. distance, x', ylab='P(detect)',type='l')
#' @export
hr2.to.p=function(x,b,w){
  fName='hr2.to.p'
  if(length(b)!=2) {
    cat(b,"\n")
    stop("b must be vector of length 2.")
  }
  theta=exp(b)
  return(1-exp(-theta[1]/(theta[2]+1)*x^(-(theta[2]+1))))
}

#'@title Waiting distance pdf
#'
#'@description Calculates the pdf of the 'waiting distance' \eqn{f(y,x)=h(y,x)*\exp(-\int_y^{ystart} h(t,x) dt)}.
#'
#'@param y scalar or vector; forward distance
#'@param x scale or vector; perp. distance
#'@param b two-element vector of hazard rate parameters, some of whihc may be logged
#'@param hr hazard rate function 
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param nint number of intervals in numerical integration.
#'@return pdf of waiting distance at x,y
#'@details Need to ensure the hazard function has decayed to (very close to) zero by \code{ystart}.
#'@examples
#'w=1; ystart=4
#'gridx=seq(0,w,length=50); gridy=seq(0,ystart,length=50)
#'f=outer(gridy,gridx,FUN=fyx,b=c(-0.2876821, -2.3025851),hr=h2,ystart)
#'persp(gridx,gridy,t(f),theta=45,phi=35,zlab="f(y|x)")
#'@export
fyx=function(y,x,b,hr,ystart,nint=100)
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  n=length(x)
  f=intval=rep(NA,n)
  hr=match.fun(hr)
#  ylo=1e-10  # set to avoid evaluating hr at y=0, which gives Inf
  for(i in 1:n) {
#    y0=max(y[i],ylo)
#    dy=(ystart-y0)/nint/2                           # for crude integration
#    yy=seq(y0,ystart,length=(nint+1))[-(nint+1)]+dy # for crude integration
    dy=(ystart-y[i])/nint/2                           # for crude integration
    yy=seq(y[i],ystart,length=(nint+1))[-(nint+1)]+dy # for crude integration
    h=hr(yy,rep(x[i],nint),b)
    int=sum(hr(yy,rep(x[i],nint),b)*dy*2)  # crude integration
    intval[i]=exp(-int)
    #    int=integrate(f=hr,lower=max(y[i],ylo),upper=ystart,x=x[i],b=b)
    #    intval[i]=exp(-int$value)
  }
  hrval=hr(y,x,b)
  bads=which(hrval>=.Machine$double.xmax) # identify infinite hazards
  if(length(bads)>0) { # infinite hazard so p(detect)=0
    f[bads]=.Machine$double.xmax
    f[-bads]=hr(y[-bads],x[-bads],b)*intval[-bads]
  }else{
    f=hr(y,x,b)*intval
  }
return(f)
}

#'@title Numerical calculation of perpendicular detection function from a hazard
#'
#'@description Calculates the perpendicular detection function, \eqn{p(x)}, for a given hazard.
#'
#'@param x scale or vector; perp. distance
#'@param b two-element vector of hazard rate parameters, some of whihc may be logged
#'@param hr hazard rate function 
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param nint number of intervals in numerical integration.
#'@return probability of detection at x
#'@examples
#'gridx=seq(0,1,length=50)
#'p.x=px(gridx,b=c(-0.2876821, -2.3025851),
#'  hr=h2,ystart=4,nint=100)
#'plot(gridx,p.x,type="l",ylim=c(0,max(p.x)),
#' xlab="prep. distance, x",ylab="p(x)")
#'@export
#px=function(x,b,hr,ystart,nint=100)
#{
#  n=length(x)
#  p=int=rep(NA,n)
#  dy=(ystart-0)/nint/2                          # for crude integration
#  y=seq(0,ystart,length=(nint+1))[-(nint+1)]+dy # for crude integration
#  for(i in 1:n) {
#    p[i]=min(1,sum(fyx(y,rep(x[i],nint),b,hr,ystart)*dy*2)) # constrain to max of 1
#    #    int=integrate(f=fyx,lower=0,upper=ystart,x=x[i],b=b,hr=hr,ystart=ystart)
#    #    p[i]=int$value
#  }
#  return(p)
#}

px=function(x,b,hr,ystart,nint=100) 1-Sy(x,rep(0.0001,length(x)),ystart,b,hr)


#'@title Product of p(x) and pi(x)
#'
#'@description Returns product of perp dist det prob px() and animal dbn pi.x
#'
#'@param y scalar or vector; forward distance observations
#'@param x scale or vector; perp. distance observations
#'@param pars c(b,logphi); hazard rate and density log-parameters in a vector (see details).  
#'@param hr hazard rate function 
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param pi.x perpendicular distance density distribution
#'@param w perpendicular truncation distance.
#'@param length.b length of the hazard rate parameter vector
#'@return negative log likelihood for forward distance, \code{y} and perpendicular distance \code{x}.
#'@export
#'@details 
#'Must to ensure the hazard function has decayed to (very close to) zero by \code{ystart}.
#'The parameter vector, \code{pars}, must be passed in with two parameters for the hazard rate first, 
#'then two parameters for perpendicular density gradient \eqn{\pi(x)} i.e. \code{c(b,logphi)}. 
#'@examples
#'p.pi.x(x,b,hr,ystart,pi.x,logphi,w) 
p.pi.x=function(x,b,hr,ystart,pi.x,logphi,w) 
  return(px(x,b,hr,ystart)*pi.x(x,logphi,w))

#F.x=function(x,b,hr,ystart,pi.x,logphi,w) return((1-px(x,b,hr,ystart))*pi.x(x,logphi,w))

#'@title Negative log-likelihood for forward distance and perpendicular distance
#'
#'@description Calculates the negative log-likelihood for forward distance, \code{y}, and 
#'perpendicular distance, \code{x}, for a given hazard and perpendicular density distribution.
#'
#'@param y scalar or vector; forward distance observations
#'@param x scale or vector; perp. distance observations
#'@param pars c(b,logphi); hazard rate and density log-parameters in a vector (see details).  
#'@param hr hazard rate function 
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param pi.x perpendicular distance density distribution
#'@param w perpendicular truncation distance.
#'@param length.b length of the hazard rate parameter vector
#'@return negative log likelihood for forward distance, \code{y} and perpendicular distance \code{x}.
#'@details 
#'Must to ensure the hazard function has decayed to (very close to) zero by \code{ystart}.
#'The parameter vector, \code{pars}, must be passed in with two parameters for the hazard rate first, 
#'then two parameters for perpendicular density gradient \eqn{\pi(x)} i.e. \code{c(b,logphi)}. 
#'@examples
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y 
#'pars=c(b,logphi)
#'negloglik.yx(y,x,pars,hr,ystart,pi.x,w)
#'@seealso \code{\link{simXY}}
#'@export
negloglik.yx=function(pars,y,x,hr,ystart,pi.x,w,length.b=2,debug=FALSE)
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  if(debug) print(pars)
  n=length(y)
  # unpack parameters *** need to change if hr and pi.x don't have 2 pars each
  b=pars[1:length.b]
  if(as.character(substitute(pi.x))=="pi.const") logphi=NULL else logphi=pars[(1+length.b):length(pars)]
#  hr=match.fun(hr)
#  pi.x=match.fun(pi.x)
  llik=rep(NA,n)
  # caluclate numerator:
  num=sum(log(fyx(y,x,b,hr,ystart)) + log(pi.x(x,logphi,w)))
  # calculate denominator:
  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  #  F.x=function(x,b,hr,ystart,pi.x,logphi,w) return((1-px(x,b,hr,ystart))*pi.x(x,logphi,w))
  #  int=integrate(f=F.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  denom=log(int$value)
  # likelihood:
  llik=num-n*denom
  ##llik[is.nan(llik)]=-9e37
  if(debug) print(llik)
  #message(-llik)
  return(-llik)
}


#' Simulate sightings given a known perpendicular density distribution and hazard function
#' @description Simulates sightings from a known population given a density distribution
#' and hazard rate.  This function has been replaced by \code{\link{sim.n}}.
#' @param N animal population
#' @param pi.x function describing the perpendicular distance density distribution
#' @param logphi parameters for pi.x (some maybe logged)
#' @param hr function describing the hazard rate
#' @param b hazard rate parameter vector
#' @param w truncation distance
#' @param ystart max forward distance at which could possibly detect animal (see details).
#' @param xSampL length of x-dimension vector to sample perpendicular distances from.
#' @param discardNotSeen boolean; discard individuals not detected.  See details
#' @param ... arguments to be passed into \code{\link{simnhPP}}
#' @details if \code{discardNotSeen=FALSE} individuals that are not detected are 
#' assigned y-dimension distances = -999, otherwise \code{discardNotSeen=TRUE} 
#' invididuals are removed and not returned
#' @return list of \code{$locs} x and y coordinates for simulated sightings and 
#' \code{$settings} simulation settings.
#' @export
#'@examples
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y 
#'
simXY=function(N,pi.x,logphi,hr,b,w,ystart,xSampL=1e3,discardNotSeen=TRUE,...)
{
  xV=seq(0,w,length=xSampL) 
  x=sample(x=xV,size=N,replace=TRUE,
           prob=pi.x(x=xV,logphi=logphi,w=w))
  y=simnhPP(x,b,ystart,hr,...)
  if(discardNotSeen){
    keep=which(y>=0)
    n=length(keep)
    x=x[keep]
    y=y[keep]}
  return(list(locs=cbind.data.frame(x,y),settings=
                list(N=N,pi.x=pi.x,logphi=logphi,
                     hr=hr,b=b,w=w,ystart=ystart,
                     discardNotSeen=discardNotSeen)))
}

#'@title Plot the simulated positions
#' 
#' @description Plots the simulated data.  
#' 
#' @param simDat object from a call of \code{\link{simXY}}
#' @examples 
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=500 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#' plotSimNotUsed(simDat)
#' @export
plotSimNotUsed=function(simDat){
  x=simDat$locs$x; y=simDat$locs$y 
  pi.x=simDat$settings$pi.x
  par(mfrow=c(2,2),mar=c(3,3,3,3),mgp=c(1.5,0.5,0))
  gridx=seq(0,simDat$settings$w,length=1000)
  adbn=pi.x(gridx,logphi=simDat$settings$logphi,w=simDat$settings$w)
  plot(gridx,adbn,type="l",
       xlim=c(0,simDat$settings$w),
       xlab="perp. dist (x)",ylab="pi(x)",
       main='Perp. density function')
  rug(x,ticksize=0.1)
  plot(x,y,
       xlim=c(0,simDat$settings$w),
       ylim=c(0,simDat$settings$ystart),
       cex=0.6,
       main='Sighting locations',
       xlab='perp. distance (x)',
       ylab='forward dist. (y)')
  mtext(paste('N=',simDat$settings$N,'; n=',nrow(simDat$locs)))
  hist(x,freq=FALSE,xlab="perp. dist. (x)",
       main="perp. dist. (x)",
       xlim=c(0,simDat$settings$w))
  hist(y,freq=FALSE,xlab="forward dist. (y)",
       main='Forward dist. (y)',
       xlim=c(0,simDat$settings$ystart))
}

histline <-
  function(height,breaks,lineonly=FALSE,outline=FALSE,fill=FALSE,xlim=range(breaks),
           ylim=range(height),xlab="x",ylab="y",
           transpose=FALSE,...){
    #-------------------------------------------------------------------------------------
    # Takes bar heights (height) and cutbpoints (breaks), and constructs a line-only 
    # histogram from them using the function plot() (if lineonly==FALSE) or lines()
    # (if lineonly==TRUE). 
    # If fill==TRUE, uses polygon() to fill bars
    # If fill==TRUE, valid arguments to plot() or lines() are passed via argument(s) "..."
    # If outline==TRUE, only outline of histogram is plotted
    # If fill!=TRUE, valid arguments to polygon() are passed via argument(s) "..."
    #
    # DLB 2009
    #-------------------------------------------------------------------------------------
    
    n=length(height)
    if(length(breaks)!=(n+1)) stop("breaks must be 1 longer than height")
    if(outline) {
      y=c(0,rep(height,times=rep(2,n)),0)
      x=rep(breaks,times=rep(2,(n+1)))
    }   else {
      y=rep(0,4*n)
      x=rep(0,4*n+2)
      for(i in 1:n) {
        y[((i-1)*4+1):(i*4)]=c(0,rep(height[i],2),0)
        x[((i-1)*4+1):(i*4)]=c(rep(breaks[i],2),rep(breaks[i+1],2))
      }
      x=x[1:(4*n)]
    }
    if(transpose)
    {
      xstore=x
      x=y
      y=xstore
      xlimstore=xlim
    }
    if(lineonly) {
      if(!fill) lines(x,y,...)
      else polygon(x,y,...)
    } else {
      if(!fill) plot(x,y,type="l",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab,...)
      else {
        plot(x,y,type="n",xlim=xlim,ylim=ylim,xlab=xlab,ylab=ylab)
        polygon(x,y,...)
      }
    }
  }


#'@title Plot the simulated positions
#' 
#' @description Plots the simulated data.
#' 
#' @param simDat object from a call of \code{\link{simXY}}
#' @param nclass number of bins in x and y histograms
#' @param image Boolean \code{TRUE} image background on sightings scatter plot.
#' @param xlab x-axis label for sightings scatter plot
#' @param ylab y-axis label for sightings scatter plot
#' @param ... other arguments to be passed into \code{\link{plotfit.y}}
#' @seealso \code{\link{simXY}} \code{\link{plotSim}} \code{\link{plotfit.y}}
#'@export
#'@examples
#'n=100;ymin=0.01;ymax=5;W=1
#' b=log(c(0.75,1));logphi=c(0.5,log(0.2))
#' simDat=sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#'plotSim(simDat=simDat,nclass=12)
plotSim = function(simDat, nclass=10,xlab="", ylab="",image=FALSE,...){
  b=simDat$settings$b; hr=simDat$settings$hr; ystart=simDat$settings$ystart
  pi.x=simDat$settings$pi.x;logphi=simDat$settings$logphi; w=simDat$settings$w
  x=simDat$locs$x; y=simDat$locs$y 
  pi.x=simDat$settings$pi.x
  
  zones=matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
  layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5))
  xhist = hist(x, plot=FALSE,breaks=nclass)
  yhist = hist(y, plot=FALSE,breaks=nclass)
  top = max(c(xhist$counts, yhist$counts))
  par(mar=c(3,3,1,1))
  gridx=seq(0,w,length=50); gridy=seq(0,ystart,length=50)
  #f=outer(gridy,gridx,FUN=fyx,b=b,hr=hr,ystart=ystart) * pi.x(gridx,logphi=logphi,w=w)
  #gridx=seq(0,w,length=100)#c(0+dx,xhist$mids,w-dx)
  adbn=pi.x(gridx,logphi=logphi,w=w)
  
  est=list(b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  ly=plotfit.y(y,x,est=est,nclass=nclass,plot=FALSE,lineonly=FALSE,nint=50,...)
  yd=ly$fy.[floor(seq(1,length(ly$fy.),length.out=50))]
  if(image){
    image(x=gridx,y=gridy,
          t(matrix(rep(adbn,length(gridy)),nrow=length(gridy),byrow=T)*yd),
          xlab=xlab,ylab=ylab)
    points(x,y,cex=0.4,pch=19)}
  if(!image)
    plot(x,y,xlim=c(0,w),ylim=c(0,ystart),xlab=xlab,ylab=ylab,pch=19,cex=0.4)
  par(mar=c(0,3,1,1))
  histline(height=xhist$density,breaks=xhist$breaks,
           xlim=c(0,w),xaxt='n',ylab='Density')
  #dx=mean(diff(xhist$mids))/2
  #adbn=adbn/sum(adbn)
  lines(gridx,adbn,lwd=2)#*length(x),lwd=2)
  
  par(mar=c(3,0,1,1))
  
  histline(height=yhist$density,breaks=yhist$breaks,transpose=TRUE,xlab='Density',
           xlim=c(0,max(ly$fy.)),ylim=c(0,ystart),yaxt='n',ylab='')
  #deny<<-ly$fy./sum(ly$fy.)
  #lines(ly$fy.,ly$gridy,lwd=2)
  lines(yd,gridy,lwd=2)
  par(oma=c(3,3,0,0))
  mtext(xlab, side=1, line=1, outer=TRUE, adj=0, 
        at=.8 * (mean(x) - min(x))/(max(x)-min(x)))
  mtext(ylab, side=2, line=1, outer=TRUE, adj=0, 
        at=(.8 * (mean(y) - min(y))/(max(y) - min(y))))
}

negloglik.yx2=function(y,x,ps,hr,b,ys,pi.x,logphi,w)
  #-------------------------------------------------------------------------------
# Returns negative log likelihood for forward dist y and perp. dist. x.
# Inputs:
#  y       : forward distances (scalar or vector)
#  x       : perp. distances (scalar or vector)
#  pars    : c(b,logphi); haz rate and density log-parameters in a vector 
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance
#  pi.x    : name of animal density function to use.
#  w       : perp. truncation dist.
# ------
#  NOTE: code at *** must be changed if hr and pi.x don't have 2 pars each
# ------
#-------------------------------------------------------------------------------
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  ystart=ys
  hr=match.fun(hr)
  pi.x=match.fun(pi.x)
  n=length(y)
  # unpack parameters *** need to change if hr and pi.x don't have 2 pars each
  #b=pars[1:2]
  #logphi=pars[3:4]
  llik=rep(NA,n)
  # caluclate numerator:
  num=sum(log(fyx(y,x,b,hr,ystart)) + log(pi.x(x,logphi,w)))
  # calculate denominator:
  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  #  F.x=function(x,b,hr,ystart,pi.x,logphi,w) return((1-px(x,b,hr,ystart))*pi.x(x,logphi,w))
  #  int=integrate(f=F.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  denom=log(int$value)
  # likelihood:
  llik=num-n*denom
  #message(-llik)
  return(-llik)
}

#'@title Maximum likelihood estimation for unknown hazard and perpendicular distance distribution
#'
#'@description Uses \code{\link{optim}} to obtain a MLE for the hazard function and animal
#'perpendicular distance distribution.  Functional forms for the hazard and perpendicular distance 
#'distribution must be specified.
#'
#'@param y scalar or vector; forward distance observations
#'@param x scale or vector; perp. distance observations
#'@param b two-element vector of hazard rate parameters, some of whihc may be logged  
#'@param hr hazard rate function 
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param pi.x perpendicular distance density distribution
#'@param logphi parameters for pi.x (some maybe logged)
#'@param w perpendicular truncation distance.
#'@param control see \code{\link{optim}} control
#'@param hessian return hessian.  See also \code{\link{optim}}.
#'@param corrFlag=0.7 Absolute parameter correlation value above which a warning is issued.
#'@param ... arguments to be passed into \code{\link{optim}}
#'@return 
#'\code{\link{optim}} fit object and \cr
#'\code{$hr} = hazard rate function used.\cr
#'\code{$pi.x} = perpendicular distance function used.\cr
#'\code{$ystart} = ystart max forward distance detection used.\cr
#'\code{$w} = perpendicular truncation distance used.\cr
#'\code{$b} = estimated hazard parameters\cr
#'\code{$dat} = data frame with data (\code{$x} and \code{$y})\cr
#'\code{$logphi} \cr
#'\code{AIC} AIC value\cr
#'And if \code{hessian=TRUE}:\cr
#'\code{vcov} variance covariance matrix.  Will warn if there is a problem inverting 
#'the hessian.\cr
#'\code{CVpar} Coefficient of variation for each paramter estimate. \cr
#'\code{error} Boolean, \code{TRUE} if convergence!=0 or problem inverting the hessian, 
#'or parameter correlation is exceeded.\cr
#'@details Must to ensure the hazard function has decayed to (very close to) zero by 
#'\code{ystart}.
#'
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y 
#'fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'}
#'@seealso \code{\link{negloglik.yx}}
#'@export
fityx=function(y,x,b,hr,ystart,pi.x,logphi,w,control=list(),hessian=FALSE,corrFlag=0.7,debug=FALSE,...)
{
  if(is.function(hr)){
    eval(parse(text=fNameFinder(hr)))
    hrname=fName  
  }else{
    hrname=hr
    hr=match.fun(hr)
  }
  
  if(is.function(pi.x)){
    eval(parse(text=fNameFinder(pi.x)))
    piname=fName  
  }else{
    piname=pi.x
    pi.x=match.fun(pi.x)
  }
  
  if(piname=="pi.const") pars=b else pars=c(b,logphi)
  length.b=length(b)
  fit=optim(par=pars,fn=negloglik.yx,y=y,x=x,hr=hr,ystart=ystart,pi.x=pi.x,w=w,
            length.b=length.b,
            hessian=hessian,debug=debug,control=control,...)
  fit$error=FALSE
  if(fit$convergence!=0){
    warning('Convergence issue (code = ', fit$convergence,') . Check optim() help.')
    fit$error=TRUE
  }
  fit$hr=hrname
  fit$pi.x=piname
  fit$ystart=ystart
  fit$w=w
  fit$b=fit$par[1:length.b]  # ***
  if(length.b!=length(pars)){
    fit$logphi=fit$par[(1+length.b):length(pars)]
  }else{
      fit$logphi=NA
  }    # ***
  fit$AIC=2*fit$value+2*length(fit$par)
  fit$dat=data.frame(x=x,y=y) # attach data to fitted object
  if(hessian){
    mNames=paste('b',1:length.b,sep='')
    if(!all(is.na(fit$logphi))) mNames=c(mNames,paste('logphi',1:length(fit$logphi),sep=''))
    fit$vcov=solve(fit$hessian)
    if(any(diag(fit$vcov)<=0)){
      warning('Failed to invert hessian.  Model covergance problem in fityx?')
      fit$error=TRUE
      fit$CVpar=rep(NA,length(fit$par))} else {
        fit$CVpar=sqrt(diag(solve(fit$hessian)))/abs(fit$par)}
    fit$corr=cov2cor(fit$vcov)
    row.names(fit$corr)=mNames
    colnames(fit$corr)=mNames
    corr=fit$corr
    corr[upper.tri(corr,diag=TRUE)]=NA
    corrIND=which(abs(corr)>corrFlag,arr.ind=T)
    if(nrow(corrIND)){
      warning('absolute correlation exceeds ',corrFlag,' in parameter estimates: ',
              paste(paste(mNames[corrIND[,1]],mNames[corrIND[,2]],sep=' to '),collapse='; '))
      fit$error=TRUE}
  }
  return(fit)
}

negloglik.yx.w=function(y,x,pars,hr,ystart,pi.x,logphi,w)
  #-------------------------------------------------------------------------------
# Returns negative log likelihood for forward dist y and perp. dist. x, taking
# pi(x) as known.
# Corresponding to Hayes and Buckland (1983) k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2
# on p37.
# Note: I use x for perp. dist., they use y.
# Inputs:
#  y       : forward distances (scalar or vector)
#  x       : perp. distances (scalar or vector)
#  pars    : b; haz rate log-parameters in a vector 
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance
#  pi.x    : name of animal density function to use.
#  logphi  : vector of animal density function parameters (some logged)
#  w       : perp. truncation dist.
# ------
#  NOTE: code at *** must be changed if hr and pi.x don't have 2 pars each
# ------
#-------------------------------------------------------------------------------
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  
  hr=match.fun(hr)
  pi.x=match.fun(pi.x)
  n=length(y)
  llik=rep(NA,n)
  # caluclate numerator:
  num=sum(log(fyx(y,x,b=pars,hr,ystart)) + log(pi.x(x,logphi,w)))
  # calculate denominator:
  int=integrate(f=p.pi.x,lower=0,upper=w,b=pars,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  #  F.x=function(x,b,hr,ystart,pi.x,logphi,w) return((1-px(x,b,hr,ystart))*pi.x(x,logphi,w))
  #  int=integrate(f=F.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  denom=log(int$value)
  # likelihood:
  llik=num-n*denom
  
  return(-llik)
}



fityx.w=function(y,x,b,hr,ystart,pi.x,logphi,w,control)
{
  pars=b
  fit=optim(par=pars,fn=negloglik.yx.w,y=y,x=x,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w,hessian=FALSE,control=control)
  fit$hr=hr
  fit$pi.x=pi.x
  fit$ystart=ystart
  fit$w=w
  fit$b=fit$par
  fit$logphi=logphi
  return(fit)
}




negloglik.x=function(x,pars,hr,ystart,pi.x,logphi,w,nint=100)
  #-------------------------------------------------------------------------------
# Returns negative log likelihood for perp. dist. x. given dbn pi.x and logphi
# Corresponding to Hayes and Buckland (1983) k(r,y)=a*sqrt(r^2-y^2)/r^(b+1); b>2
# on p37.
# Note: I use x for perp. dist., they use y.
# Inputs:
#  x       : perp. distances (scalar or vector)
#  pars    : b; haz rate log-parameters in a vector 
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance
#  pi.x    : name of animal density function to use.
#  logphi  : vector of animal density function parameters (some logged)
#  w       : perp. truncation dist.
#-------------------------------------------------------------------------------
{
  hr=match.fun(hr)
  pi.x=match.fun(pi.x)
  n=length(x)
  # unpack parameters *** need to change if hr and pi.x don't have 2 pars each
  b=pars
  llik=rep(NA,n)
  # caluclate numerator: 
  num=sum(log(px(x,b,hr,ystart,nint)) + log(pi.x(x,logphi,w)))
  # calculate denominator:
  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  #  F.x=function(x,b,hr,ystart,pi.x,logphi,w) return((1-px(x,b,hr,ystart))*pi.x(x,logphi,w))
  #  int=integrate(f=F.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)
  denom=log(int$value)
  # likelihood:
  llik=num-n*denom
  
  return(-llik)
}


fitx=function(x,b,hr,ystart,pi.x,logphi,w,control)
{
  pars=b
  fit=optim(par=pars,fn=negloglik.x,x=x,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w,hessian=FALSE,control=list(trace=5))
  fit$hr=hr
  fit$pi.x=pi.x
  fit$ystart=ystart
  fit$w=w
  fit$b=fit$par           # ***
  fit$logphi=logphi
  return(fit)
}




Eyx=function(y,x,b,hr,ystart)
  #Eyx=function(y,x,b,hr,ystart,nint=500)
  #-------------------------------------------------------------------------------
# Returns Expectation \int_y^ystart h(t,x) dt.
# Inputs:
#  y       : forward dist. (scalar or vector)
#  x       : perp. dist. (scalar or vector)
#  b: log(theta), where theta is vector of hazard rate parameters
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance
#-------------------------------------------------------------------------------
{
  if(length(y)!=length(x)) stop("Lengths of x and y must be the same.")
  n=length(x)
  int=rep(NA,n)
  hr=match.fun(hr)
  ylo=1e-5  # set to avoid evaluating hr at y=0, which gives Inf
  for(i in 1:n) {
    y0=max(y[i],ylo)
    #    dy=(ystart-y0)/nint/2                           # for crude integration
    #    yy=seq(y0,ystart,length=(nint+1))[-(nint+1)]+dy # for crude integration
    #    int[i]=sum(hr(yy,rep(x[i],nint),b)*dy*2)  # crude integration
    int[i]=integrate(f=hr,lower=max(y[i],ylo),upper=ystart,x=x[i],b=b)$value
  }
  return(int)
}



#' @title Simulate forward distances
#' 
#' @description Simulate non-homogeneous Poisson Process data by solving inverse CDF.  
#' This function has been replaced by \code{\link{sim.n}}.
#' 
#'@param x scale or vector; perpendicular distance observations
#'@param b two-element vector of hazard rate parameters, some of whihc may be logged  
#'@param ystart max forward distance at which could possibly detect animal (see details).
#'@param hr hazard rate function 
#'@param miss If TRUE, allows animals not detected by y=0 to have no y.
#            Misses are indicated by y==-999
#'@param ylo minimum forward distance.
#'@return vector of y distances if \code{miss=TRUE} missed animals represented by -999, otherwise
#'no animals are missed.
#'@details \code{miss} allows a fixed number of animal to be simulated i.e. known $N$ in strip, 
#'or $n$ where animals may have been missed.
#' @examples
#'ystart=4
#'hr=h2; b=log(c(0.75,1))
#'X=runif(50,0,1)
#'Y=simnhPP(x=X,b=b,ystart=ystart,
#'hr=h2,miss=TRUE)
#'length(Y[Y!=-999])
#'Y=simnhPP(x=X,b=b,ystart=ystart,
#'hr=h2,miss=FALSE)
#'length(Y[Y!=-999])
#' @seealso \code{\link{simXY}}
simnhPP=function(x,b,ystart,hr,miss=TRUE,ylo=1e-5)
  #-------------------------------------------------------------------------------
# Simulate non-homogeneous Poisson Process data by solving inverse CDF.
# Inputs:
#  x       : perp. dist. (scalar or vector)
#  b: log(theta), where theta is vector of hazard rate parameters
#  hr      : name of hazard rate function to use.
#  ystart  : max forward distance at which could possibly detect animal.
#            NB: need to ensure hazard has decayed to (very close to) zero by
#            this distance.
#  miss    : If TRUE, allows animals not detected by y=0 to have no y.
#            Misses are indicated by y==-999.
#-------------------------------------------------------------------------------
{
  obj=function(y,x,b,hr,ystart,u) return(((1-exp(-Eyx(y,x,b,hr,ystart)))-u)^2)
  n=length(x)
  u=runif(n)
  y=rep(-999,n)
  for(i in 1:n) {
    if(!miss) {
      while(y[i]<0) {
        if(u[i]>(1-exp(-Eyx(ylo,x[i],b,hr,ystart)))) u[i]=runif(1) # y<0, so missed: try again
        else {
          ymin=optimize(f=obj,interval=c(ylo,ystart),x[i],b,hr,ystart,u[i])
          y[i]=ymin$minimum
        }
      }
    } else if(u[i]<=(1-exp(-Eyx(ylo,x[i],b,hr,ystart)))) {
      ymin=optimize(f=obj,interval=c(ylo,ystart),x[i],b,hr,ystart,u[i])
      y[i]=ymin$minimum
    }
  }
  return(y)
}



#'@title Plot fitted hazard and perpendicular denisty distribution
#'
#'@description Plot fitted hazard and perpendicular denisty distribution 
#'resulting from a call of \code{\link{fityx}}.  When simulating, true hazard functions
#'and perpendicular density distribution can also be added to the plot.

#'@param x perpendicular distance observations
#'@param fit return from a call of \code{\link{fityx}}
#'@param nclass number of histogram classes
#'@param nint number of intervals in numerical integration
#'@param plot boolean, plot results
#'@param addTruth boolean add true hazard and perp. density when simulating
#'@param true.pi.x true perpendicular density distribution function used when simulating
#'@param true.logphi true perpendicular density distribution function 
#'parameters used when simulating
#'@param true.hr true hazard rate function used when simulating.
#'@param true.b true hazard rate function parameters used when simulating.
#'@param N true number of animals in simulated distribution, used to calculate bias (see details).
#'@param true.legend If true (and addTruth) plots legend for true functions in bottom left.
#'
#'@details When \code{N} is specified, bias in estimated number of animals ,\eqn{\hat N}, 
#' is calculated.
#'
#'@return
#'list with:
#'\code{$gridx} = x values used in plotting
#'\code{$p.xpifit} = product of perpendicular distance det probability p(x) 
#'and perpendicular animal  distribution \eqn{\pi(x)}.
#'\code{$mufit} = effective strip width \eqn{\hat p}
#'\code{$f.xfit} =
#'\code{$p.xfit} =
#'\code{$ptot} =
#'\code{$p.xfit.std} =
#'\code{$adbn} =
#'\code{$N} = true number of animals in population
#'\code{$n} = number of seen animals
#'\code{$Nhat} = estimated number of animals.
#'\code{$bias} = \eqn{\hat N} bias.

#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y 
#'est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'plotdat.yx=plotfit.x(x,est.yx,addTruth=TRUE,true.logphi=logphi,true.b=b,N=N)
#'}
#'@seealso \code{\link{fityx}}
#'@export
plotfit.x=function(x,est,nclass=10,nint=100,
                   plot=TRUE,dotitle="FALSE",
                   addTruth=FALSE,
                   true.pi.x=NULL,
                   true.logphi=NULL,
                   true.hr=NULL,
                   true.b=NULL,
                   N=NULL,...)
{
  Nhat.yx=bias=NULL
  b=est$b; hr=match.fun(est$hr); ystart=est$ystart; pi.x=match.fun(est$pi.x)
  logphi=est$logphi; w=est$w
  x=x[x>=0 & x<=w]
  f.x=p.x.std=adbnTRUE=0
  # calculate stuff to plot:
  gridx=seq(1e-10,w,length=100)
  # first do f(x)
  p.xpifit=p.pi.x(gridx,b,hr,ystart,pi.x,logphi,w)
  mufit=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr
                  ,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)$value 
  f.xfit=p.xpifit/mufit
  p.xfit=px(gridx,b,hr,ystart,nint=nint)
  ptot=integrate(f=px,lower=0,upper=w,b=b,hr=hr,ystart=ystart)$value
  p.xfit.std=p.xfit/ptot
  adbn=pi.x(gridx,logphi,w)
  
  if(addTruth) {
    if(!is.null(true.pi.x)) pi.x=true.pi.x
    if(!is.null(true.logphi)) logphi=true.logphi
    if(!is.null(true.hr)) hr=true.hr
    if(!is.null(true.b)) b=true.b
    p.xpi=p.pi.x(gridx,b,hr,ystart,pi.x,logphi,w)
    mu=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)$value 
    f.x=p.xpi/mu
    adbnTRUE=pi.x(gridx,logphi,w)
  }
  if(plot){
    breaks=seq(0,w,length=(nclass+1))
    hx=hist(x,breaks=breaks,plot=FALSE) # get hist bar heights
    ymax=max(f.xfit,p.xfit.std,adbn,f.x,p.x.std,adbnTRUE,hx$density)
    main=""
    if(dotitle) main="Fitted curves"
    if(addTruth) main="Fitted curves (grey=true)"
    hx=hist(x,breaks=breaks,freq=FALSE,ylim=c(0,ymax),
            main=main,xlab="perpendicular distance (x)",ylab="pdf")
    lines(gridx,f.xfit,lwd=1)
    # overlay p(x), scaled to have area=1
    lines(gridx,p.xfit.std,lty=2,col="black",lwd=2)
    # overlay animal pdf:
    lines(gridx,adbn,lty=3,col="black",lwd=2)
    
    if(addTruth) legend("topright",title="Estimated",legend=c("f(x)","p(x)","pi(x)"),
           col=c("black","black","black"),lwd=c(2,2,2),lty=c(1,2,3))
    else legend("topright",legend=c("f(x)","p(x)","pi(x)"),
                col=c("black","black","black"),lwd=c(2,2,2),lty=c(1,2,3))
    if(addTruth){
      lines(gridx,f.x,col="grey",lwd=2)
      p.x=px(gridx,b,hr,ystart,nint=nint)
      ptot=integrate(f=px,lower=0,upper=w,b=b,hr=hr,ystart=ystart)$value
      #      p.x.std=p.xfit.std=p.xfit/ptot
      #      lines(gridx,p.x*p.x.std,col="grey",lty=2,lwd=2)
      p.x.std=p.x/ptot
      lines(gridx,p.x.std,col="grey",lty=2,lwd=2)
      
      lines(gridx,adbnTRUE,col="grey",lty=3,lwd=2)
    }
  }
  if(!is.null(N)){
    n=length(x)
    Nhat.yx=n/mufit
    bias=(Nhat.yx/N-1)*100
    msg=paste("N=",N,"; n=",n,"; Nhat.yx=",signif(Nhat.yx,3),";``bias''=",signif(bias,3),"%\n",sep="")
    message(msg)
    if(plot) {
      mtext(msg,cex=0.8)
    }
  }else{N=N;n=NULL;Nhat=NULL;bias=NULL}
  invisible(list(gridx=gridx,p.xpifit=p.xpifit,mufit=mufit,
                 f.xfit=f.xfit,p.xfit=p.xfit,ptot=ptot,p.xfit.std=p.xfit.std,adbn=adbn,
                 N=N,n=n,Nhat=Nhat.yx,bias=bias))
}


#'@title Plot fitted f(y) and forward distance distribution
#'
#'@description Plot f(y) and forward distance distribution 
#'resulting from a call of \code{\link{fityx}}.
#'
#'@param y forward distance observations (if NULL, uses est$dat$y)
#'@param x perpendicular distance observations (if NULL, uses est$dat$x)
#'@param est return from a call of \code{\link{fityx}}
#'@param nclass number of histogram classes
#'@param breaks break points passed to hist (overrides nclass if not NULL)
#'@param plot boolean, plot results
#'@param lineonly if TRUE plots only f(y), else plots histogram of forward distances too
#'@param nint number of intervals to use in calculating f(y)
#'@param max.obs If TRUE, plots only up to maximum observed forward distance, else plots
#'up to est$ystart (the forward distance beyond which detection is assumed impossible).
#'@param add if TRUE, adds line to existing plot, else creates new plot. Only applicable 
#'if lineonly==TRUE.
#'@param ... other parameters passed to \code{hist} and \code{plot} (need to separate these two!)
#'
#'@details Plot f(y) and forward distance distribution resulting from a call of 
#'\code{\link{fityx}}, optionall with overlaid histogram of foward distances. 
#'Invisibly returns various f(y)'s, as detailed below.
#'
#'@return
#'Invisibly returns a list with these elements
#'\code{$gridx} = x values used in plotting
#'\code{$fy.x} = Unscaled f(y|x) for all the xs observed
#'\code{$fy.} = Unscaled mean of f(y|x) for all the xs observed
#'\code{$scaled. fy.} = Mean of f(y|x), scaled to integrate to 1
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y 
#'est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'plotfit.y(y,x,est.yx,nclass=10)
#'}
#'@seealso \code{\link{fityx}}
#'@export
plotfit.y=function(y=NULL,x=NULL,est,nclass=10,breaks=NULL,plot=TRUE,dotitle=FALSE,
                   lineonly=FALSE,nint=100,max.obs=TRUE,add=FALSE,...)
{
  b=est$b; hr=match.fun(est$hr); ystart=est$ystart; pi.x=match.fun(est$pi.x)
  logphi=est$logphi; w=est$w
  if(is.null(y)) y=est$dat$y
  if(is.null(x)) x=est$dat$x
  # calculate stuff to plot:
  n=length(y)
  res=100
  gridy=seq(1e-10,ystart,length=res)
  fy.x=matrix(rep(NA,n*res),nrow=n)
  for(i in 1:n) {
    fy.x[i,]=fyx(gridy,rep(x[i],res),b,hr,ystart,nint=nint)
  }
  
  fy.=apply(fy.x,2,mean)
  fy.area=sum((fy.[-1]+fy.[-length(fy.)])/2*diff(gridy))
  scaled.fy.=fy./fy.area
  if(plot){
    ymax=ystart
    if(max.obs) ymax=max(y)
    if(is.null(breaks)) breaks=seq(min(y,1e-10),ymax,length=(nclass+1)) 
    fy.area=sum((fy.[-1]+fy.[-length(fy.)])/2*diff(gridy))
    scaled.fy.=fy./fy.area
    if(lineonly) {
      if(add) lines(gridy,scaled.fy.,...)
      else plot(gridy,scaled.fy.,ylim=c(0,max(scaled.fy.)),type="l",
                xlab="forward distance (y)",ylab="f(y)",...)
    }
    else {
      # hst=hist(y,plot=FALSE)
      # hist(y,freq=FALSE,xlab="forward distance (y)",nclass=nclass,ylim=c(0,max(hst$intensities,fy.)))
      hst=hist(y,breaks=breaks,plot=FALSE,...)
      # cat("hist area=",hst$desity*diff(hst$breaks),"\n")
      hmax=max(scaled.fy.,hst$density)
      if(dotitle) hist(y,freq=FALSE,xlab="forward distance (y)",breaks=breaks,ylim=c(0,hmax),...)
      else hist(y,freq=FALSE,xlab="forward distance (y)",breaks=breaks,ylim=c(0,hmax),main="",...)
      lines(gridy,scaled.fy.,...)
      # cat("fy area=",sum((fy.[-1]+fy.[-length(fy.)])/2*diff(gridy)),"\n")
    }}
  
  invisible(list(gridy=gridy,fy.x=fy.x,fy.=fy.,scaled.fy.=scaled.fy.))
}



#'@title Calculate coverage probabilities of \eqn{\hat p} for simulated data
#'
#'@description Calculate coverage probabilities of \eqn{\hat p} using the delta method and assuming a 
#'log-normal error distribution.
#'
#'@param fit object resulting from a call of \code{\link{fityx}} (see details)
#'@param interval the interval used to determine coverage probability
#'@param true.hr true form of hazard rate function 
#'@param true.b true values of hazard rate parameters
#'@param true.pi.x true form of perpendicular density distribution pi(x)
#'@param true.logphi true values for pi(x) perpendicular density distribution parameters
#'@param verbose boolean.  FALSE only covered returned.  TRUE a data frame of covered and associated calculations returned (see returns).
#'@param type \code{LOGNORM} log-normal confidence intervals; \code{NORM} normal confidence intervals.
#'@details In the call of \code{\link{fityx}} must have \code{hessian=TRUE}
#'@return boolean; \code{TRUE} if within log-normal confidence interval, otherwise \code{FALSE}.
#'\code{verbose=TRUE}   a data frame of p phat, and CV[phat], interval, lower bound returned.
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y 
#'est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'coveragep(fit=est.yx,true.hr=hr,true.b=b,interval=0.95,
#'  true.pi.x=pi.x,true.logphi=logphi,verbose=TRUE)
#'}     
#'@seealso \code{\link{phat}} \code{\link{fityx}}
#'@export
coveragep=function(fit,true.hr,true.b,true.pi.x,true.logphi,type='LOGNORM',
                   interval=0.95,verbose=FALSE){
  lnci.nmin=function(stat,cv,stat.min=0,interval=interval){
    q=Mod(qnorm((1-interval)/2,0,1))
    varNhat=(stat*cv)^2
    cfactor=exp(q*sqrt(log(1+varNhat/(stat-stat.min)^2)))
    lower=stat.min+(stat-stat.min)/cfactor
    upper=stat.min+(stat-stat.min)*cfactor
    return(list(lower=lower,upper=upper))
  }
  
  if(!'hessian' %in% names(fit))
    stop('fit ARG must include a hessian matrix')
  pars=fit$par; 
  hr=match.fun(fit$hr); b=fit$b;
  ystart=fit$ystart; w=fit$w
  pi.x=fit$pi.x; logphi=fit$logphi
  #true p
  p=phat(w=w,hr=true.hr,b=true.b,ystart=ystart,pi.x=true.pi.x,logphi=true.logphi)
  #estimated p
  p.hat=phat(w=w,hr=hr,b=b,ystart=ystart,pi.x=pi.x,logphi=logphi)
  #calc. variance-covariance
  vcov=solve(fit$hessian)
  #Implement the delta method:
  #numerical differentiation
  dbyd=numericDeriv(quote(phat(w=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi)), c("b","logphi"))  
  dbyd=as.vector(slot(dbyd,'gradient'))
  var.p.hat=as.vector(t(dbyd)%*%vcov%*%dbyd) #$var[hat{p}(\hat{\Beta})$
  if(type=='LOGNORM')
    bounds=lnci.nmin(stat=p.hat,cv=sqrt(var.p.hat)/p.hat,interval=interval)
  if(type=='NORM'){
    bounds=list(lower=qnorm((1-interval)/2,p.hat,sqrt(var.p.hat)),
                upper=qnorm(interval+(1-interval)/2,p.hat,sqrt(var.p.hat)))
  }
  if(any(sapply(bounds,is.nan))){
    warning('One or both p.hat bounds NaN')
    return(data.frame(covered=NA,p=p,phat=p.hat,CV.phat=NA,
                      interval=interval,lower.bound=bounds$lower,upper.bound=bounds$upper))
  }
  covered=TRUE
  if(p < bounds$lower | p > bounds$upper) covered=FALSE
  if(!verbose)
    return(covered)
  return(data.frame(covered=covered,p=p,phat=p.hat,CV.phat=sqrt(var.p.hat)/p.hat,
                    interval=interval,lower.bound=bounds$lower,upper.bound=bounds$upper))
}

#'@title Calculate effective strip width
#'
#'@description Calculate effective strip width, \eqn{\hat p} for a given hazard rate function and 
#'perpendicular density distribution.
#'
#' @param w truncation distance
#' @param hr function describing the hazard rate
#' @param b hazard rate parameter vector
#' @param ystart maximum forward distance
#' @param pi.x function describing the perpendicular distance density distribution
#' @param logphi parameters for pi.x (some maybe logged)
#' @param fit=NULL alternatively just pass in an object resulting from a call of \link{fityx}.
#'@return Effective strip widith \eqn{\hat p}
#'@examples phat(w=1,hr=h2,b=log(c(0.75,1)),ystart=4,pi.x=pi.norm,logphi=c(0.5,log(0.2)))
#'@export
phat=function(fit=NULL,w=NULL,hr=NULL,b=NULL,ystart=NULL,pi.x=NULL,logphi=NULL)
{
  if(!is.null(fit)){
    #f=fit$p.pi.x; 
    upper=fit$w; b=fit$b; hr=match.fun(fit$hr)
    ystart=fit$ystart;pi.x=match.fun(fit$pi.x);logphi=fit$logphi;w=fit$w
  }
  int=integrate(f=p.pi.x,lower=0,upper=w,b=b,hr=hr,
                ystart=ystart,pi.x=pi.x,logphi=logphi,w=w)$value
  return(int)
}

#--------------- Functions added by DLB 24/7/14 ------------------------------

#' @title Simulate n sightings from NHPP
#' 
#' @description Simulates n sighting locations (x,y) given a perp.dist 
#' distribution model and detection hazard model
#' 
#' @param n sample size
#' @param ymin smallest forward distance
#' @param ymax largest forward distance
#' @param W perpendicular truncation distance
#' @param hfun detection hazard function
#' @param b vector of detection hazard function parameters
#' @param pi.x perpendicular distance distribution function
#' @param logphi vector with log of pi.x parameters
#' @param fix.n if TRUE sample size of exactly n is generated, else sample is
#'        generated from model with expected sample size n
#' @param intscale amount by which to multiply detection location pdf f(x,y)
#'        in order to get required sample size. Either an object of class 
#'        "ppscale" output by \code{\link{calc.lpars}} or NULL (in which case
#'        \code{\link{calc.lpars}} is called inside \code{sim.n}.
#' @param nbuffer amount by which to multiply the expected n by (given all model
#'        parameters and intscale) to reduce probability that generated n is less 
#'        that n on first call to NHPP generating funciton rpoispp. If NULL, 
#'        it is set to 1.25 inside \code{sim.n}.
#' @details Uses the \code{spatstat} function \code{rpoispp} to 
#' generate detections from a NHPP, with intensity parameter such that the 
#' expected (if fix.n is FALSE) or actual (if fix.n is TRUE) sample size is n.
#' 
#' @return a list object with element 1 a data.frame of simulated \code{x} and \code{y} sightings locations; element 2 \code{spatstat} object of class "ppp" with x- and y-coordinates 
#' of detections and element 2 simulation settings, comprising of \code{n}
#'@examples
#' \dontrun{
#' # simulate with fixed n:
#' n=100;ymin=0.01;ymax=5;W=2
#' b=log(c(0.75,1));logphi=c(0.5,log(0.3))
#' dat=sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#' dat$locs$n
#' 
#' plot(density(dat))
#' contour(density(dat),add=TRUE)
#' plot(dat,pch="+",cex=0.75,add=TRUE)
#' hist(abs(dat$y),xlab="Perpendicular distance",main="")
#' hist(dat$x,xlab="Forward distance",main="")
#' # do same with random n:
#' dat=sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi,fix.n=FALSE)
#' 
#' # compare time taken if calculate intscale on the run vs pass it:
#' # first calculate each time:
#' system.time(for(i in 1:20) dat<-sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi))
#' # then calculate once and pass:
#' intscale=calc.lpars(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#' system.time(for(i in 1:20) dat<-sim.n(n,ymin,ymax,W,h2,b,pi.norm,logphi,intscale=intscale))
#' }
#' @export
sim.n=function(n,ymin,ystart,w,hr,b,pi.x,logphi,fix.n=TRUE,intscale=NULL,nbuffer=NULL){
  simDat=list()
  simDat[[3]]=list(n=n,ymin=ymin,ystart=ystart,b=b,hr=hr,pi.x=pi.x,logphi=logphi,w=w,
                   fix.n=fix.n,intscale=intscale,nbuffer=nbuffer)
  names(simDat)[3]='settings'
  #I'd like to keep the function ARGS to keep them in line with the other functions-
  #- but I don't want to mess with the args inside the function.  So:
  # change args to David's
  ymax=ystart; W=w; hfun=hr
  
  # calculate scaling needed to give E[n]=n
  if(is.null(intscale)) 
    intscale=calc.lpars(n,ymin,ymax,W,hfun,b,pi.x,logphi)
  if(class(intscale)!="ppscale") stop("intscale must be class `ppscale' (output of calc.lpars).")
  window=owin(c(0,ymax),c(-W,W)) # create observation window
  lmax=intscale$lmax # maximum value of intensity function on grid used by calc.lpars
  lscale=intscale$lscale # multiplier required to get intensity with expected sampls size n
  if(is.null(nbuffer)){
    # increase E[n] by 25% to reduce prob that sample size is < n:
    nbuffer=ifelse(fix.n,1.25,1)
  }
  pp=rpoispp(poisint,lmax,window,ymin=ymin,ymax=ymax,b=b,hfun=hfun,pi.x=pi.x,logphi=logphi,W=W,lscale=lscale*nbuffer)
  if(fix.n) {
    while(pp$n<n){ # crude way of generating big enough sample size:
      pp2=rpoispp(poisint,lmax,window,ymin=ymin,ymax=ymax,b=b,hfun=h2,pi.x=pi.norm,logphi=logphi,W=W,lscale=lscale*nbuffer)
      pp$n=pp$n+pp2$n
      pp$x=c(pp$x,pp2$x)
      pp$y=c(pp$y,pp2$y)
    }
    pp$n=n
    pp$x=pp$x[1:n]
    pp$y=pp$y[1:n]
  }
  x=abs(pp$y);y=pp$x
  simDat[[1]]=data.frame(x=x,y=y);names(simDat)[1]='locs'
  simDat[[2]]=pp;names(simDat)[2]='pp'
  return(simDat)
}




#' @title Calculates scaling parameters required by \code{sim.n}
#' 
#' @description Calculates scaling parameters required by \code{sim.n} in 
#'              order to generate sample of given size.
#'              
#' @param n expected sample size reqired from NHPP
#' @param ymin smallest forward distance
#' @param ymax largest forward distance
#' @param W perpendicular truncation distance
#' @param hfun detection hazard function
#' @param b vector of detection hazard function parameters
#' @param pi.x perpendicular distance distribution function
#' @param logphi vector with log of pi.x parameters
#' @param nx number of intervals on x-axis at which to calculate intensity
#' @param ny number of intervals on y-axis at which to calculate intensity
#' @param inflate: multiplier by which to increase max intensity in region; this is 
#'        just a safeguard against the max at some point not calculated being
#'        greater than the points at which it was calculated. (The sample 
#'        generator in \code{rpoispp} uses rejection sampling so needs the
#'        global max.)
#'        
#' @details Calls \code{\link{poisint}} on grid of (x,y) points, calculates total
#' intensity in area and then calculates scaling needed to make this equal to 
#' n, and the maximum intensity in the region.
#' 
#' Output from this function can be passed as the \code{lscale} argument to
#' \code{\link{sim.n}}.
#' 
#' @return Object of class "ppscale", being alist with two parameters: 
#' \code{$lscale} is the multiplier needed to make the total intensity equal 
#' to n, while \code{$lmax} is maximum lintensity in the region. 
#' 
#'@examples
#' \dontrun{
#' n=100;ymin=0.01;ymax=5;W=2
#' b=log(c(0.75,1));logphi=c(0.5,log(0.3))
#' intscale=calc.lpars(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#' intscale
#' }
#' @export
calc.lpars=function(n,ymin,ymax,W,hfun,b,pi.x,logphi,nx=100,ny=100,inflate=1.05){
  x=seq(-W,W,length=nx)
  y=seq(ymin,ymax,length=ny)
  a=diff(x)[1]*diff(y)[1] # grid cell area
  lambda=outer(y,x,poisint,ymin=ymin,ymax=ymax,hfun=hfun,b=b,pi.x=pi.x,logphi=logphi,W=W)
  En=sum(lambda*a)
  lscale=n/En
  lmax=max(lambda*lscale)*inflate # bigger than observed in case higher between grid cells
  outlist=list(lscale=lscale,lmax=lmax)
  class(outlist)="ppscale"
  return(outlist)
}


#' @title NHPP intensity calculation in region
#' 
#' @description Calculates NHPP intensity in given region, using 
#'              given detection hazard and perp. dist. distribution.
#' @param y forward distances at which to calculate intensities (vector)
#' @param x perpendicular distances at which to calculate intensities (must
#'          be same length as y).
#' @param ymin smallest forward distance
#' @param ymax largest forward distance
#' @param W perpendicular truncation distance
#' @param hfun detection hazard function
#' @param b vector of detection hazard function parameters
#' @param pi.x perpendicular distance distribution function
#' @param logphi vector with log of pi.x parameters
#' @param lscale output of \code{\link{calc.lpars}} (object of class ``ppscale'')
#' 
#' @details Calculates survival model pdf of forward distance \code{y}, given 
#' perpendicular distance \code{x} (\code{f(y|x).}), multiplies this by the 
#' perpendicular distance pdf \code{pi.x} and then scales it by multiplying
#' by lscale (in order to get some total intensity - typically that to 
#' generate some expected sample size). See \code{\link{calc.lpars}} for 
#' details of the scaling.
#' 
#' This function is called by \code{rpoispp} inside \code{sim.n} to generate 
#' samples from NHPPs.
#' 
#' @return NHPP intensities at all \code{(x,y)}s input. 
#' 
#'@examples
#' \dontrun{
#' n=100;ymin=0.01;ymax=5;W=2
#' b=log(c(0.75,1));logphi=c(0.5,log(0.3))
#' nf=100
#' ys=seq(ymin,ymax,length=nf)
#' intscale=calc.lpars(n,ymin,ymax,W,h2,b,pi.norm,logphi)
#' f=poisint(ys,rep(0,nf),ymin=ymin,ymax=ymax,hfun=h2,b=b,pi.x=pi.norm,logphi=logphi,W=W,lscale=intscale$lscale)
#' plot(ys,f,type="l",xlab="Forward distance (y)",ylab="f(y)")
#' }
#' @export
poisint=function(y,x,ymin,ymax,hfun,b,pi.x,logphi,W,lscale=1){
  h=match.fun(hfun)
  pix=match.fun(pi.x)
  nx=length(x)
  ax=abs(x)
  if(length(y)!=nx) stop("Lengths of x and y must be the same.")
  pxx=pix(ax,logphi,W)
  f=p=p0=rep(NA,nx)
  for(i in 1:nx){
    #    p0[i]=1-Sy(ax[i],ymin,ymax,b,hfun)
    p[i]=1-Sy(ax[i],y[i],ymax,b,hfun)
    f[i]=h(y[i],ax[i],b)*(1-p[i])
  }
  #  return(list(f=f,p=p,p0=p0))
  return(f*pxx*lscale)
}


#' @title Calculates survivor function
#' 
#' @description Calculates survivor function to forward distance \code{y}, 
#' for given perpendicular distance \code{x} and given forward distance range.
#' 
#' @param x perpendicular distance (scalar)
#' @param y forward distance  (scalar)
#' @param ymax largest forward distance
#' @param hfun detection hazard function
#' @param b vector of detection hazard function parameters
#' @details Calculates probability of making it from forward distance 
#' \code{ymax} to \code{y} without being detected.
#' 
#' @return Probability of making it from forward distance 
#' \code{ymax} to \code{y} without being detected.
#' 
#'@examples
#' \dontrun{
#' ymax=5
#  b=log(c(0.75,1))
#' Sy(0,0.1,ymax,b,h2)
#' }
#' @export
Sy=function(x,y,ymax,b,hfun) {
  n=length(x)
  if(length(y)!=n) stop("Lengths of x and y must be the same.")
  pS=rep(NA,n)
  if(is.character(hfun)) {
    fName=hfun
    h=match.fun(hfun)
  }else {
    eval(parse(text=fNameFinder(hfun)))#this will create an object with the name fName that contains-
    #-the name of the hazard rate function passed into the Sy hfun ARG. The fName object should-
    #-only be available within the scope of the Sy function.
    #hchar=as.character(substitute(hfun))
    h=hfun
  }
  
  if(fName=="h1") { # Hayes & Buckland hazard rate model, so can do analytically
    hmax=h(y,x,b)
    for(i in 1:n){
      if(y[i]==0 | hmax[i]>1e10){ # computer will think integral divergent for large hmax
        pS[i]=1-HBhr(x[i],h1.to.HB(b))
      } else {
        pS[i]=exp(-integrate(match.fun(hfun),y[i],ymax,x=x[i],b=b,subdivisions = 1000L)$value)
        pS<<-pS[i]
      }
    }
  } else { # Not Hayes & Buckland hazard rate model, so can't do analytically
    for(i in 1:n){
      pS[i]=exp(-integrate(match.fun(hfun),y[i],ymax,x=x[i],b=b,subdivisions = 1000L)$value)
      pS2<<-pS[i]
    }
  }
  return(pS)
}
  

#' AIC-based model selection for models fitted using \link{fitxy}
#'
#'AIC-based model selection of a list of models fitted using \link{fitxy}.  Parameter 
#'estimates for the hazard and perpendicular distribution functions along with 
#'associated coefficients of variation are also provided. 
#'Models are ranked in order of increasing dAIC.  The model list can also be 
#'returned in a data frame suitable making a tex-type table.
#'@param modList a list of \link{fitxy}-type models.
#'@param modNames a vector of model names.
#'@param tab boolean - return a data frame suitable for generating a table for report or paper.
#'@param digits - number of digits of rounding (see \link{round}) in the table
#'@param maxblength=NULL; maximum number of parameters for the hazard function.  If NULL 
#'the maximum number of parameters for the models in modList is used.
#'@param maxlogphilength=NULL; maximum number of parameters for the perpendicular density function.  If NULL 
#'the maximum number of parameters for the models in modList is used.
#'@export
#'@return tab=FALSE list; element 1 =data frame of AIC-ranked models, with parameter estimates, 
#'coefficients of variation and AIC; element 2= vector of the order of AIC-ranked models. tab=TRUE, a three element list with AIC-ranked model order and two tables, the first as tab=FALSE, the 
#'second, a table for use in reports or manuscripts.
#'@seealso \link{fitxy}
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y 
#'fit1=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'fit2=fityx(y,x,b=log(c(0.001,1)),hr=h1,ystart,pi.x,logphi,w)
#'modSelect(modList=list(fit1,fit2),modNames=c('h1','h2'),tab=TRUE)
#'}
modSelect=function(modList,modNames=NULL,tab=FALSE,digits=2,
                   maxblength=NULL,maxlogphilength=NULL)
{
  if(is.null(maxblength))
    maxblength=max(sapply(modList,function(x) length(x$b)))
  if(is.null(maxlogphilength))
    maxlogphilength=max(sapply(modList,function(x) length(x$logphi)))
  
  fitVal=function(fit,maxblength=maxblength,
                  maxlogphilength=maxlogphilength)
  {
    bhat=CVbhat=rep(NA,maxblength)
    logphi=CVlogphi=rep(NA,maxlogphilength)
    bhat[1:length(fit$b)]=fit$b[1:length(fit$b)]
    CVbhat[1:length(fit$b)]=fit$CVpar[1:length(fit$b)]
    if(all(!is.na(fit$logphi)))
    {
      valLoc=which(!is.na(fit$logphi))
      logphi[valLoc]=fit$logphi[valLoc]
      CVlogphi[valLoc]=fit$CVpar[length(fit$b)+valLoc]
    }
    out=c(bhat,CVbhat,logphi,CVlogphi,length(fit$par),fit$value,fit$AIC)
    names(out)=c(paste('b',1:maxblength,sep=''),paste('CVb',1:maxblength,sep=''),
                 paste('logphi',1:maxlogphilength,sep=''),paste('CVlogphi',1:maxlogphilength,sep=''),
                 'n','logLik','AIC')
    return(out)
  }
  tab1=t(sapply(modList,fitVal,maxblength=maxblength,
                maxlogphilength=maxlogphilength))
  out1=data.frame(tab1)
  names(out1)=colnames(tab1)
  rownames(out1)=modNames
  minAICLoc=which.min(out1$AIC)
  out1$dAIC=out1$AIC-out1$AIC[minAICLoc]
  AICorder=order(out1$dAIC)
  out1=out1[AICorder,]
  ww=exp( -0.5 * out1$dAIC)
  out1$w=ww/sum(ww)
  if(tab)
  {
    bhat=logphihat=vector(mode='character',length=nrow(out1))
    for(i in 1:nrow(out1)){
      bhat[i]=paste(paste(round(tab1[i,1:maxblength],digits),'(',round(tab1[i,(maxblength+1):(2*maxblength)],digits),')',sep=''),collapse='; ')
      logphihat[i]=paste(paste(round(tab1[i,(2*maxblength+1):(2*maxblength+maxlogphilength)],digits),
                               '(',round(tab1[i,(2*maxblength+maxlogphilength+1):(2*maxblength*maxlogphilength)],digits),')',sep=''),collapse='; ')
    }
    out2=data.frame(bhat=bhat,logphihat=logphihat)
    out2=out2[AICorder,]
    out2$n=out1$n; out2$logLik=out1$logLik; out2$AIC=out1$AIC
    rownames(out2)=rownames(out1)
    out2$dAIC=out1$dAIC; out2$w=out1$w
    out2[,c('logLik','AIC','dAIC','w')]=round(out2[,c('logLik','AIC','dAIC','w')],digits)
    return(list(AICorder=AICorder,res=out1,tab=out2))
  }
  return(list(AICorder=AICorder,res=out1))
}

#'@title Calculates coverage probabilities of \eqn{\hat p}
#'
#'@description Calculate coverage probabilities of \eqn{\hat p} using the delta method and assuming a 
#'log-normal error distribution.
#'
#'@param fit object resulting from a call of \code{\link{fityx}} (see details)
#'@param interval the interval used to determine coverage probability
#'@param type \code{LOGNORM} log-normal confidence intervals; \code{NORM} normal confidence intervals.
#'@details In the call of \code{\link{fityx}} in the \code{fit} argument must have \code{hessian=TRUE}
#'@return a data frame of p phat, and CV[phat], interval, lower bound returned.
#'@examples
#'\dontrun{
#'ystart=4;w=1
#'hr=h2; b=log(c(0.75,1))
#'pi.x=pi.norm; logphi=c(0.5,log(0.2))
#'N=50 #true number of animals
#'#generate some observations
#'simDat=simXY(N=N,pi.x=pi.x,logphi=logphi,
#'hr=hr,b=b,w=w,ystart=ystart)
#'x=simDat$locs$x; y=simDat$locs$y 
#'est.yx=fityx(y,x,b,hr,ystart,pi.x,logphi,w)
#'phatInterval(fit=est.yx,interval=0.95)
#'}     
#'@seealso \code{\link{phat}} \code{\link{fityx}}
#'@export
phatInterval=function(fit,type='LOGNORM',
                      interval=0.95){
  lnci.nmin=function(stat,cv,stat.min=0,interval=interval){
    q=Mod(qnorm((1-interval)/2,0,1))
    varNhat=(stat*cv)^2
    cfactor=exp(q*sqrt(log(1+varNhat/(stat-stat.min)^2)))
    lower=stat.min+(stat-stat.min)/cfactor
    upper=stat.min+(stat-stat.min)*cfactor
    return(list(lower=lower,upper=upper))
  }
  
  if(!'hessian' %in% names(fit))
    stop('fit ARG must include a hessian matrix')
  pars=fit$par; 
  hr=match.fun(fit$hr); b=fit$b;
  ystart=fit$ystart; w=fit$w
  pi.x=match.fun(fit$pi.x); logphi=fit$logphi
  #estimated p
  p.hat=phat(fit=fit)
  #variance-covariance matrix
  vcov=fit$vcov
  #Implement the delta method:
  #numerical differentiation
  if(!is.numeric(fit$logphi))
  {dbyd=numericDeriv(quote(phat(w=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi)), c("b"))  } else {  
    dbyd=numericDeriv(quote(phat(w=w,b=b,hr=hr,ystart=ystart,pi.x=pi.x,logphi=logphi)), c("b","logphi"))  
  }
  dbyd=as.vector(slot(dbyd,'gradient'))
  var.p.hat=as.vector(t(dbyd)%*%vcov%*%dbyd) #$var[hat{p}(\hat{\Beta})$
  if(type=='LOGNORM')
    bounds=lnci.nmin(stat=p.hat,cv=sqrt(var.p.hat)/p.hat,interval=interval)
  if(type=='NORM'){
    bounds=list(lower=qnorm((1-interval)/2,p.hat,sqrt(var.p.hat)),
                upper=qnorm(interval+(1-interval)/2,p.hat,sqrt(var.p.hat)))
  }
  if(any(sapply(bounds,is.nan))){
    warning('One or both p.hat bounds NaN')
    return(data.frame(covered=NA,p=p,phat=p.hat,CV.phat=NA,
                      interval=interval,lower.bound=bounds$lower,upper.bound=bounds$upper))
  }
  
  return(data.frame(phat=p.hat,CV.phat=sqrt(var.p.hat)/p.hat,
                    interval=interval,lower.bound=bounds$lower,upper.bound=bounds$upper))
}

#'Calculate \eqn{\hat p} and optionally \eqn{\hat N} for a list of models
#'
#'Calculate \eqn{\hat p} along with variance, \eqn{Var[\hat p]}, using the delta method.  Optionally \eqn{\hat N} can also be calculated.
#'@param modList list object of models created by \link{fityx}.
#'@param n=NULL number of animals detected. If n!=NULL \eqn{\hat N} is calculated
#'@param tab boolean - return a data frame suitable for generating a table for report or paper.
#'@param digits - number of digits of rounding (see \link{round}) in the table
#'@param ... arguments to be passed into \link{phatInterval}
#'@export
#'@return data frame with:
#'\code{phat} estimate of \eqn{\hat p}    
#'\code{CV.phat} estimate of \eqn{CV[\hat p]} 
#'\code{interval} confidence interval specified in the \code{interval} argument in \link{phatInterval}
#'\code{lower.bound} lower bound of \eqn{\hat p}
#'\code{upper.bound} upper bound of \eqn{\hat p}      
#'and optionally if n!=NULL
#'\code{n} number of detected animals
#'\code{Nhat} estimated number of animals in covered region.
#'\code{NhatLower} lower bound of \eqn{\hat N} 
#'\code{NhatUpper} upper boudn of \eqn{\hat N}
#'@export
#'@seealso \link{fityx} \link{phatInterval}

phatModels=function(modList,n=NULL,tab=FALSE,digits=2,...)
{
  phatTab=sapply(modList,function(x) phatInterval(fit=x,...),...)
  colName=row.names(phatTab)
  phatTab=as.data.frame(t(matrix(as.numeric(phatTab),nrow(phatTab),ncol(phatTab))),
                        row.names=colnames(phatTab))
  names(phatTab)=colName
  if(!is.null(n))
  {
    phatTab$n=n
    phatTab$Nhat=n/phatTab$phat
    phatTab$NhatLower=n/phatTab$upper.bound
    phatTab$NhatUpper=n/phatTab$lower.bound
  }
  if(!tab){
    return(phatTab)}else{
      phatV=vector(length=length(modList))
      for(i in 1:length(modList))
        phatV[i]=paste(paste(round(phatTab$phat[i],digits),
                             '(',round(phatTab$CV.phat[i],digits),')',sep=''),collapse='; ')
      tab=data.frame(phat=phatV,row.names=row.names(phatTab))
      if(!is.null(n)){
        NhatV=vector(length=length(modList))
        for(i in 1:length(modList))
          NhatV[i]=paste(paste(round(phatTab$Nhat[i],0),
                               '(',round(phatTab$NhatLower[i],0),',',
                               round(phatTab$NhatUpper[i],0),
                               ')',sep=''),collapse='; ')
        tab=cbind.data.frame(tab,Nhat=NhatV)
      }
      return(list(res=phatTab,tab=tab))
    }
}

#' Plot the 2D fit of a model
#' 
#' Plot the 2D fit of a model resulting from a call of \link{fityx}
#' @param fit object resulting from a call of \link{fityx}
#' @param ... other parameters passed into \link{plotSim}
#' @details This function is a wrapper for \link{plotSim}.
#' @seealso \link{plotSim} \link{fityx}
plotFit=function(fit,...){
  obj=list(locs=data.frame(x=fit$dat$x,y=fit$dat$y),
           settings=list(pi.x=match.fun(fit$pi.x),
                         logphi=fit$logphi,
                         hr=fit$hr,
                         b=fit$b,
                         w=fit$w,
                         ystart=fit$ystart))
  plotSim(simDat=obj, nclass=10,xlab="perpendicular distance (x)", 
          ylab="forward distance (y)",image=TRUE,...)      
  
}

#'Find a function name when the function is passed as an argument into other functions.
#'
#'@param x an R function/method that has not been evaluated.  
#'@param fSearchSting='fName' function name search string
#'@details The function passed into x the x argument must have the object \code{fSearchSting} 
#'within its body.  The fSearchSting must be the only code on a line within the function,
#'but can be placed anywhere in the function body, e.g. \code{fName='f1'}. See examples.
#'@returns a character string e.g. "fName='f1'" that can be evaluated using 
#'\code{eval(parse(text="fName='f1'") } either outside of a function, see example or within the scope of a function.
#'@export
#'@examples
#'##Example 1 - check this approach works with nested functions
#'f1=function(x) {
#'fName='f1'
#'x**2}
#'
#'f2=function(f,y){
#'  res=f(x=y)
#'  funcName1=fNameFinder(f)
#'  return(list(res,funcName1))}
#'
#'f3=function(f,y){
#'  res=f2(f,y)
#'  return(res)}
#'
#'eval(parse(text=f3(f=f1,y=2)[[2]]))
#'
fNameFinder=function(x,fSearchSting='fName'){
  x=as.character(attributes(x)[[1]]) #get attributes of function x
  x=gsub(" ", "",x[grep(fSearchSting,x)], fixed = TRUE) #search for line in function body with 
  return(x)}




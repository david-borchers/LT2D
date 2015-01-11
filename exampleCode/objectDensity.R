#'Object density at a given distance and angle from the observer
#'
#'This function calculates the density of objects at a given distance and angle from the observer.
#'
#'Object density can be modelled using the following parametric forms
#'(\code{f} argument): and parameters passed to the function using the
#'\code{pars} argument.
#'The expoential \code{f='Exp'} parameters \code{pars} is a three element vector defining the 
#'exponential relationship \eqn{\beta[0]} and \eqn{\beta[1]}, and the power by which 
#'radial distance \code{r} is raised to. Typically radial distance is squared 
#'giving \code{pars[3]=2}.
#'@param r radial distance from observer to detected object.
#'@param theta angle from observer to detected object. Not used.
#'@param f parametric form of the density function.
#' Currently only exponential c('Exp') is implemented. See details.
#'@param pars vector of detection function parameters.
#'@return Density of objects at range \code{r} and angle \code{theta} from an observer.

#'@examples
#'objectDensity(r=50,theta=NULL,pars=c(0.09,0.03,2))
#' @export
objectDensity=function(r,theta=NULL,f='Exp',pars){
  switch(f,
         Exp=exp(pars[1]+pars[2]*r*pars[3]))
}
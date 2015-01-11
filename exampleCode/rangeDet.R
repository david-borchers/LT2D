#'Probability of detection at a given range
#'
#'This function calculates the probability of detecting an object at a given distance from the observer.
#'
#'The probability of detecting a target at a given radial distance from an
#'observer can be modelled using the following parametric forms
#'(\code{f} argument): and parameters passed to the function using the
#'\code{pars} argument.
#'The half-normal \code{f='hn'} parameters \code{pars} is a single value of sigma.
#'The uniform detectability has a single value, the constant 
#'probability of detection. Perfect detectability is specified by \code{pars=1}.
#'@param r radial distance from observer to detected object.
#'@param f parametric form of the detection function. See details.
#'@param pars vector of detection function parameters.
#'@return Detection probability at range \code{r}.
#'@section Reference: Buckland, S.T., Anderson,D. R., Burnham, K. P., Laake, J.
#'L., Borchers,D. L. and Thomas, L. (2001) Introduction to Distance Sampling.
#'Oxford: Oxford University Press, page 47.
#'@examples
#'rangeDet(r=50,pars=25)
#'rangeDet(r=50,f='unif',pars=1)
#' @export
rangeDet=function(r,f='hn',pars=NULL){
  switch(f,
         hn=exp(-r^2/(2*pars[1]^2)),
         unif=rep(pars[1],length(r)))
}
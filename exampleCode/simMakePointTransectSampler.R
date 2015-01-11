library(sp)
library(sampSurf)
library(rgeos)
library(raster)
library(polyclip)
library(animation)
#' Set-up a point transect sampler
#' Set-up a single point transect sampler as a SpatialPolygons object using \code{\link{}}
#' 
#' @param w numeric; truncation distance
#' @param loc numeric vector; x,y coordinates
#' @param ... other arguments to be passed to \code{\link{SpatialPolygons}} and \code{\link{spCircle}}
#'@seealso \code{\link{spCircle}}
#'@return a circle as a SpatialPolygons object centered on \code{loc}, of radius \code{w}
#'@examples
#'g=simSetupRegion(cells.dim=c(20,20))
#'pt=simMakePointTransectSampler(w=5,loc=c(10,10))
#'plot(g)
#'plot(pt,add=TRUE,col='purple')
#'@export
#'
simMakePointTransectSampler=function(w,loc,...){
  TMPcircle=spCircle(radius=w,centerPoint=c(x=loc[1],y=loc[2]),...)
  circle=slot(TMPcircle$spCircle,'polygons')$pgsCircle
return(SpatialPolygons(list(circle),...))
}
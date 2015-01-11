#'Set-up a survey region to simulate from
#'
#'This function sets up a square or rectangular survey region to simulate from.
#'
#'The survey region is simulated using a spatial grid, so can be assigned a coordinate system.
#'@seealso \code{\link{sp}} \code{\link{GridTopology}}
#'@param cellsize numeric; vector with the cell size in each dimension (x,y)
#'@param cells.dim numeric; vector with the cell size in each dimension
#'@param ... other argumments passed into \code{\link{GridTopology}} 
#'@return spatial polygon grid object.
#'@examples
#'simSetupRegion(cells.dim=c(30,20))
#'@export
simSetupRegion=function(cellsize=c(1,1),cells.dim,...){
  BaseGrid = GridTopology(cellcentre.offset=cellsize/2, 
                          cellsize=cellsize, 
                          cells.dim=cells.dim,...)
  SpP_grd <- as.SpatialPolygons.GridTopology(BaseGrid,...)
  return(SpP_grd)
}

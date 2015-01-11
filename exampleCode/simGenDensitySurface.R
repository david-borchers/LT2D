#' Create a simulated density surface within a survey region
#' 
#' This function is currently set-up to generated a density surface for point transects.  
#' The density varies from the observer.  This code is currently used to simulate confounding
#' between the density surface and the detection function.
#' 
#' @param simRegion SpatialPplygon; Simulated survey region typically created using \code{\link{simSetupRegion}}
#' @param samplerLoc real, two-element vector; x,y coordinates of the centre of a point transect sampler.
#' @param denType character, 
#' @param ... other arguments for \code{\link{objectDensity}}. Not used.
simGenDensitySurface=function(simRegion,samplerLoc,denType='Exp',denPars,...)


makeTrapSpPolyList=function(x,y,r){
  #depends on spCircle in sampSurf package
  x.store=x;y.store=y
  if(length(r)==1)
    r=rep(r,length(x))
  circleL=list()
  for(i in 1:length(x.store))
  {
    x=x.store[i];y=y.store[i]
    TMPcircle=spCircle(radius=r[i],centerPoint=c(x=x,y=y))
    circleL[[i]]=slot(TMPcircle$spCircle,'polygons')$pgsCircle
  } #end of i for loop
  return(SpatialPolygons(circleL,1:i))
} #end of makeTrapSpPolyListed from top left to bottom right.



dDF=data.frame(survey.overlap = rep(0, length(SpP_grd)),
               camera.overlap=rep(0, length(SpP_grd)),
               camera.cell.prop=rep(0, length(SpP_grd)),
               svy.overlap=rep(0, length(SpP_grd)),
               svy.cell.prop=rep(0, length(SpP_grd)),
               row.names=sapply(slot(SpP_grd, "polygons"), function(i) slot(i, "ID")))

g <- SpatialPolygonsDataFrame(SpP_grd,dDF, match.ID = TRUE)
@
We can now plot (\ref{fig:one}) the survey grid and add the unique polygon (grid cell) identifier:
  <<label=plotGrid,eval=FALSE>>=
  plot(g)
text(coordinates(SpP_grd), sapply(slot(SpP_grd, "polygons"), 
                                  function(i) slot(i, "ID")), cex=0.5)
@
\begin{figure} [ht!]
\begin{center}
<<label=fig1,fig=TRUE,echo=FALSE>>=
  <<plotGrid>>
  @
\end{center}
\caption{Survey grid with polygon (cell) identifiers.}
\label{fig:one}
\end{figure}  

\section{Add survey boundary and camera locations to survey grid}
In this section we create a function, \texttt{makeTrapSpPolyList} to set up spatial objects defining camera locations (x,y) and their sampling radii (r). Camera locations are specified as a vcetor of x locations and a vector of y locations.  Sampling radii can be specified as a single common value (across all cameras) or a vector of individual camera sampling radii. \texttt{makeTrapSpPolyList()} returns a spatial list object containing circles for each camera location: 
  <<label=makeTrapSpPolyListFunc,tidy=TRUE>>=
  makeTrapSpPolyList=function(x,y,r){
    #depends on spCircle in sampSurf package
    x.store=x;y.store=y
    if(length(r)==1)
      r=rep(r,length(x))
    circleL=list()
    for(i in 1:length(x.store))
    {
      x=x.store[i];y=y.store[i]
      TMPcircle=spCircle(radius=r[i],centerPoint=c(x=x,y=y))
      circleL[[i]]=slot(TMPcircle$spCircle,'polygons')$pgsCircle
    } #end of i for loop
    return(SpatialPolygons(circleL,1:i))
  } #end of makeTrapSpPolyListed from top left to bottom right.
@
Using \texttt{makeTrapSpPolyList()} we will create a spatial polygons object describing the camera locations.  

\section{Assign animal home ranges}
Before we set up the sampler locations, we will look at setting up spatial polygon offset from the survey boundary to help us account for edge effects.
<<offsetSvyEdge1>>=
  unitsToOffset=0 #negative denotes inside the survey boundary.  
@
Next we will create a polygon based on the survey boundary object:
  <<offsetSvyEdge2>>=
  offsetObj<-list(list(x=c(extent(g)@xmin,extent(g)@xmax,extent(g)@xmax,extent(g)@xmin,extent(g)@xmin),
                       y=c(extent(g)@ymin,extent(g)@ymin,extent(g)@ymax,extent(g)@ymax,extent(g)@ymin)))

Cobj <- polyoffset(offsetObj, unitsToOffset, jointype="round")

CP=Polygon(coords=cbind(c(Cobj[[1]]$x,Cobj[[1]]$x[1]),c(Cobj[[1]]$y,Cobj[[1]]$y[1])))
bndSPPolyClip<-SpatialPolygons(list(Polygons(list(CP),ID=1)))
@
Now we have set-up the survey boundary, we will randomly allocate home ranges for $N =$ \Sexpr{N} animals.
<<homeRng>>=
  animalHome=spsample(x=bndSPPolyClip,n=N,type='random')

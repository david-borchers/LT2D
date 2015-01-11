\documentclass{article}

\begin{document}
\SweaveOpts{concordance=TRUE}
Load packages into \texttt{R} workspace:
  <<loadPackages>>=
  @
  Start by analysing Wallaby creek data found in the files 'WC86TN' for
transects and 'WC86IS' for sightings:
  <<readDat>>=
  transects=read.csv('c:\\Users\\martin_cox\\Documents\\2Ddsitext\\dat\\WC86TN.csv')
sightings=read.csv('c:\\Users\\martin_cox\\Documents\\2Ddsitext\\dat\\WC86IS.csv')
@
  and source functions:
  <<sourceFuncation>>=
  source('c:\\Users\\martin_cox\\Documents\\2Ddsitext\\scatterhist.R')
@
  
  Examine transect data:
  <<transectDat>>=
  summary(transects)
@
  Have a look at transect lengths:
  <<tNbr>>=
  hist(transects$TLGT)
@
  and transect durations
<<transectDuration>>=
  hist(transects$TDUR)
with(transects,plot(TDUR,TLGT))
@
  and the direction of transects:
  <<transectDir>>=
  table(transects$TBRG)
@
  merge the sightings and transect data:
  <<mergeDat>>=
  sightings=merge(sightings,transects,'TNNU')
@
  \subsection{Sightings data}
<<sightingsData>>=
  summary(sightings)
@
  get a y-coordinate for each sightigns:
  <<Ycoord>>=
  sightings$Y=sqrt(sightings$RADL**2-sightings$PERP**2)
@
  
  <<allSightings>>=
  with(sightings, scatterhist(PERP, Y, xlab="X",
                              ylab="Y",col=as.numeric(as.factor(TDRN)),pch=19))
@
  
  \end{document}


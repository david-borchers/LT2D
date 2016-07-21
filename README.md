
<!-- README.md is generated from README.Rmd. Please edit that file -->
### LT2D: An R package for analysing line transect sampling data in two dimensions

By formulating distance sampling models as survival models, the `2DLT` package uses time to first detection in addition to perpendicular distance (line transect surveys) or radial distance (point transect surveys) allows estimation of detection probability, and hence density, when animal distribution in the vicinity of lines or points is not uniform and is unknown.

Conventional distance sampling (CDS) methods assume that animals are uniformly distributed in the vicinity of lines or points. But when animals move in response to observers before detection, or when lines or points are not located randomly, this assumption may fail. By formulating distance sampling models as survival models, we show that using time to first detection in addition to perpendicular distance (line transect surveys) or radial distance (point transect surveys) allows estimation of detection probability, and hence density, when animal distribution in the vicinity of lines or points is not uniform and is unknown.

-   The package incorporates times to detection and can provide information about failure of the conventional distance sampling (CDS) assumption that detection probability is 1 at distance zero.

-   The maximum likelihood estimator (MLE) of line transect survey detection probability and effective strip half-width using times to detection has been implemented in the package (see for example `LT2D::fityx`).

-   The properties of this MLE have been investigated by simulation in situations where animals are nonuniformly distributed and their distribution is unknown. The MLE is found to perform well when detection probability at distance zero is 1.

-   The MLE allows unbiased estimates of density to be obtained in this case from surveys in which there has been responsive movement prior to animals coming within detectable range.

-   The perforamnce of the MLE is illustrated using two real datasets: firstly by estimating primate density from a line transect survey in which animals are known to avoid the transect line, and secondly using data from a shipboard survey of dolphins that are attracted to the ship.

### Installing LT2D

You can install the `LT2D` package directly from github using `devtools` (Wickham & Winston, 2016):

``` r
if(!"devtools" %in% rownames(installed.packages())) {install.packages("devtools")}
devtools::install_github('david-borchers/LT2D')
library(LT2D) 
```

Please cite \`LT2D' when you use it:

``` r
citation('LT2D')
#> 
#> To cite 2DLT in publications use:
#> 
#>   Borchers DL and Cox MJ, (2015). Distance sampling detection
#>   functions: 2D or not 2D? Biometrics, -:--, 10.1111/---
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {Distance sampling detection functions: 2D or not 2D?},
#>     author = {David L Borchers and Martin J Cox},
#>     journal = {Biometrics},
#>     year = {2016},
#>     volume = {-},
#>     number = {-},
#>     doi = {10.1111/---},
#>   }
```

### Example anlayses and R code

In the following two subsections we provide example code to run some the analyses presented in our paper. There is a section for each of the datasets, primate that aovid the observer and dolphins that are attracted to the observer. The data are distributed with the package, see:

``` r
help(primate.dat)
help(dolphin.dat)
```

### Primate data analysis

We are grateful to Matthew Nowak from the Sumatran Orangutan Conservation Programme (SOCP) for allowing us to use the primate survey data from the Jantho Reintroduction Station. The initial survey was developed by Matthew Nowak and Serge Wich (Liverpool John Moores University) and then undertaken by the SOCP with funding from Chester Zoo.

The primate line transect data were used to illustrate joint estimation of a hazard function, \(h\), and prependicular density distribution, \(\pi\), when animals exhibit avoidance behaviour. In the code chuck below we plot the sightings \(x,y\) positions and the distribution in the \(x\) (perpendicular distance) and \(y\) (forward distance).

A perpendicular truncation distance, \(w\), of and a maximum forward distance, \(y_{max}\), of were used in our primate data analyses. For line transect notation see Fig. 1 in the paper.

``` r
data(primate.dat)
x=primate.dat$x
y=primate.dat$y

openGraph(h=3,w=12)
par(mfrow=c(1,3))
pdlab="Perpendicular distance"
fdlab="Distance along transect"
plot(jitter(y,1,0),jitter(x),pch="+",ylab=pdlab,xlab=fdlab,main="")
hist(y,breaks=seq(0,max(na.omit(y)),length=16),xlab=fdlab,main="")
hist(x,breaks=seq(0,max(na.omit(x)),length=12),xlab=pdlab,main="")
```

We selected the hazard and perpendicular density functions using AIC. The following code chuck is code to run the estimate parameters for the selected model that has a normally distributed prependicular density function, \(\pi_N\), and an inverse power hazard function, \(h_{IP0}\), in which the scale parameter fixed, \(\beta_2\), at 1.

``` r
# Normal bump with ip0 hazard function (SELECTED MODEL):
b=c(5, 8)
logphi=c(0.02, -4.4)
fit.n.ip0=fityx(y[x<=w],x[x<=w],b=b,hr=ip0,ystart=ystart,
                pi.x=pi.norm,logphi=logphi,w=w,hessian=TRUE,control=list(trace=5,maxit=1000))
```

The following code chuck will reproduce the primate fits from the paper (Fig. 3) with the top row being the perpendicular distance, \(x\), the bottom row, forward distance, \(y\). Histograms of observations are given in the first column along with fitted PDFs as solid lines. The right column are Q-Q plots in each dimension and Cramer-von Mises and Kolmogorov–Smirnov goodness-of-fit statistic p-values are available from calls of `GoFx` and `GoFy`.

``` r
openGraph();par(mfrow=c(2,2))
plotfit.x(x[x<=w],fit.n.ip0,nclass=20);rug(x[x<=w])

GoFx(fit.n.ip0,plot=TRUE)$pvals

plotfit.smoothfy(fit.n.ip0,nclass=32);rug(x=y[x<=w])

GoFy(fit.n.ip0,plot=TRUE)$pvals
```

The effective strip half-width, \(\tilde{w}\), is estimated using `phatInterval`:

``` r
phatInterval(fit.n.ip0)*w
```

The proportion of primates detected at distance 0, \(p(0)\), is calculated by:

``` r
1-Sy(0,0,ystart,fit.n.ip0$b,ip0)
```

### Dolphin data analysis

We are gerateful to North Atlantic Marine Mammal Commission (NAMMCO) and the Faroese Museum of Natural History for allowing us to use the dolphin survey data from 1995 North Atlantic Sightings Survey (NASS95). These data are analysed using mark-recapture distance sampling methods by Ca~nadas et al. (2004).

The dolphin line transect data were used to illustrate joint estimation of a hazard function, \(h\), and prependicular density distribution, \(\pi\), when animals exhibit attraction . In the code chuck below we plot the sightings \(x,y\) positions and the distribution in the \(x\) (perpendicular distance) and \(y\) (forward distance).

A perpendicular truncation distance, \(w\), of and a maximum forward distance, \(y_{max}\), of were used in our primate data analyses.

``` r
data(dolphin.dat)
xd=dolphin.dat$x
yd=dolphin.dat$y

openGraph(h=3,w=12)
par(mfrow=c(1,3))
pdlab="Perpendicular distance"
fdlab="Forward distance"
plot(jitter(yd[xd<=wd],1,0),jitter(xd[xd<=wd],1,0),pch="+",ylab=pdlab,xlab=fdlab,main="")
hist(yd[xd<=wd],breaks=seq(0,max(yd[xd<=wd]),length=26),xlab=fdlab,main="")
hist(xd[xd<=wd],breaks=seq(0,max(xd[xd<=wd]),length=19),xlab=pdlab,main="")
```

Again, we selected the hazard and perpendicular density functions using AIC. The following code chuck is code to run the estimate parameters for the selected model that has a half-normal prependicular density function, \(\pi_HN\), and the hazard function of Hayes and Buckland (1983), \(h_{HB}\).

``` r
b=c(-7.3287948, 0.9945317)
logphi=-0.4811025
dfit.hn=fityx(yd[xd<=wd],xd[xd<=wd],b=b,hr=h1,ystart=ystartd,
              pi.x=pi.hnorm,logphi=logphi,w=wd,hessian=TRUE,control=list(trace=5))
```

The following code chuck will reproduce the dolphin fits from the paper (Fig. 4) with the top row being the perpendicular distance, \(x\), the bottom row, forward distance, \(y\). Histograms of observations are given in the first column along with fitted PDFs as solid lines. The right column are Q-Q plots in each dimension and Cramer-von Mises and Kolmogorov–Smirnov goodness-of-fit statistic p-values are available from calls of `GoFx` and `GoFy`.

``` r
openGraph();par(mfrow=c(2,2))
plotfit.x(xd[xd<=wd],dfit.hn,nclass=20);rug(xd[xd<=wd])
GoFx(dfit.hn,plot=TRUE)$pvals
plotfit.smoothfy(dfit.hn,nclass=32);rug(x=yd[xd<=wd])
GoFy(dfit.hn,plot=TRUE)$pvals
```

The effective strip half-width, \(\tilde{w}\), is estimated using `phatInterval`:

``` r
phatInterval(dfit.hn)
phatInterval(dfit.hn)*wd
```

We estimated dolphin density using \(\hat{D}= n \times \hat{p}/2\times w \times L\) :

``` r
n=length(dfit.hn$dat$x)
L=1672.77 #nautical miles from Canadas et al. (2004)
Dhat=(n/phatInterval(dfit.hn)$phat)/(2*wd*L)
Dhat
```

R code to run all the model fits in the paper is available in the file 'FitsForPaper.r'.

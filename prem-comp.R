## Make plot of PREM and calculated wavespeeds with liquid metal model
##   G. Helffrich, 15 Oct. 2015

source('premfunc.R')
source('yokomodel.R')

if (!exists('zvals')) zvals <- seq(0,3480-1221) ## Default resolution for plot

gam.vals <- c(1.3,1.8)      ## Range of gamma values to explore

zY <- 1.649; mwY <- 48.8    ## Yokoyama parameters
Tcmb <- 4300                ## Estimated CMB T
T<-Tcmb * tad(zvals)        ## Temperature increase below CMB

## Function for depth-dependent z (valence electron count) that matches PREM
##   exactly.

zf <- function(h){
   rep(zY,length(h)) +
   (1.353967e-02 + 8.051915e-05*h - 2.214906e-07*h^2 + 9.533911e-11*h^3) /
      (1 - 4.369232e-03*h + 7.911845e-06*h^2 - 1.693164e-09*h^3)
}

## Comment out redefinition of zf below if you want to see the exact match.
## Otherwise, the plot shows the simple fit to PREM by varying
##   Tcmb, z, and mol. wt.

zf <- function(h)rep(zY, length(h))

spar<-par('mar'); par(mar=0.1+c(5,5,4,2)) ## Room for fancy label on left margin

plot(3480-zvals,vprem(3480-zvals),type='l',lty=5,lwd=2,
   xlab='R (km)',
   ylab=expression(list(V[P]~(km~s^-1)))
)

dtvals <- c(0,-200,200)         ## Temperature varies this much from Tcmb

for(dt in dtvals){
   lines(3480-zvals,
      sapply(zvals,
         function(h){vyoko(pgpa(3480-h),(Tcmb+dt)*tad(h),rhoprem(3480-h),
	    z=zf(h),sp='Fe',mw=mwY)$c}
      )/1000,
      lty=1
   )
   if (dt == 0) {
      ## Change in gamma does nothing due to low T dependence
      par(cex=1.25)
      text(1900, 8.50, 'PREM',adj=c(0,0))
      text(1450, 8.35, expression(T[cmb] == 4300~K),adj=c(0,0))
      text(1450, 8.2, expression(paste(1.3 <= gamma) <= 1.8),adj=c(0,0))
      par(cex=1.0)
      lines(c(1450,1850),rep(8.53,2),lty=5,lwd=2)
   }
}

par(mar=spar)

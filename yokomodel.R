## R routines to implement Yokoyama liquid metal wavespeed model as described
## in Helffrich (2015), The hard sphere view of the outer core,
##   Earth, Planets, Space, v.67, 73-84.

N0 <- 6.0221413e23

spdat <- list(          ## List of metal properties for reproduction of
                        ## values in Yokoyama (2001), Mat. Trans. v42, 2021ff.
   Fe=list(mw=55.845,    rho=7010, Tm=1833, zm=1.33),
   Ni=list(mw=58.6394,   rho=7720, Tm=1773, zm=1.30),
   Ti=list(mw=47.867,    rho=4150, Tm=1973, zm=1.47),
   V= list(mw=50.9415,   rho=5360, Tm=2173, zm=1.49),
   Cr=list(mw=51.9961,   rho=6270, Tm=2173, zm=1.46),
   Mn=list(mw=54.938045, rho=5970, Tm=1533, zm=1.25),
   Co=list(mw=58.933195, rho=7700, Tm=1823, zm=1.32)
)

Bh <- function(z, rhom, Tm, mw=55.847, xim=0.463,
   R=8.3144, l=5.2917725e-11, Ryd=2.1798741e-18){    ## Ryd = N0*Rydberg energy
   ## xim is like the Protopapas et al. melting value, but by Yokoyama (1999)
   if(z<=0) return(0)
   n <- N0 * rhom*1e3/mw
   a <- (3/(4*pi*n))^(1/3)/l; zz <- z^(1/3)
   a^6/6 * ( Ryd*(0.031*z + z*(0.916*zz + 1.8*z)/a - 4.42*zz^5/a^2) -
      (1 + xim + xim^2 - xim^3)/(1 - xim)^3 * 3*1.3806488e-23*Tm
   )
}

zmod <- function(T, rho, P, mw=55.847, z=1.33, sp=spdat[['Fe']],
   R=8.3144, l=5.2917725e-11, Ryd=2.1798741e-18){    ## Ryd = N0*Rydberg energy
   n <- rho*1e3/mw
   a <- (3/(4*pi*n*N0))^(1/3)/l; zz <- z^(1/3)
   P*1e9/(n*R*T) + 1/(3*1.3806488e-23*T) * (Ryd*(0.031*z + z*(0.916*zz + 1.8*z)/a -
      4.42*zz^5/a^2) - 6/a^6*Bh(z, sp$rho, sp$Tm, mw=sp$mw))
}

cmod <- function(T, rho, P, mw=55.847, z=1.33, sp=spdat[['Fe']],
   R=8.3144, l=5.2917725e-11, Ryd=2.1798741e-18){     ## Ryd = N0*Rydberg energy
   n <- rho*1e3/mw
   a <- (3/(4*pi*n*N0))^(1/3)/l; zz <- z^(1/3)
   (Ryd*(-0.031*z - 4/3*z*(0.916*zz + 1.8*z)/a + 22.1/3*zz^5/a^2) +
      18/a^6*Bh(z, sp$rho, sp$Tm, mw=sp$mw)) / (3*1.3806488e-23*T)
}

dsig <- function(phi, T, Tm=1833){
   ## sig <- 6*phi/pi; sig0 <- 1.0878
   ## -0.056*(sig0/sig)^(1/3)*sqrt(T/Tm)
   -0.056*sqrt(T/Tm)/(1.112*0.888*(1-0.112*sqrt(T/Tm)))
}

Phi <- function(P, T, rho, z=1.33, mw=55.847, sp=spdat[['Fe']], ds=0.1){
   zp <- function(phi){(1 + phi + phi^2 - phi^3)/(1-phi)^3}
   zm <- zmod(T, rho, P, z=z, mw=mw, sp=sp)
   if (zm < 1) return(1e-10)
   for(s in seq(0,1,ds)){
      if (zm - zp(s) < 0) break
   }
   ## cat('z,eta:',zm,s,'\n')
   uniroot(function(sig){zm-zp(sig)},c(s-ds,s),tol=1e-6)$root
}

vyoko <- function(P, T, rho, z=NA, mw=NA, sp='Fe', R=8.3144){
   ## P in GPa, T in K, rho in kg/m^3
   ## Liquid metal species assumed to be Fe
   ## if z, mw not explicitly given, these default to species properties 
   if(is.na(z)) z <- spdat[[sp]]$zm
   if(is.na(mw)) mw <- spdat[[sp]]$mw

   Z <- zmod(T, rho, P, z=z, mw=mw, sp=spdat[[sp]])
   if (Z<=1){
      phi <- 0.463
   } else {
      phi <- Phi(P, T, rho, z=z, mw=mw, sp=spdat[[sp]])
   }
   p <- (1 + phi + phi^2 - phi^3)/(1 - phi)^3
   pp <- 2*(2 + 2*phi - phi^2)/(1 - phi)^4
   dsigdt <- dsig(phi,T,Tm=spdat[[sp]]$Tm)
   vel <- sqrt(1000/mw * R*T*(
      cmod(T, rho, P, z=z, mw=mw, sp=spdat[[sp]]) +
      p + phi*pp + 2/3*(p + 3*phi*pp*dsigdt)^2)
   )
   list(c=vel, eta=phi, Z=Z)
}

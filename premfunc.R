## Various routines encoding properties of the PREM model and other
##   seismic wavespeed models (iasp91, AK135) for comparison.
## G. Helffrich, 2008.
##   last change 16 Oct. 2016

## Reference:  Dziewonski, A. and Anderson, D. (1981), Preliminary Reference
##   Earth Model, PEPI v. 25, 297-356.

rhoprem<-function(r,re=6371){
   ## Outer core density values from PREM
   eval<-function(x,a,b,c,d) {((d*x + c)*x + b)*x + a}
   rn <- r/re
   1000*eval(rn,12.5815,-1.2638,-3.6426,-5.5281)
}

vprem<-function(r){
   ## PREM wavespeeds in outer core.
   g<-function(x,a,b,c,d){((d*x + c)*x + b)*x + a}
   sapply(r/6371, g, a=11.0487, b=-4.0362, c=4.8023, d=-13.5732)
}

g<-function(r){
   ## Gravity in the PREM model.
   ##   r is in km; return g in m/s^2. Fit to PREM density model.
   rep(2.114073e-03,length(r)) + 3.818921e-03*r -2.131285e-07*r^2
}

pgpa<-function(r){
   ## Pressure in the PREM model.
   ##   r is in km; return P in GPa. Fit to PREM density model.
   rep(3.853676e+02,length(r)) - 3.083633e-02*r - 1.184717e-05*r^2
}

vkhoc<-function(r,rcmb=3480){
   ## KHOC wavespeeds in the outer core
   lr <- length(r)
   rb<-rcmb-c(150,250,450); sb<-0.01*c(+1,-1,-1); pb<-c(3.910645, 2.768334, 3.409311)
   dvp<-rep(0,lr)
   for(j in 1:lr){
      for(i in 1:3) {
         dvp[j]<-dvp[j]+pb[i]*sb[i]*(max(0,r[j]-rb[i])/(rcmb-rb[i]))^2
      }
   }
#  for(i in 1:3) {
#     dvp<-dvp+pb[i]*sb[i]*(pmax(rep(0,lr),r-rep(rb[i],lr))/(rcmb-rb[i]))^2
#  }
   dvp+vprem(r)
}

vsr<-function(r,ricb=1221.5){
   ## Souriau & Roudil (1995) wavespeeds in the outer core
   ##   Reference:  Souriau, A., Roudil, P. (1995).  Attenuation in the
   ##   uppermost inner core from broad-band GEOSCOPE PKP data.  GJI v. 123,
   ##   572-587.
   if(!exists('vsr.bnd')){
      vsr.bnd <<- uniroot(function(r){10.27-vprem(r)},c(ricb+100,ricb+200))$root
   }
   ifelse(r<vsr.bnd,rep(10.27,length(r)),vprem(r))
}

vspr<-function(r,ricb=1221.5){
   ## Ohtaki et al. (2012) south polar wavespeed in the outer core.
   ## Reference:  Ohtaki, T., Kaneshima, S., Kanjo, K. (2012).  Seismic
   ##   structure near the inner core boundary in the south polar region.
   ##   JGR v. 117, B03312, doi:10.1029/2011JB008717.
   if(!exists('vspr.bnd')){
      vspr.bnd <<- uniroot(
         function(r,val=vprem(ricb)-0.04){vprem(r)-val},c(ricb+50,ricb+100)
      )$root
   }
   ifelse(r<vspr.bnd,rep(vprem(ricb)-0.04,length(r)),vprem(r))
}

viasp91<-function(r, xn=6371){
   ## iasp91 in the outer core.
   x <- r/xn
   sapply(x, function(x){10.03904 + x*(3.75665 - x*13.67046)}
   )
}
## data p/11.24094, 10.03904, 14.49470, 25.1486, 25.969838, 29.38896,
##        30.78765, 25.41389, 8.785412, 6.5   , 5.8   , 0., 0.
##        0.   ,    3.75665,  -1.47089,-41.1538,-16.934118,-21.40656,
##        -23.25415,-17.69722,-0.7495294, 4*0.,
##        -4.09689, -13.67046, 0.0,       51.9932, 9*0.,
##         0.     ,  0.      , 0.     ,  -26.6083, 9*0./

vak135<-function(r){
   ## AK135 in the outer core.
   if(!exists('ak135.interp')) {
      dat <- matrix(c(
       # R/km     Vp        Vs     ignore
       2891.50,  9.9145,  8.0000,  0.0000,
       2939.33,  9.9942,  8.0382,  0.0000,

       2989.66, 10.0722,  8.1283,  0.0000,
       3039.99, 10.1485,  8.2213,  0.0000,
       3090.32, 10.2233,  8.3122,  0.0000,
       3140.66, 10.2964,  8.4001,  0.0000,
       3190.99, 10.3679,  8.4861,  0.0000,
       3241.32, 10.4378,  8.5692,  0.0000,
       3291.65, 10.5062,  8.6496,  0.0000,
       3341.98, 10.5731,  8.7283,  0.0000,
       3392.31, 10.6385,  8.8036,  0.0000,
       3442.64, 10.7023,  8.8761,  0.0000,

       3492.97, 10.7647,  8.9461,  0.0000,
       3543.30, 10.8257,  9.0138,  0.0000,
       3593.64, 10.8852,  9.0792,  0.0000,
       3643.97, 10.9434,  9.1426,  0.0000,
       3694.30, 11.0001,  9.2042,  0.0000,
       3744.63, 11.0555,  9.2634,  0.0000,
       3794.96, 11.1095,  9.3205,  0.0000,
       3845.29, 11.1623,  9.3760,  0.0000,
       3895.62, 11.2137,  9.4297,  0.0000,
       3945.95, 11.2639,  9.4814,  0.0000,

       3996.28, 11.3127,  9.5306,  0.0000,
       4046.62, 11.3604,  9.5777,  0.0000,
       4096.95, 11.4069,  9.6232,  0.0000,
       4147.28, 11.4521,  9.6673,  0.0000,
       4197.61, 11.4962,  9.7100,  0.0000,
       4247.94, 11.5391,  9.7513,  0.0000,
       4298.27, 11.5809,  9.7914,  0.0000,
       4348.60, 11.6216,  9.8304,  0.0000,
       4398.93, 11.6612,  9.8682,  0.0000,
       4449.26, 11.6998,  9.9051,  0.0000,

       4499.60, 11.7373,  9.9410,  0.0000,
       4549.93, 11.7737,  9.9761,  0.0000,
       4600.26, 11.8092, 10.0103,  0.0000,
       4650.59, 11.8437, 10.0439,  0.0000,
       4700.92, 11.8772, 10.0768,  0.0000,
       4751.25, 11.9098, 10.1095,  0.0000,
       4801.58, 11.9414, 10.1415,  0.0000,
       4851.91, 11.9722, 10.1739,  0.0000,
       4902.24, 12.0001, 10.2049,  0.0000,
       4952.58, 12.0311, 10.2329,  0.0000,

       5002.91, 12.0593, 10.2565,  0.0000,
       5053.24, 12.0867, 10.2745,  0.0000,
       5103.57, 12.1133, 10.2854,  0.0000,
       5153.50, 12.1391, 10.2890,  0.0000 
      ),
      ncol=4, byrow=TRUE)
      ak135.interp <<- splinefun(6371-dat[,1],dat[,3],method="natural")
      ak135.cmbicb <<- 6371-dat[c(1,nrow(dat)),1]
   }
   ak135.interp(r)
}

tad<-function(h,rc=3480,gam=1.52){
   ## Adiabatic temperature profile calculation
   ##    h is depth (km) below reference point (rc)
   ##    Integrates gamma*g(r)/V(r)^2
   ##    Returns multiplier to temperature at reference level (rc),
   ##    e.g. 4300*tad(z) is T at z km below rc if T(rc) = 4300 K.
   fac <- sapply(h,
      function(x)integrate(function(z){g(rc-z)/vprem(rc-z)^2*1e-3},
      0,x)$value
   )
   exp(gam * fac)
}

## Make comparison plot of Yokoyama model predicted and calculated wavespeed
##   data.
## G. Helffrich 16 Oct. 2015

if (!exists('vyoko')) source('yokomodel.R')

dat <- list(
   ##   obs.   off  calc.
   Ti=c(5230, -1.5, 4960),
   V= c(5150, -1.5, 4959),
   Cr=c(4520, -1.5, 4864),
   Mn=c(3710, -1.1, 4044),
   Fe=c(4375, -1.8, 4307),
   Co=c(4060, -1.1, 4172),
   Ni=c(4029, -1.1, 4122)
)

plot(NA, NA,
   xlab=expression(z[Gamma]~(effective~~valence)),
   ylab=expression(c~(m~s^{-1})),
   xlim=c(1.2,1.5), ylim=c(3500,5500)
)

for(n in names(dat)){
   Teta <- uniroot(
      function(T){0.463-vyoko(0, T, spdat[[n]]$rho, sp=n)$eta},
      c(0.5,2)*spdat[[n]]$Tm
   )$root
   stuff <- vyoko(0, spdat[[n]]$Tm, spdat[[n]]$rho, sp=n)
   points(spdat[[n]]$zm, dat[[n]][3], pch=19)
   points(spdat[[n]]$zm, dat[[n]][1], pch=21, col='black')
   points(
      spdat[[n]]$zm, stuff$c,
      pch=22, col='black', cex=1.5
   )
   text(spdat[[n]]$zm, dat[[n]][3], n, adj=c(0.5,dat[[n]][2]))
   with(stuff, cat(paste(n,'c, Z, eta:', c, Z, eta, spdat[[n]]$Tm),'\n'))
}

points(1.35, 4000, pch=19)
   text(1.355, 4000, 'c(calc) Yokoyama (2001) Table 1', adj=c(0,0.5))
points(1.35, 3800, pch=22, col='black')
   text(1.355, 3800, 'Calculated', adj=c(0,0.5))
points(1.35, 3900, pch=21, col='black')
   text(1.355, 3900, 'c(meas) Yokoyama (2001) Table 1', adj=c(0,0.5))

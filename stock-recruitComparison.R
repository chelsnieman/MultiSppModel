rm(list=ls())
# working on setting appropriate parms for recruitment models

#Ricker model

N=seq(0,5000, length.out = 100)

r.a=20
r.b=0.001

R.r=r.a*N*exp(-r.b*N)
plot(N,R.r,pch=16)


# Beverton Holt model

N=seq(0,5000,length.out = 100)

bh.a=5000
bh.b=5e2

R.bh=bh.a*N/(bh.b+N)
plot(N,R.bh,pch=16)


#plotting both together
plot(N,R.bh,type="l",lwd=2, ylab = "Recruits", xlab = "Spawners", ylim = c(0,max(c(R.bh,R.r))))
lines(N,R.r, lwd=2, col="grey")
legend("right",legend = c("Bev-Holt", "Ricker"), lwd=2, col=c("black","grey"), bty = "n")

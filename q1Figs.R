#script to hold code that we're messing around with to make figures before adding it to the RMD doc. 

# trying to make some heat maps to explore how stocking can maintain sp1 dominance
library(deSolve)
library(ggplot2)
library(ggpubr)
rm(list=ls())
source('q2Func.R')


#### BASE MODEL BEHAVIOR, FIG 1 ####

#THIS IS THE CODE FROM THE Q2 FIG, FIND A NEW WAY TO SHOW THIS
store=data.frame(qEs=seq(0,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(100,10,0,0)
tstep=1:300
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(1.8,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(s1=0.1,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.002,v1=1,f1=2,
      s2=0.1,cJ2A2=0.001,cJ2A1=0.5,cJ2J1=0.002,v2=1,f2=2)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
store2=data.frame(qEs=seq(0,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(10,100,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(1.8,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(s1=0.1,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.002,v1=1,f1=2,
      s2=0.1,cJ2A2=0.001,cJ2A1=0.5,cJ2J1=0.002,v2=1,f2=2)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}

par(mfcol=c(2,1), mar=c(1,1,3,1), oma=c(4,4,1,1))
plot(store$qEs,store$A1,lwd=3,type='l',ylim=c(0,max(store[,2:3],na.rm = T)),ylab = "",xlab = "", main = "Sp1 > Sp2")
lines(store$qEs,store$A2,lwd=3,col='grey')
legend("topright",legend = c("sp 1", "sp 2"), lty=1, lwd=2, col = c("black","grey"),bty="n")
plot(store2$qEs,store2$A1,lwd=3, type='l', ylim = c(0, max(store2[,2:3],na.rm = T)),ylab = "", xlab = "", main = "Sp2 > Sp1")
lines(store2$qEs,store2$A2,lwd=3,col='grey')
legend("topright",legend = c("sp 1", "sp 2"), lty=1, lwd=2, col = c("black","grey"),bty="n")
mtext("Abundance", side = 2, outer = T, line = 2.5)
mtext("Species 1 Harvest Rate (q*E)", side = 1, outer = T, line = 2.5)

#trying to find a new way to show base model behavior for q1

y0=c(100,10,0,0)
tstep=1:300
qE1Fun=approxfun(x=tstep,y=c(rep(1.8,length(tstep))))
qE2Fun=approxfun(x=tstep,y=rep(1.8,length(tstep)))
h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))

p=c(s1=0.1,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.002,v1=1,f1=2,
    s2=0.1,cJ2A2=0.001,cJ2A1=0.5,cJ2J1=0.002,v2=1,f2=2)
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

plot(tstep,  sim[,2], type='l')
lines(tstep, sim[,5]/sim[,2], col='green')

plot(sim[,2], sim[,4]/sim[,2])
plot(sim[,2], sim[,5]/sim[,2])
points(sim[,3], sim[,4]/sim[,3],col='red')

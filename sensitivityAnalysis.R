## 12.30.2020
## Sensitivity analysis for model parms; based off 'BiggsModParmDiffs.R' script

rm(list=ls())
library(deSolve)
library(ggplot2)
library(ggpubr)
#model function
# qE - harvest, species specific, this is what is controlled by 'regulations'
# s - juvenile overwinter survival, species specific
# m - adult natural mortality rate, species specific
# cJA - effect of adults of a given species on juveniles of a given species (cover cannibalism or interspecific predation, both happen in foraging arena)
# cJJ - effect of juveniles of one species on juveniles of the other (can be predation or competition)
# h - rate at which juveniles leave foraging arena for refuge, species specific
# v - rate at which juveniles enter foraging arena from refuge, species specific
# stock1 - annual stocked num spp1
# stock2 - annual stocked num spp2
simBiggsQ2<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(5000*A1/(500+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

#### SURVIVAL LOOP ####
times=1:300

ss=seq(0.05,.9,length.out = 10)
store=data.frame(qEs=rep(seq(.05,8,length.out=30),10),A1=0,A2=0,J1=0,J2=0,ss=rep(0,300))
for(f in 1:length(ss)){
  
  y0=c(10,100,0,0)
  for(i in 1:nrow(store)){
    p=c(c(qE1=store$qEs[i],s1=ss$ss[f],cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=1,f1=1),
        c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=20,f2=1))
    sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
    store$A1[i]=sim[nrow(sim),2]
    store$A2[i]=sim[nrow(sim),3]
    store$J1[i]=sim[nrow(sim),4]
    store$J2[i]=sim[nrow(sim),5]
  }
  store2=data.frame(qEs=seq(.05,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
  y0=c(100,10,40,40)
  for(i in 1:nrow(store)){
    p=c(c(qE1=store2$qEs[i],s1=ss$ss[f],cJ1A1=0.002,cJ1A2=0.002,cJ1J2=0.001,v1=1,h1=6,f1=1),
        c(qE2=1.8,s2=0.5,cJ2A2=0.002,cJ2A1=0.002,cJ2J1=0.001,v2=1,h2=6,f2=1))
    sim=ode(y=y0,times=times,func=simBiggs3,parms=p)
    store2$A1[i]=sim[nrow(sim),2]
    store2$A2[i]=sim[nrow(sim),3]
    store2$J1[i]=sim[nrow(sim),4]
    store2$J2[i]=sim[nrow(sim),5]
  }
  if(any(abs(store$A2-store2$A2)>1) & any(abs(store$A1-store2$A1)>1)){ss$hyst[f]=1}
  
}
plot(ss$ss,ss$hyst,pch=16,xlab = "Survival",ylab = "Hysteresis Present?")
plot(ss$ss-0.5,ss$hyst,pch=16,xlab = "Diff Survival",ylab = "Hysteresis Present?")

## CD 12.22.2020
## Supplementary information figures
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

#single model run, describe these dynamics in the beginning of the results.This fig may go in supplemental
# slow increase in harvest of species 1 brings its abund down, system flips to sp2 over time. This is in a system where all else is the same. Stochasticity and other compeititive imbalances speed this flip up.
tstep=1:300
qE1Fun=approxfun(x=tstep,y=c(seq(0,25,length.out=length(tstep))))
qE2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
y0=c(100,10,0,0)

p=c(s1=0.1,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1)
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
sf1=as.data.frame(rbind(cbind(sim[,c(1,2)],rep("A1",nrow(sim))),cbind(sim[,c(1,3)],rep("A2",nrow(sim)))))
colnames(sf1)=c("Time","Abund","sp")
sf1$Abund=as.numeric(sf1$Abund);sf1$Time=as.numeric(sf1$Time)

figS1=ggplot(data = sf1,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  labs(x='Time',y='Adult Abundance')+
  theme(legend.position = 'bottom')
figS1


#same thing to demonstrate the effect of refuge loss, describe in methods, fig may go in supplement
# no harvest so the system doesn't flip but the decline in habitat brings down sp1. Since juves share same habitat, the decline effects them both. If you add even a little harvest to sp1 here it flips immediately.
tstep=1:300
qE1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
qE2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
h1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
h2Fun=approxfun(x=tstep,y=c(rep(20,length(tstep))))
st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
y0=c(100,10,0,0)

p=c(s1=0.1,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1)
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
sf2=as.data.frame(rbind(cbind(sim[,c(1,2)],rep("A1",nrow(sim))),cbind(sim[,c(1,3)],rep("A2",nrow(sim)))))
colnames(sf2)=c("Time","Abund","sp")
sf2$Abund=as.numeric(sf2$Abund);sf2$Time=as.numeric(sf2$Time)

figS2=ggplot(data = sf2,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  labs(x='Time',y='Adult Abundance')+
  theme(legend.position = 'bottom')
figS2


#flip scenarios

tstep=1:300

h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
minDiff=100

#### FLIP W/O sp2 harv ####

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
df2=expand.grid(X=qEs,Y=sto)
df2$A1=numeric(nrow(df2))
df2$A2=numeric(nrow(df2))
df2$J1=numeric(nrow(df2))
df2$J2=numeric(nrow(df2))


for(i in 1:nrow(df2)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(df2$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(df2$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df2$A1[i]=sim[nrow(sim)-1,2]
  df2$A2[i]=sim[nrow(sim)-1,3]
  df2$J1[i]=sim[nrow(sim)-1,4]
  df2$J2[i]=sim[nrow(sim)-1,5]
}

df2$diff=df2$A1-df2$A2
df2$outcome=ifelse(df2$A1-df2$A2 < minDiff,"darkgreen", "darkred")


#### FLIP W/ sp2 harv ####

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
dfwo2=expand.grid(X=qEs,Y=sto)
dfwo2$A1=numeric(nrow(dfwo2))
dfwo2$A2=numeric(nrow(dfwo2))
dfwo2$J1=numeric(nrow(dfwo2))
dfwo2$J2=numeric(nrow(dfwo2))

for(i in 1:nrow(dfwo2)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfwo2$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(2, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfwo2$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo2$A1[i]=sim[nrow(sim)-1,2]
  dfwo2$A2[i]=sim[nrow(sim)-1,3]
  dfwo2$J1[i]=sim[nrow(sim)-1,4]
  dfwo2$J2[i]=sim[nrow(sim)-1,5]
}

dfwo2$diff=dfwo2$A1-dfwo2$A2
dfwo2$outcome=ifelse(dfwo2$A1-dfwo2$A2 < minDiff,"darkgreen", "darkred")
df2$mod=rep("Flip to 1, ignore 2",nrow(df2))
dfwo2$mod=rep("Flip to 1, harv 2",nrow(dfwo2))
allSen=rbind(df2,dfwo2)
vzI=ggplot(data=allSen, aes(x=allSen$X,y=allSen$Y,linetype=allSen$mod))+theme_classic()+
  geom_contour(aes(z=allSen$diff),breaks = c(minDiff), color='black', size=2)+
  labs(x="Species 1 Harvest Rate", y="Species 1 Stocking",linetype="Scenario")+
  theme(legend.position = 'bottom')
vzI

#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. Flip scenario instead of maintain scenario 

minDiff=100
qEs=rep(seq(0,8,length.out = 30),3)
sto=rep(seq(0,2000, length.out = 30),3)
dfT=expand.grid(X=qEs,Y=sto)
dfT$A1=numeric(nrow(dfT))
dfT$A2=numeric(nrow(dfT))
dfT$J1=numeric(nrow(dfT))
dfT$J2=numeric(nrow(dfT))
dfT$sp1H=c(rep(2,30),rep(4,30),rep(6,30))

for(i in 1:nrow(dfT)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfT$sp1H[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(dfT$X[i], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfT$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfT$A1[i]=sim[nrow(sim)-1,2]
  dfT$A2[i]=sim[nrow(sim)-1,3]
  dfT$J1[i]=sim[nrow(sim)-1,4]
  dfT$J2[i]=sim[nrow(sim)-1,5]
}

dfT$diff=dfT$A1-dfT$A2
dfT$qNorm=(dfT$X-min(dfT$X))/(max(dfT$X)-min(dfT$X))
dfT$sNorm=(dfT$Y-min(dfT$Y))/(max(dfT$Y)-min(dfT$Y))
dfT$sp1Norm=((dfT$sp1H-min(dfT$X))/(max(dfT$X)-min(dfT$X))) #putting sp1 harv on same scale as sp2

vzT=ggplot(data=dfT, aes(x=dfT$X,y=dfT$Y,linetype=as.factor(dfT$sp1H)))+theme_classic()+
  geom_contour(aes(z=dfT$diff),breaks = c(minDiff), color='black', size=1)+
  labs(x="Species 2 Harvest Rate", y="Species 1 Stocking",linetype="Species 1 Harvest")+
  theme(legend.position = 'bottom')+
  xlim(0,8)+
  ylim(0,2000)
vzT
## code to think about new ways to vizualise some of our results

rm(list=ls())
library(deSolve)

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
plot(sim[,1],sim[,2],type = 'l',lwd=2,ylim = c(0,max(sim[,2:3],na.rm = T)))
lines(sim[,1],sim[,3],col='gray',lwd=2)


#same thing to demonstrate the effect of refuge loss, describe in methods, fig may go in supplement
# no harvest so the system doesn't flip but the decline in habitat brings down sp1. Since juves share same habitat, the decline effects them both. If you add even a little harves to sp1 here it flips immediately.
tstep=1:300
qE1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
qE2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
h1Fun=approxfun(x=tstep,y=c(seq(8,0,length.out=100),rep(0,200)))
h2Fun=h1Fun
st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
y0=c(100,10,0,0)

p=c(s1=0.1,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1)
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
plot(sim[,1],sim[,2],type = 'l',lwd=2,ylim = c(0,max(sim[,2:3],na.rm = T)))
lines(sim[,1],sim[,3],col='gray',lwd=2)


#demonstrate alternative stable states
store=data.frame(qEs=seq(0,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(100,10,0,0)
tstep=1:300
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(s1=0.1,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1)
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
  
  p=c(s1=0.1,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}

par(mfcol=c(2,1), mar=c(1,1,3,1), oma=c(4,4,1,1))
plot(store$qEs,store$A1,lwd=3,type='l',ylim=c(0,max(store[,2:3],na.rm = T)),ylab = "",xlab = "")
lines(store$qEs,store$A2,lwd=3,col='grey')
legend("topright",legend = c("sp 1", "sp 2"), lty=1, lwd=2, col = c("black","grey"),bty="n")
plot(store2$qEs,store2$A1,lwd=3, type='l', ylim = c(0, max(store2[,2:3],na.rm = T)),ylab = "", xlab = "")
lines(store2$qEs,store2$A2,lwd=3,col='grey')
#legend("topright",legend = c("sp 1", "sp 2"), lty=1, lwd=2, col = c("black","grey"),bty="n")
mtext("Abundance", side = 2, outer = T, line = 2.5)
mtext("Species 1 Harvest Rate (q*E)", side = 1, outer = T, line = 2.5)

##### HEATMAPS #####

tstep=1:300

h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))

##### MAINTAIN W/O sp2 harv ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
df=expand.grid(X=qEs,Y=sto)
df$A1=numeric(nrow(df))
df$A2=numeric(nrow(df))
df$J1=numeric(nrow(df))
df$J2=numeric(nrow(df))
minDiff=100

for(i in 1:nrow(df)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(df$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(df$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(1000,100,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df$A1[i]=sim[nrow(sim)-1,2]
  df$A2[i]=sim[nrow(sim)-1,3]
  df$J1[i]=sim[nrow(sim)-1,4]
  df$J2[i]=sim[nrow(sim)-1,5]
}

df$diff=df$A1-df$A2
df$outcome=ifelse(df$A1-df$A2 < minDiff,"darkgreen", "darkred")


#### FLIP W/O sp2 harv ####

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
df2=expand.grid(X=qEs,Y=sto)
df2$A1=numeric(nrow(df2))
df2$A2=numeric(nrow(df2))
df2$J1=numeric(nrow(df2))
df2$J2=numeric(nrow(df2))


for(i in 1:nrow(df)){
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


##### MAINTAIN W/ sp2 harv ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
dfwo=expand.grid(X=qEs,Y=sto)
dfwo$A1=numeric(nrow(dfwo))
dfwo$A2=numeric(nrow(dfwo))
dfwo$J1=numeric(nrow(dfwo))
dfwo$J2=numeric(nrow(dfwo))

for(i in 1:nrow(dfwo)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfwo$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(2, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfwo$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(1000,100,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo$A1[i]=sim[nrow(sim)-1,2]
  dfwo$A2[i]=sim[nrow(sim)-1,3]
  dfwo$J1[i]=sim[nrow(sim)-1,4]
  dfwo$J2[i]=sim[nrow(sim)-1,5]
}

dfwo$diff=dfwo$A1-dfwo$A2
dfwo$outcome=ifelse(dfwo$A1-dfwo$A2 < minDiff,"darkgreen", "darkred")


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


#trying some new plots

#plot isoclines for all scenarios on one plot with x=harvest, y= stocking, line where 1>2
df$mod=rep("Maintain 1, ignore 2",nrow(df))
#df2$mod=rep("Flip to 1, ignore 2",nrow(df2))
dfwo$mod=rep("Maintain 1, harv 2",nrow(dfwo))
#dfwo2$mod=rep("Flip to 1, harv 2",nrow(dfwo2))
allSen=rbind(df,dfwo)
vzI=ggplot(data=allSen, aes(x=allSen$X,y=allSen$Y,linetype=allSen$mod))+theme_classic()+
  geom_contour(aes(z=allSen$diff),breaks = c(minDiff), color='black', size=2)+
  labs(x="Species 1 Harvest Rate", y="Species 1 Stocking",linetype="Scenario")+
  theme(legend.position = 'bottom')
vzI
#comparing with and without hysteresis plots together 
#dev.new()
#ggarrange(a,b,c,d,labels = c("A","B","C","D"), nrow = 2,ncol = 2, common.legend = T, legend = 'bottom')


#new figure 3, isoclines at different harvests

#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. 

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
  y0=c(1000,100,0,0)
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

vzT=ggplot(data=dfT, aes(x=dfT$qNorm,y=dfT$sNorm,linetype=as.factor(dfT$sp1Norm)))+theme_classic()+
  geom_contour(aes(z=dfT$diff),breaks = c(minDiff), color='black', size=1)+
  labs(x="Species 2 Harvest Rate", y="Species 1 Stocking",linetype="Species 1 Harvest")+
  theme(legend.position = 'bottom')+
  xlim(0,1)+
  ylim(0,1)
vzT

#flipping sp1 and sp2 harvest in the plot requires rerunning loop
minDiff=100
qEs=rep(seq(0,8,length.out = 30),3)
sto=rep(seq(0,2000, length.out = 30),3)
dfT2=expand.grid(X=qEs,Y=sto)
dfT2$A1=numeric(nrow(dfT2))
dfT2$A2=numeric(nrow(dfT2))
dfT2$J1=numeric(nrow(dfT2))
dfT2$J2=numeric(nrow(dfT2))
dfT2$sp2H=c(rep(2,30),rep(4,30),rep(6,30))

for(i in 1:nrow(dfT2)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfT2$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(dfT2$sp2H[i], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfT2$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(1000,100,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfT2$A1[i]=sim[nrow(sim)-1,2]
  dfT2$A2[i]=sim[nrow(sim)-1,3]
  dfT2$J1[i]=sim[nrow(sim)-1,4]
  dfT2$J2[i]=sim[nrow(sim)-1,5]
}

dfT2$diff=dfT2$A1-dfT2$A2
dfT2$qNorm=(dfT2$X-min(dfT2$X))/(max(dfT2$X)-min(dfT2$X))
dfT2$sNorm=(dfT2$Y-min(dfT2$Y))/(max(dfT2$Y)-min(dfT2$Y))
dfT2$sp2Norm=((dfT2$sp2H-min(dfT2$X))/(max(dfT2$X)-min(dfT2$X))) #putting sp1 harv on same scale as sp2

vzT2=ggplot(data=dfT2, aes(x=dfT2$qNorm,y=dfT2$sNorm,linetype=as.factor(dfT2$sp2Norm)))+theme_classic()+
  geom_contour(aes(z=dfT2$diff),breaks = c(minDiff), color='black', size=1)+
  labs(x="Species 1 Harvest Rate", y="Species 1 Stocking",linetype="Species 2 Harvest")+
  theme(legend.position = 'bottom')+
  xlim(0,1)+
  ylim(0,1)
vzT2


# similar plot where we try to maintain sp2 instead of sp1

minDiff=100
qEs=rep(seq(0,8,length.out = 30),3)
sto=rep(seq(0,2000, length.out = 30),3)
dfT3=expand.grid(X=qEs,Y=sto)
dfT3$A1=numeric(nrow(dfT3))
dfT3$A2=numeric(nrow(dfT3))
dfT3$J1=numeric(nrow(dfT3))
dfT3$J2=numeric(nrow(dfT3))
dfT3$sp1H=c(rep(2,30),rep(4,30),rep(6,30))

for(i in 1:nrow(dfT3)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfT3$sp1H[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(dfT3$X[i], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfT3$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfT3$A1[i]=sim[nrow(sim)-1,2]
  dfT3$A2[i]=sim[nrow(sim)-1,3]
  dfT3$J1[i]=sim[nrow(sim)-1,4]
  dfT3$J2[i]=sim[nrow(sim)-1,5]
}

dfT3$diff=dfT3$A2-dfT3$A1
dfT3$qNorm=(dfT3$X-min(dfT3$X))/(max(dfT3$X)-min(dfT3$X))
dfT3$sNorm=(dfT3$Y-min(dfT3$Y))/(max(dfT3$Y)-min(dfT3$Y))
dfT3$sp1Norm=((dfT3$sp1H-min(dfT3$X))/(max(dfT3$X)-min(dfT3$X))) #putting sp1 harv on same scale as sp2

vzT3=ggplot(data=dfT3, aes(x=dfT3$qNorm,y=dfT3$sNorm,linetype=as.factor(dfT3$sp1Norm)))+theme_classic()+
  geom_contour(aes(z=dfT3$diff),breaks = c(minDiff), color='black', size=1)+
  labs(x="Species 2 Harvest Rate", y="Species 1 Stocking",linetype="Species 1 Harvest")+
  theme(legend.position = 'bottom')+
  xlim(0,1)+
  ylim(0,1)
vzT3

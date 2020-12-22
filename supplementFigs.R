## CD 12.22.2020
## Supplementary information figures


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

allSen=rbind(df2,dfwo2)
vzI=ggplot(data=allSen, aes(x=allSen$X,y=allSen$Y,linetype=allSen$mod))+theme_classic()+
  geom_contour(aes(z=allSen$diff),breaks = c(minDiff), color='black', size=2)+
  labs(x="Species 1 Harvest Rate", y="Species 1 Stocking",linetype="Scenario")+
  theme(legend.position = 'bottom')
vzI
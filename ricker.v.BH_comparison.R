#comparing different recruitment models in the q2 function to see what the effects are, if any.
rm(list=ls())
library(deSolve)
rickerQ2<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(20*A1*exp(-0.001*A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(20*A2*exp(-0.001*A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}
BHQ2<-function(t,y,params){
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


#### Q2 DRAFT FIG. 1 - RICKER ####
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
  sim=ode(y=y0,times=tstep,func=rickerQ2,parms=p)
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
  sim=ode(y=y0,times=tstep,func=rickerQ2,parms=p)
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
title(main = "Ricker")

#### Q2 DRAFT FIG. 1 - BH ####
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
  sim=ode(y=y0,times=tstep,func=BHQ2,parms=p)
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
  sim=ode(y=y0,times=tstep,func=BHQ2,parms=p)
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
title(main = "BH")

#### Q2 DRAFT FIG.2 - RICKER ####

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
  sim=ode(y=y0,times=tstep,func=rickerQ2,parms=p)
  df$A1[i]=sim[nrow(sim)-1,2]
  df$A2[i]=sim[nrow(sim)-1,3]
  df$J1[i]=sim[nrow(sim)-1,4]
  df$J2[i]=sim[nrow(sim)-1,5]
}

df$diff=df$A1-df$A2
df$outcome=ifelse(df$A1-df$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = df)+theme_classic()
a=vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(df$diff), 1, max(df$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "")+
  geom_contour(aes(x=df$X, y=df$Y, z=df$diff), breaks = c(minDiff))

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
  sim=ode(y=y0,times=tstep,func=rickerQ2,parms=p)
  df2$A1[i]=sim[nrow(sim)-1,2]
  df2$A2[i]=sim[nrow(sim)-1,3]
  df2$J1[i]=sim[nrow(sim)-1,4]
  df2$J2[i]=sim[nrow(sim)-1,5]
}

df2$diff=df2$A1-df2$A2
df2$outcome=ifelse(df2$A1-df2$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = df2)+theme_classic()
b=vz+geom_tile(aes(x=df2$X, y=df2$Y, fill=df2$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(df2$diff), 1, max(df2$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "")+
  geom_contour(aes(x=df2$X, y=df2$Y, z=df2$diff), breaks = c(minDiff))

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
  sim=ode(y=y0,times=tstep,func=rickerQ2,parms=p)
  dfwo$A1[i]=sim[nrow(sim)-1,2]
  dfwo$A2[i]=sim[nrow(sim)-1,3]
  dfwo$J1[i]=sim[nrow(sim)-1,4]
  dfwo$J2[i]=sim[nrow(sim)-1,5]
}

dfwo$diff=dfwo$A1-dfwo$A2
dfwo$outcome=ifelse(dfwo$A1-dfwo$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfwo)+theme_classic()
c=vz+geom_tile(aes(x=dfwo$X, y=dfwo$Y, fill=dfwo$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(dfwo$diff), 1, max(dfwo$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "")+
  geom_contour(aes(x=dfwo$X, y=dfwo$Y, z=dfwo$diff), breaks = c(minDiff))

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
  sim=ode(y=y0,times=tstep,func=rickerQ2,parms=p)
  dfwo2$A1[i]=sim[nrow(sim)-1,2]
  dfwo2$A2[i]=sim[nrow(sim)-1,3]
  dfwo2$J1[i]=sim[nrow(sim)-1,4]
  dfwo2$J2[i]=sim[nrow(sim)-1,5]
}

dfwo2$diff=dfwo2$A1-dfwo2$A2
dfwo2$outcome=ifelse(dfwo2$A1-dfwo2$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfwo2)+theme_classic()
d=vz+geom_tile(aes(x=dfwo2$X, y=dfwo2$Y, fill=dfwo2$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(dfwo2$diff), 1, max(dfwo2$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "")+
  geom_contour(aes(x=dfwo2$X, y=dfwo2$Y, z=dfwo2$diff), breaks = c(minDiff))

#comparing with and without hysteresis plots together using cowplot package
#dev.new()
ggarrange(a,b,c,d,labels = c("A","B","C","D"), nrow = 2,ncol = 2, common.legend = T, legend = 'bottom')

#### Q2 DRAFT FIG.2 - BH ####

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
  sim=ode(y=y0,times=tstep,func=BHQ2,parms=p)
  df$A1[i]=sim[nrow(sim)-1,2]
  df$A2[i]=sim[nrow(sim)-1,3]
  df$J1[i]=sim[nrow(sim)-1,4]
  df$J2[i]=sim[nrow(sim)-1,5]
}

df$diff=df$A1-df$A2
df$outcome=ifelse(df$A1-df$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = df)+theme_classic()
a=vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(df$diff), 1, max(df$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "")+
  geom_contour(aes(x=df$X, y=df$Y, z=df$diff), breaks = c(minDiff))

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
  sim=ode(y=y0,times=tstep,func=BHQ2,parms=p)
  df2$A1[i]=sim[nrow(sim)-1,2]
  df2$A2[i]=sim[nrow(sim)-1,3]
  df2$J1[i]=sim[nrow(sim)-1,4]
  df2$J2[i]=sim[nrow(sim)-1,5]
}

df2$diff=df2$A1-df2$A2
df2$outcome=ifelse(df2$A1-df2$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = df2)+theme_classic()
b=vz+geom_tile(aes(x=df2$X, y=df2$Y, fill=df2$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(df2$diff), 1, max(df2$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "")+
  geom_contour(aes(x=df2$X, y=df2$Y, z=df2$diff), breaks = c(minDiff))

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
  sim=ode(y=y0,times=tstep,func=BHQ2,parms=p)
  dfwo$A1[i]=sim[nrow(sim)-1,2]
  dfwo$A2[i]=sim[nrow(sim)-1,3]
  dfwo$J1[i]=sim[nrow(sim)-1,4]
  dfwo$J2[i]=sim[nrow(sim)-1,5]
}

dfwo$diff=dfwo$A1-dfwo$A2
dfwo$outcome=ifelse(dfwo$A1-dfwo$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfwo)+theme_classic()
c=vz+geom_tile(aes(x=dfwo$X, y=dfwo$Y, fill=dfwo$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(dfwo$diff), 1, max(dfwo$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "")+
  geom_contour(aes(x=dfwo$X, y=dfwo$Y, z=dfwo$diff), breaks = c(minDiff))

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
  sim=ode(y=y0,times=tstep,func=BHQ2,parms=p)
  dfwo2$A1[i]=sim[nrow(sim)-1,2]
  dfwo2$A2[i]=sim[nrow(sim)-1,3]
  dfwo2$J1[i]=sim[nrow(sim)-1,4]
  dfwo2$J2[i]=sim[nrow(sim)-1,5]
}

dfwo2$diff=dfwo2$A1-dfwo2$A2
dfwo2$outcome=ifelse(dfwo2$A1-dfwo2$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfwo2)+theme_classic()
d=vz+geom_tile(aes(x=dfwo2$X, y=dfwo2$Y, fill=dfwo2$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(dfwo2$diff), 1, max(dfwo2$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "")+
  geom_contour(aes(x=dfwo2$X, y=dfwo2$Y, z=dfwo2$diff), breaks = c(minDiff))

#comparing with and without hysteresis plots together using cowplot package
#dev.new()
ggarrange(a,b,c,d,labels = c("A","B","C","D"), nrow = 2,ncol = 2, common.legend = T, legend = 'bottom')

#### Q2 DRAFT FIG. 3 - RICKER ####
#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. 

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
dfT=expand.grid(X=qEs,Y=sto)
dfT$A1=numeric(nrow(dfT))
dfT$A2=numeric(nrow(dfT))
dfT$J1=numeric(nrow(dfT))
dfT$J2=numeric(nrow(dfT))

for(i in 1:nrow(dfT)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(4, length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(dfT$X[i], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfT$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=rickerQ2,parms=p)
  dfT$A1[i]=sim[nrow(sim)-1,2]
  dfT$A2[i]=sim[nrow(sim)-1,3]
  dfT$J1[i]=sim[nrow(sim)-1,4]
  dfT$J2[i]=sim[nrow(sim)-1,5]
}

dfT$diff=dfT$A1-dfT$A2
dfT$qNorm=(dfT$X-min(dfT$X))/(max(dfT$X)-min(dfT$X))
dfT$sNorm=(dfT$Y-min(dfT$Y))/(max(dfT$Y)-min(dfT$Y))
#dfT$outcome=ifelse(dfT$A1-dfT$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfT)+theme_classic()
vz+geom_tile(aes(x=dfT$qNorm, y=dfT$sNorm, fill=dfT$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 \nDominance", breaks=c(min(dfT$diff), 1, max(dfT$diff)), labels=c("sp2", "even", "sp1"))+
  labs(x="Sp2 harvest Rate", y="Sp1 stocked fish")+
  ggtitle("Ricker")
#geom_contour(aes(x=dfT$X, y=dfT$Y, z=dfT$diff))+
#geom_label_contour(aes(x=dfT$X, y=dfT$Y, z=dfT$diff))

#### Q2 DRAFT FIG. 3 - BH ####
#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. 

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
dfT=expand.grid(X=qEs,Y=sto)
dfT$A1=numeric(nrow(dfT))
dfT$A2=numeric(nrow(dfT))
dfT$J1=numeric(nrow(dfT))
dfT$J2=numeric(nrow(dfT))

for(i in 1:nrow(dfT)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(4, length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(dfT$X[i], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfT$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=BHQ2,parms=p)
  dfT$A1[i]=sim[nrow(sim)-1,2]
  dfT$A2[i]=sim[nrow(sim)-1,3]
  dfT$J1[i]=sim[nrow(sim)-1,4]
  dfT$J2[i]=sim[nrow(sim)-1,5]
}

dfT$diff=dfT$A1-dfT$A2
dfT$qNorm=(dfT$X-min(dfT$X))/(max(dfT$X)-min(dfT$X))
dfT$sNorm=(dfT$Y-min(dfT$Y))/(max(dfT$Y)-min(dfT$Y))
#dfT$outcome=ifelse(dfT$A1-dfT$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfT)+theme_classic()
vz+geom_tile(aes(x=dfT$qNorm, y=dfT$sNorm, fill=dfT$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 \nDominance", breaks=c(min(dfT$diff), 1, max(dfT$diff)), labels=c("sp2", "even", "sp1"))+
  labs(x="Sp2 harvest Rate", y="Sp1 stocked fish")+
  ggtitle("B-H")
#geom_contour(aes(x=dfT$X, y=dfT$Y, z=dfT$diff))+
#geom_label_contour(aes(x=dfT$X, y=dfT$Y, z=dfT$diff))

#### Q2 DRAFT FIG. 4 - RICKER ####
#figure on delaying a transition

#running to sp1 dominating before habitat decline
tstep=1:200
h1Fun=approxfun(x=tstep,y=c(rep(10,200)))
h2Fun=approxfun(x=tstep,y=c(rep(10,200)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
qE2Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
y0=c(500,100,0,0)
simPre=ode(y=y0,times=tstep,func=rickerQ2,parms=p)

#habitat decline and stocking
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0, length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,10, length.out = 100),rep(10,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)))
st1Fun=approxfun(x=tstep,y=c(rep(0,500)))
st2Fun=approxfun(x=tstep,y=c(rep(4,500)))#this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amount of harvest that declines to 0
p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
y0=simPre[199,2:5]
sim=ode(y=y0,times=tstep,func=rickerQ2,parms=p)

#stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0, length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,10, length.out = 100),rep(10,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)))
st1Fun=approxfun(x=tstep,y=c(rep(0,100),rep(80,400)))
st2Fun=approxfun(x=tstep,y=c(rep(4,500))) #this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amount of harvest that declines to 0

simS=ode(y=y0,times=tstep,func=rickerQ2,parms=p)

#harvesting delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0,length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,10,length.out = 100),rep(10,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,.5,length.out = 200), rep(.5,200)))
st1Fun=approxfun(x=tstep,y=c(rep(0,500)))
st2Fun=approxfun(x=tstep,y=c(rep(5,500))) #this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amount of harvest that declines to 0

simH=ode(y=y0,times=tstep,func=rickerQ2,parms=p)

#harvesting and stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0,length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,10,length.out = 100),rep(10,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,.5,length.out = 200), rep(.5,200)))
st1Fun=approxfun(x=tstep,y=c(rep(0,100),rep(80,400)))
st2Fun=approxfun(x=tstep,y=c(rep(5,500))) #this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amount of harvest that declines to 0

simB=ode(y=y0,times=tstep,func=rickerQ2,parms=p)

#plotting
par(mfcol=c(2,2), mar=c(1,1,1,1), oma=c(4,4,1,1))
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim[,1],sim[,3],col='black',lwd=2)
text(x=375,y=1200,labels = "No action taken")
text(x=max(sim[,1],na.rm = T), y=max(sim[,2:3],na.rm = T),labels='A',font = 2)
legend("topleft", legend = c("sp1","sp2"), col = c('grey','black'),lty = c(1,1),lwd=2,bty = "n")
plot(simS[,1],simS[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(simS[,1],simS[,3],col='black',lwd=2)
text(x=375, y=1200, labels = "Stocking sp1")
text(x=min(simS[,1],na.rm=T), y=max(simS[,2:3],na.rm = T),labels='B',font = 2)
plot(simH[,1],simH[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(simH[,1],simH[,3],col='black',lwd=2)
text(x=375, y=1200, labels = "Harvesting sp2")
text(x=min(simH[,1],na.rm = T), y=max(simH[,2:3],na.rm = T),labels='C',font = 2)
plot(simB[,1],simB[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(simB[,1],simB[,3],col='black',lwd=2)
text(x=375, y=1200, labels = "Harvesting sp2 \n Stock sp1")
text(x=min(simB[,1],na.rm = T), y=max(simB[,2:3],na.rm = T),labels="D",font = 2)
mtext("Year", side = 1, outer = T, line = 2.5)
mtext("Population Size", side = 2, outer = T, line = 2.5)

#### Q2 DRAFT FIG. 4 - BH ####
#figure on delaying a transition

#running to sp1 dominating before habitat decline
tstep=1:200
h1Fun=approxfun(x=tstep,y=c(rep(10,200)))
h2Fun=approxfun(x=tstep,y=c(rep(10,200)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
qE2Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
y0=c(500,100,0,0)
simPre=ode(y=y0,times=tstep,func=BHQ2,parms=p)

#habitat decline and stocking
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0, length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,10, length.out = 100),rep(10,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)))
st1Fun=approxfun(x=tstep,y=c(rep(0,500)))
st2Fun=approxfun(x=tstep,y=c(rep(4,500)))#this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amount of harvest that declines to 0
p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
y0=simPre[199,2:5]
sim=ode(y=y0,times=tstep,func=BHQ2,parms=p)

#stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0, length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,10, length.out = 100),rep(10,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)))
st1Fun=approxfun(x=tstep,y=c(rep(0,100),rep(80,400)))
st2Fun=approxfun(x=tstep,y=c(rep(4,500))) #this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amount of harvest that declines to 0

simS=ode(y=y0,times=tstep,func=BHQ2,parms=p)

#harvesting delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0,length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,10,length.out = 100),rep(10,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,.5,length.out = 200), rep(.5,200)))
st1Fun=approxfun(x=tstep,y=c(rep(0,500)))
st2Fun=approxfun(x=tstep,y=c(rep(5,500))) #this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amount of harvest that declines to 0

simH=ode(y=y0,times=tstep,func=BHQ2,parms=p)

#harvesting and stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0,length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,10,length.out = 100),rep(10,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,.5,length.out = 200), rep(.5,200)))
st1Fun=approxfun(x=tstep,y=c(rep(0,100),rep(80,400)))
st2Fun=approxfun(x=tstep,y=c(rep(5,500))) #this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amount of harvest that declines to 0

simB=ode(y=y0,times=tstep,func=BHQ2,parms=p)

#plotting
par(mfcol=c(2,2), mar=c(1,1,1,1), oma=c(4,4,1,1))
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim[,1],sim[,3],col='black',lwd=2)
text(x=375,y=1200,labels = "No action taken")
text(x=max(sim[,1],na.rm = T), y=max(sim[,2:3],na.rm = T),labels='A',font = 2)
legend("topleft", legend = c("sp1","sp2"), col = c('grey','black'),lty = c(1,1),lwd=2,bty = "n")
plot(simS[,1],simS[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(simS[,1],simS[,3],col='black',lwd=2)
text(x=375, y=1200, labels = "Stocking sp1")
text(x=min(simS[,1],na.rm=T), y=max(simS[,2:3],na.rm = T),labels='B',font = 2)
plot(simH[,1],simH[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(simH[,1],simH[,3],col='black',lwd=2)
text(x=375, y=1200, labels = "Harvesting sp2")
text(x=min(simH[,1],na.rm = T), y=max(simH[,2:3],na.rm = T),labels='C',font = 2)
plot(simB[,1],simB[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(simB[,1],simB[,3],col='black',lwd=2)
text(x=375, y=1200, labels = "Harvesting sp2 \n Stock sp1")
text(x=min(simB[,1],na.rm = T), y=max(simB[,2:3],na.rm = T),labels="D",font = 2)
mtext("Year", side = 1, outer = T, line = 2.5)
mtext("Population Size", side = 2, outer = T, line = 2.5)


#### MIXED RICKER & BH RECRUITMENT ####

# Here I'm making sp1 (walleye) rectuitment ricker, and sp2 (bass) BH which might be the closest to reality

mixedQ2<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(20*A1*exp(-0.001*A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

#tuning to make sure roughly the same maximum number of recruits are produced
ds=seq(0,5000,by=500)
lmb=(5000*ds/(500+ds))
wly=(20*ds*exp(-1e-3*ds))
plot(ds, lmb, type = 'l', lwd=2, col='grey', ylim=c(0,max(c(lmb,wly))))
lines(ds, wly, lwd=2, col='black')

#### Q2 DRAFT FIG. 1 - MIXED ####
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
  sim=ode(y=y0,times=tstep,func=mixedQ2,parms=p)
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
  sim=ode(y=y0,times=tstep,func=mixedQ2,parms=p)
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
title(main = "mixed")

#### Q2 DRAFT FIG. 3 - MIXED ####
#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. 

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
dfT=expand.grid(X=qEs,Y=sto)
dfT$A1=numeric(nrow(dfT))
dfT$A2=numeric(nrow(dfT))
dfT$J1=numeric(nrow(dfT))
dfT$J2=numeric(nrow(dfT))

for(i in 1:nrow(dfT)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(4, length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(dfT$X[i], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfT$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=mixedQ2,parms=p)
  dfT$A1[i]=sim[nrow(sim)-1,2]
  dfT$A2[i]=sim[nrow(sim)-1,3]
  dfT$J1[i]=sim[nrow(sim)-1,4]
  dfT$J2[i]=sim[nrow(sim)-1,5]
}

dfT$diff=dfT$A1-dfT$A2
dfT$qNorm=(dfT$X-min(dfT$X))/(max(dfT$X)-min(dfT$X))
dfT$sNorm=(dfT$Y-min(dfT$Y))/(max(dfT$Y)-min(dfT$Y))
#dfT$outcome=ifelse(dfT$A1-dfT$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfT)+theme_classic()
vz+geom_tile(aes(x=dfT$qNorm, y=dfT$sNorm, fill=dfT$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 \nDominance", breaks=c(min(dfT$diff), 1, max(dfT$diff)), labels=c("sp2", "even", "sp1"))+
  labs(x="Sp2 harvest Rate", y="Sp1 stocked fish")+
  ggtitle("Mixed")
# trying to make some heat maps to explore how stocking can maintain sp1 dominance
library(deSolve)
library(ggplot2)
library(ggpubr)
library(metR)
rm(list=ls())
source('q2Func.R')

tstep=1:300


####  BASIC MODEL BEHAVIOR ####

#approx funs for changing parms through time
#qE1Fun=approxfun(x=tstep,y=c(seq(0,by=0.1111,length.out = 113),rep(0,187)))
# qE2Fun=approxfun(x=tstep,y=seq(0,5,length.out = length(tstep)))
qE1Fun=approxfun(x=tstep,y=rep(1.8, length(tstep)))
qE2Fun=approxfun(x=tstep,y=rep(1.8, length(tstep)))

# h1Fun=approxfun(x=tstep,y=seq(0,10,length.out = length(tstep)))
# h2Fun=approxfun(x=tstep,y=seq(0,10,length.out = length(tstep)))
h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))

st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))

#single model run
p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.005,cJ1J2=0.003,v1=1,f1=2,
    s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.004,cJ2J1=0.003,v2=1,f2=2)
y0=c(10,10,0,0)
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
sim[max(tstep),]
plot(sim[,1],sim[,2], type = 'l', ylim=c(0,max(sim[,2:3], na.rm = T)))
lines(sim[,1],sim[,3],col ='green')

#### OUTCOMES W/HYSTERESIS ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,8,length.out = 15)
sto=seq(0,200, length.out = 15)
df=expand.grid(X=qEs,Y=sto)
df$A1=numeric(nrow(df))
df$A2=numeric(nrow(df))
df$J1=numeric(nrow(df))
df$J2=numeric(nrow(df))


for(i in 1:nrow(df)){
  tstep=1:300
  qE1Fun=approxfun(x=tstep,y=rep(df$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(df$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.005,cJ1J2=0.003,v1=1,f1=2,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.004,cJ2J1=0.003,v2=1,f2=2)
  y0=c(10,10,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df$A1[i]=sim[nrow(sim)-1,2]
  df$A2[i]=sim[nrow(sim)-1,3]
  df$J1[i]=sim[nrow(sim)-1,4]
  df$J2[i]=sim[nrow(sim)-1,5]
}

vz=ggplot(data = df)+theme_classic()
a=vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$A1))+scale_fill_gradient(low ="blue", high = "yellow", name="")+
  labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")
vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$A2))+scale_fill_gradient(low ="blue", high = "yellow", name="")+
  labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")
vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$Y/df$A1))+scale_fill_gradient(low ="blue", high = "yellow", name="")+
  labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")


##### OUTCOMES W/O HYSTERESIS ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,8,length.out = 15)
sto=seq(0,200, length.out = 15)
dfwo=expand.grid(X=qEs,Y=sto)
dfwo$A1=numeric(nrow(dfwo))
dfwo$A2=numeric(nrow(dfwo))
dfwo$J1=numeric(nrow(dfwo))
dfwo$J2=numeric(nrow(dfwo))

source('q2Func.R')

for(i in 1:nrow(dfwo)){
  tstep=1:300
  qE1Fun=approxfun(x=tstep,y=rep(dfwo$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfwo$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,cJ1A1=0.004,cJ1A2=0.001,cJ1J2=0.001,v1=1,f1=2,
      s2=0.1,cJ2A2=0.004,cJ2A1=0.001,cJ2J1=0.001,v2=1,f2=2)
  y0=c(10,10,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo$A1[i]=sim[nrow(sim)-1,2]
  dfwo$A2[i]=sim[nrow(sim)-1,3]
  dfwo$J1[i]=sim[nrow(sim)-1,4]
  dfwo$J2[i]=sim[nrow(sim)-1,5]
}

vzwo=ggplot(data = dfwo)+theme_classic()
b=vzwo+geom_tile(aes(x=dfwo$X, y=dfwo$Y, fill=dfwo$A1))+
  scale_fill_gradient(low ="blue", high = "yellow", name="")+
  labs(x="Harvest Rate", y="Stocked fish", title = "No Hysteresis Present")
vzwo+geom_tile(aes(x=dfwo$X, y=dfwo$Y, fill=dfwo$A2))+
  scale_fill_gradient(low ="blue", high = "yellow", name="")+
  labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")


#comparing with and without hysteresis plots together using cowplot package
dev.new()
ggarrange(a,b,labels = c("A","B"), nrow = 1,ncol = 2)


#### NON-EQUILIBRIUM TRY ####
#thinking of a new way to show the effects of hysteresis on management outcomes using plots somewhat like Biggs. Goal is to start system in an undesirable state, see what amount of stocking is necessary to push it over the edge to a desirable state. Then start system in a desirable state and see what amount of harvest can be tolerated while still staying in that state? Or what combo of regs and stocking?


#### WITH HYSTERESIS
tstep=1:300
h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
qE1Fun=approxfun(x=tstep,y=c(seq(0,by=.74,length.out = 100),rep(0,200)))
qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
p=c(s1=0.1,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
    s2=0.1,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1,f2=2)
y0=c(100,10,0,0)
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "Year",ylab = "pop. size")
lines(sim[,1],sim[,3],col='black',lwd=2)
lines(sim[,1],qE1Fun(sim[,1]),lty=2)
legend("topleft", legend = c("sp1","sp2"), col = c('grey','black'),lty = 1,lwd=2,bty = "n")

#### WITHOUT HYSTERESIS
tstep=1:300
h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
qE1Fun=approxfun(x=tstep,y=c(seq(0,by=.74,length.out = 100),rep(0,200)))
qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
p=c(s1=0.1,cJ1A1=0.004,cJ1A2=0.001,cJ1J2=0.001,v1=1,f1=2,
    s2=0.1,cJ2A2=0.004,cJ2A1=0.001,cJ2J1=0.001,v2=1,f2=2)
y0=c(100,10,0,0)
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "Year",ylab = "pop. size")
lines(sim[,1],sim[,3],col='black',lwd=2)
lines(sim[,1],qE1Fun(sim[,1]),lty=2)
legend("topleft", legend = c("sp1","sp2"), col = c('grey','black'),lty = 1,lwd=2,bty = "n")


#### figure 1 ####

# mangement intervention in time
tstep=1:300
h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
qE1Fun=approxfun(x=tstep,y=c(seq(0,by=.74,length.out = 94),rep(0,206)))
qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
p=c(s1=0.1,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
    s2=0.1,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1,f2=2)
y0=c(100,10,0,0)
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#intervention one year later has no effect
qE1Fun=approxfun(x=tstep,y=c(seq(0,by=.74,length.out = 95),rep(0,205)))
sim2=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#plotting
par(mfcol=c(2,1), mar=c(1,1,1,1), oma=c(4,4,1,1))
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim[,1],sim[,3],col='black',lwd=2)
lines(sim[,1],qE1Fun(sim[,1]),lty=2)
text(x=225, y=100, labels = "94 years of harvest \nbefore intervention", cex=0.75)
legend("topleft", legend = c("sp1","sp2", "harvest rate"), col = c('grey','black','black'),lty = c(1,1,2),lwd=2,bty = "n")
plot(sim2[,1],sim2[,2],type='l',ylim=c(0,max(sim2[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim2[,1],sim2[,3],col='black',lwd=2)
lines(sim2[,1],qE1Fun(sim2[,1]),lty=2)
text(x=225, y=100, labels = "95 years of harvest \nbefore intervention", cex=0.75)
mtext("Year", side = 1, outer = T, line = 2.5)
mtext("Population Size", side = 2, outer = T, line = 2.5)

#legend("topleft", legend = c("sp1","sp2"), col = c('grey','black'),lty = 1,lwd=2,bty = "n")

##### HEATMAPS #####

  ##### MAINTAIN W/ HYST. ####
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
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
  y0=c(1000,900,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df$A1[i]=sim[nrow(sim)-1,2]
  df$A2[i]=sim[nrow(sim)-1,3]
  df$J1[i]=sim[nrow(sim)-1,4]
  df$J2[i]=sim[nrow(sim)-1,5]
}

df$diff=df$A1-df$A2
df$outcome=ifelse(df$A1-df$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = df)+theme_classic()
a=vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="", breaks=c(min(df$diff), 1, max(df$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")+
  geom_contour(aes(x=df$X, y=df$Y, z=df$diff), breaks = c(minDiff))
  
#a=vz+geom_tile(aes(x=df$X, y=df$Y), fill=df$outcome)+
#  labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")
# vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$A2))+scale_fill_gradient(low ="blue", high = "yellow", name="")+
#   labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")
# vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$Y/df$A1))+scale_fill_gradient(low ="blue", high = "yellow", name="")+
#   labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")

#### FLIP W/ HYST. ####

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
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
  y0=c(900,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df2$A1[i]=sim[nrow(sim)-1,2]
  df2$A2[i]=sim[nrow(sim)-1,3]
  df2$J1[i]=sim[nrow(sim)-1,4]
  df2$J2[i]=sim[nrow(sim)-1,5]
}

df2$diff=df2$A1-df2$A2
df2$outcome=ifelse(df2$A1-df2$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = df2)+theme_classic()
b=vz+geom_tile(aes(x=df2$X, y=df2$Y, fill=df2$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="")+
  labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")+
  geom_contour(aes(x=df2$X, y=df2$Y, z=df2$diff), breaks = c(minDiff))
#b=vz+geom_tile(aes(x=df2$X, y=df2$Y), fill=df2$outcome)+
#  labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")

##### MAINTAIN W/O HYST. ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
dfwo=expand.grid(X=qEs,Y=sto)
dfwo$A1=numeric(nrow(dfwo))
dfwo$A2=numeric(nrow(dfwo))
dfwo$J1=numeric(nrow(dfwo))
dfwo$J2=numeric(nrow(dfwo))

source('q2Func.R')

for(i in 1:nrow(dfwo)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfwo$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfwo$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.9,cJ1A1=0.004,cJ1A2=0.001,cJ1J2=0.001,v1=1,f1=6,
      s2=0.1,m2=0.9,cJ2A2=0.004,cJ2A1=0.001,cJ2J1=0.001,v2=1,f2=6)
  y0=c(1000,900,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo$A1[i]=sim[nrow(sim)-1,2]
  dfwo$A2[i]=sim[nrow(sim)-1,3]
  dfwo$J1[i]=sim[nrow(sim)-1,4]
  dfwo$J2[i]=sim[nrow(sim)-1,5]
}

dfwo$diff=dfwo$A1-dfwo$A2
dfwo$outcome=ifelse(dfwo$A1-dfwo$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfwo)+theme_classic()
c=vz+geom_tile(aes(x=dfwo$X, y=dfwo$Y, fill=dfwo$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="")+
  labs(x="Harvest Rate", y="Stocked fish", title = "No Hysteresis Present")+
  geom_contour(aes(x=dfwo$X, y=dfwo$Y, z=dfwo$diff), breaks = c(minDiff))
# c=vzwo+geom_tile(aes(x=dfwo$X, y=dfwo$Y), fill=dfwo$outcome)+
#   scale_fill_gradient(low ="blue", high = "yellow", name="")+
#   labs(x="Harvest Rate", y="Stocked fish", title = "No Hysteresis Present")
# vzwo+geom_tile(aes(x=dfwo$X, y=dfwo$Y, fill=dfwo$A2))+
#   scale_fill_gradient(low ="blue", high = "yellow", name="")+
#   labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")

#### FLIP W/O HYST. ####

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
dfwo2=expand.grid(X=qEs,Y=sto)
dfwo2$A1=numeric(nrow(dfwo2))
dfwo2$A2=numeric(nrow(dfwo2))
dfwo2$J1=numeric(nrow(dfwo2))
dfwo2$J2=numeric(nrow(dfwo2))

source('q2Func.R')

for(i in 1:nrow(dfwo2)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfwo2$X[i], length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfwo2$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.9,cJ1A1=0.004,cJ1A2=0.001,cJ1J2=0.001,v1=1,f1=6,
      s2=0.1,m2=0.9,cJ2A2=0.004,cJ2A1=0.001,cJ2J1=0.001,v2=1,f2=6)
  y0=c(900,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo2$A1[i]=sim[nrow(sim)-1,2]
  dfwo2$A2[i]=sim[nrow(sim)-1,3]
  dfwo2$J1[i]=sim[nrow(sim)-1,4]
  dfwo2$J2[i]=sim[nrow(sim)-1,5]
}

dfwo2$diff=dfwo2$A1-dfwo2$A2
dfwo2$outcome=ifelse(dfwo2$A1-dfwo2$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfwo2)+theme_classic()
d=vz+geom_tile(aes(x=dfwo2$X, y=dfwo2$Y, fill=dfwo2$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="")+
  labs(x="Harvest Rate", y="Stocked fish", title = " No Hysteresis Present")+
  geom_contour(aes(x=dfwo2$X, y=dfwo2$Y, z=dfwo2$diff), breaks = c(minDiff))
# d=vzwo+geom_tile(aes(x=dfwo2$X, y=dfwo2$Y), fill=dfwo2$outcome)+
#   scale_fill_gradient(low ="blue", high = "yellow", name="")+
#   labs(x="Harvest Rate", y="Stocked fish", title = "No Hysteresis Present")
# vzwo+geom_tile(aes(x=dfwo2$X, y=dfwo2$Y, fill=dfwo2$A2))+
#   scale_fill_gradient(low ="blue", high = "yellow", name="")+
#   labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")

#comparing with and without hysteresis plots together
dev.new()
ggarrange(a,b,c,d,labels = c("A","B","C","D"), nrow = 2,ncol = 2, common.legend = T, legend = 'bottom')

#### DELAY A FLIP ####
#figure on delaying a transition

tstep=1:100
h1Fun=approxfun(x=tstep,y=c(seq(3,0, length.out = 50), rep(0,50)))
h2Fun=approxfun(x=tstep,y=c(seq(10,5, length.out = 50),rep(5,50)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,100)))
st1Fun=approxfun(x=tstep,y=c(rep(0,100)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
p=c(s1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
    s2=0.5,cJ2A2=0.001,cJ2A1=0.05,cJ2J1=0.003,v2=1,f2=2)
y0=c(200,100,0,0)
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#stocking delays the transition to sp2
tstep=1:100
h1Fun=approxfun(x=tstep,y=c(seq(3,0, length.out = 50), rep(0,50)))
h2Fun=approxfun(x=tstep,y=c(seq(10,5, length.out = 50),rep(5,50)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,100)))
st1Fun=approxfun(x=tstep,y=c(rep(20,100)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))

sim2=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#plotting
par(mfcol=c(2,1), mar=c(1,1,1,1), oma=c(4,4,1,1))
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim[,1],sim[,3],col='black',lwd=2)
text(x=60,y=1000,labels = "No stocking")
plot(sim2[,1],sim2[,2],type='l',ylim=c(0,max(sim2[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim2[,1],sim2[,3],col='black',lwd=2)
legend("topleft", legend = c("sp1","sp2"), col = c('grey','black'),lty = c(1,1),lwd=2,bty = "n")
text(x=60, y=1000, labels = "Stocking ps1 \n20 indv/yr")
mtext("Year", side = 1, outer = T, line = 2.5)
mtext("Population Size", side = 2, outer = T, line = 2.5)

# #### OPTIM FUNC ####
# source('q2Func.R')
# 
# optA1=function(parms){
#   tstep=1:300
#   qE1Fun=approxfun(x=tstep,y=rep(parms[1], length(tstep)))
#   qE2Fun=approxfun(x=tstep,y=rep(parms[2], length(tstep)))
#   st1Fun=approxfun(x=tstep,y=rep(parms[3],length(tstep)))
#   st2Fun=approxfun(x=tstep,y=rep(parms[4],length(tstep)))
#   p=c(s1=0.1,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
#       s2=0.1,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1,f2=2)
#   y0=c(10,10,0,0)
#   sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
#   return(sim[nrow(sim)-1,2]) #A1 abundance, what I want to maximize
# }
# 
# parms=c(1.8,0,200,0)
# out=optim(par=parms, optA1,control = list(fnscale=-1))
# out


#### Diagnostic plot ####

# running to equilibrium over a range of harvest values

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
  
  p=c(s1=0.1,m1=0.9,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.002,v1=1,f1=2,
      s2=0.1,m2=0.9,cJ2A2=0.001,cJ2A1=0.5,cJ2J1=0.002,v2=1,f2=2)
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
  
  p=c(s1=0.1,m1=0.9,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.002,v1=1,f1=2,
      s2=0.1,m2=0.9,cJ2A2=0.001,cJ2A1=0.5,cJ2J1=0.002,v2=1,f2=2)
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

##### HEATMAPS_NEW 9.28.2020 #####

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
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
  y0=c(1000,100,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df$A1[i]=sim[nrow(sim)-1,2]
  df$A2[i]=sim[nrow(sim)-1,3]
  df$J1[i]=sim[nrow(sim)-1,4]
  df$J2[i]=sim[nrow(sim)-1,5]
}

df$diff=df$A1-df$A2
df$outcome=ifelse(df$A1-df$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = df)+theme_classic()
a=vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="", breaks=c(min(df$diff), 1, max(df$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "No reduction in sp2 mortality, sp1>sp2")+
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
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df2$A1[i]=sim[nrow(sim)-1,2]
  df2$A2[i]=sim[nrow(sim)-1,3]
  df2$J1[i]=sim[nrow(sim)-1,4]
  df2$J2[i]=sim[nrow(sim)-1,5]
}

df2$diff=df2$A1-df2$A2
df2$outcome=ifelse(df2$A1-df2$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = df2)+theme_classic()
b=vz+geom_tile(aes(x=df2$X, y=df2$Y, fill=df2$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="", breaks=c(min(df2$diff), 1, max(df2$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "No reduction in sp2 mortality sp2>sp1")+
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
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
  y0=c(1000,100,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo$A1[i]=sim[nrow(sim)-1,2]
  dfwo$A2[i]=sim[nrow(sim)-1,3]
  dfwo$J1[i]=sim[nrow(sim)-1,4]
  dfwo$J2[i]=sim[nrow(sim)-1,5]
}

dfwo$diff=dfwo$A1-dfwo$A2
dfwo$outcome=ifelse(dfwo$A1-dfwo$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfwo)+theme_classic()
c=vz+geom_tile(aes(x=dfwo$X, y=dfwo$Y, fill=dfwo$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="", breaks=c(min(dfwo$diff), 1, max(dfwo$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "Sp2 CPUE=2, sp1>sp2")+
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
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo2$A1[i]=sim[nrow(sim)-1,2]
  dfwo2$A2[i]=sim[nrow(sim)-1,3]
  dfwo2$J1[i]=sim[nrow(sim)-1,4]
  dfwo2$J2[i]=sim[nrow(sim)-1,5]
}

dfwo2$diff=dfwo2$A1-dfwo2$A2
dfwo2$outcome=ifelse(dfwo2$A1-dfwo2$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfwo2)+theme_classic()
d=vz+geom_tile(aes(x=dfwo2$X, y=dfwo2$Y, fill=dfwo2$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="", breaks=c(min(dfwo2$diff), 1, max(dfwo2$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = " Sp2 CPUE=2 sp2>sp1")+
  geom_contour(aes(x=dfwo2$X, y=dfwo2$Y, z=dfwo2$diff), breaks = c(minDiff))

#comparing with and without hysteresis plots together using cowplot package
#dev.new()
ggarrange(a,b,c,d,labels = c("A","B","C","D"), nrow = 2,ncol = 2, common.legend = T, legend = 'bottom')


#### TRADEOFFS ####

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
  qE1Fun=approxfun(x=tstep,y=rep(2, length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(dfT$X[i], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfT$Y[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
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
dfT$outcome=ifelse(dfT$A1-dfT$A2 < minDiff,"darkgreen", "darkred")
vz=ggplot(data = dfT)+theme_classic()
d=vz+geom_tile(aes(x=dfT$qNorm, y=dfT$sNorm, fill=dfT$A1))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 \nDominance", breaks=c(min(dfT$diff), 1, max(dfT$diff)), labels=c("sp2", "even", "sp1"))+
  labs(x="Sp2 harvest Rate", y="Sp1 stocked fish")#+
  #geom_contour(aes(x=dfT$X, y=dfT$Y, z=dfT$diff))+
  #geom_label_contour(aes(x=dfT$X, y=dfT$Y, z=dfT$diff))
d

plot(dfT$diff,dfT$qNorm/dfT$sNorm, pch=16)
plot(dfT$qNorm/dfT$sNorm, pch=16)

# different version of the above figure
# I run model to equilibrium first over a  range of sp1 stocking and then separately over a range of sp2 harvest and then plot those against each other

qEs=seq(0,8,length.out = 30)
dfHarv2=data.frame(qEs=qEs,A1=numeric(length(qEs)),A2=numeric(length(qEs)),J1=numeric(length(qEs)),J2=numeric(length(qEs)))
sto=seq(0,2000, length.out = 30)
dfSto1=data.frame(sto=sto,A1=numeric(length(sto)),A2=numeric(length(sto)),J1=numeric(length(sto)),J2=numeric(length(sto)))

for(i in 1:nrow(dfHarv2)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(2, length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(dfHarv2$qEs[i], length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfHarv2$A1[i]=sim[nrow(sim)-1,2]
  dfHarv2$A2[i]=sim[nrow(sim)-1,3]
  dfHarv2$J1[i]=sim[nrow(sim)-1,4]
  dfHarv2$J2[i]=sim[nrow(sim)-1,5]
}

for(i in 1:nrow(dfSto1)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(2, length(tstep)))
  qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(dfSto1$sto[i],length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
  y0=c(100,1000,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfSto1$A1[i]=sim[nrow(sim)-1,2]
  dfSto1$A2[i]=sim[nrow(sim)-1,3]
  dfSto1$J1[i]=sim[nrow(sim)-1,4]
  dfSto1$J2[i]=sim[nrow(sim)-1,5]
}


dfHarv2$diff=dfHarv2$A1-dfHarv2$A2
dfHarv2$qNorm=(dfHarv2$qEs-min(dfHarv2$qEs))/(max(dfHarv2$qEs)-min(dfHarv2$qEs))
dfHarv2$outcome=ifelse(dfHarv2$A1-dfHarv2$A2 < minDiff,"darkgreen", "darkred")
dfSto1$diff=dfSto1$A1-dfSto1$A2
dfSto1$sNorm=(dfSto1$sto-min(dfSto1$sto))/(max(dfSto1$sto)-min(dfSto1$sto))
dfSto1$outcome=ifelse(dfSto1$A1-dfSto1$A2 < minDiff,"darkgreen", "darkred")
 #got this far then wasn't sure what to do, didn't work out the way I thought
vz=ggplot(data = dfHarv2)+theme_classic()
d=vz+geom_tile(aes(x=dfHarv2$X, y=dfHarv2$Y, fill=dfHarv2$diff))+
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 \nDominance", breaks=c(min(dfHarv2$diff), 1, max(dfHarv2$diff)), labels=c("sp2", "even", "sp1"))+
  labs(x="Sp2 harvest Rate", y="Sp1 stocked fish")#+
#geom_contour(aes(x=dfT$X, y=dfT$Y, z=dfT$diff))+
#geom_label_contour(aes(x=dfT$X, y=dfT$Y, z=dfT$diff))
d



#### DELAY A FLIP_new 9.28.2020 ####

#figure on delaying a transition using stocking

#running to sp1 dominating before habitat decline
tstep=1:200
h1Fun=approxfun(x=tstep,y=c(rep(10,200)))
h2Fun=approxfun(x=tstep,y=c(rep(10,200)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
qE2Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
    s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.05,cJ2J1=0.003,v2=1,f2=2)
y0=c(500,100,0,0)
simPre=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#check to make sure it's at equilibrium
# plot(simPre[,1], simPre[,2],type='l',ylim=c(0,max(simPre[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
# lines(simPre[,1],simPre[,3],col='black',lwd=2)

#habitat decline and stocking
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0, length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,20, length.out = 100),rep(20,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)))
st1Fun=approxfun(x=tstep,y=c(rep(0,500)))
st2Fun=approxfun(x=tstep,y=c(rep(0,200),rep(10,50),rep(0,250)))#this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amoutn of harvest that declines to 0
p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
    s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1,f2=2)
y0=simPre[199,2:5]
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0, length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,20, length.out = 100),rep(20,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)))
st1Fun=approxfun(x=tstep,y=c(rep(0,100),rep(80,400)))
st2Fun=approxfun(x=tstep,y=c(rep(0,200),rep(10,50),rep(0,250))) #this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amoutn of harvest that declines to 0

sim2=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#plotting
par(mfcol=c(2,1), mar=c(1,1,1,1), oma=c(4,4,1,1))
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim[,1],sim[,3],col='black',lwd=2)
text(x=375,y=1500,labels = "No stocking")
legend("topleft", legend = c("sp1","sp2"), col = c('grey','black'),lty = c(1,1),lwd=2,bty = "n")

plot(sim2[,1],sim2[,2],type='l',ylim=c(0,max(sim2[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim2[,1],sim2[,3],col='black',lwd=2)
text(x=450, y=1500, labels = "Stocking sp1 \n80 indv/yr")
mtext("Year", side = 1, outer = T, line = 2.5)
mtext("Population Size", side = 2, outer = T, line = 2.5)


#figure on delaying a transition using harvest sp2

#running to sp1 dominating before habitat decline
tstep=1:200
h1Fun=approxfun(x=tstep,y=c(rep(10,200)))
h2Fun=approxfun(x=tstep,y=c(rep(10,200)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
qE2Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
    s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.05,cJ2J1=0.003,v2=1,f2=2)
y0=c(500,100,0,0)
simPre=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#check to make sure it's at equilibrium
# plot(simPre[,1], simPre[,2],type='l',ylim=c(0,max(simPre[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
# lines(simPre[,1],simPre[,3],col='black',lwd=2)

#habitat decline and stocking
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0, length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,20, length.out = 100),rep(20,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)))
st1Fun=approxfun(x=tstep,y=c(rep(0,500)))
st2Fun=approxfun(x=tstep,y=c(rep(0,200),rep(10,50),rep(0,250)))#this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amount of harvest that declines to 0
p=c(s1=0.1,m1=0.9,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
    s2=0.1,m2=0.9,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1,f2=2)
y0=simPre[199,2:5]
sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,0, length.out = 100), rep(0,300)))
h2Fun=approxfun(x=tstep,y=c(rep(10,100),seq(10,20, length.out = 100),rep(20,300)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(4,200)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,4,length.out = 200), rep(2,200)))
st1Fun=approxfun(x=tstep,y=c(rep(0,500)))
st2Fun=approxfun(x=tstep,y=c(rep(0,200),rep(10,50),rep(0,250))) #this is to simulate increase in fecundity, could change the model to make f time dependent, or have some amoutn of harvest that declines to 0

sim2=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)

#plotting
par(mfcol=c(2,1), mar=c(1,1,1,1), oma=c(4,4,1,1))
plot(sim[,1],sim[,2],type='l',ylim=c(0,max(sim[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim[,1],sim[,3],col='black',lwd=2)
text(x=375,y=1500,labels = "No Harvest Species 2")
legend("topleft", legend = c("sp1","sp2"), col = c('grey','black'),lty = c(1,1),lwd=2,bty = "n")

plot(sim2[,1],sim2[,2],type='l',ylim=c(0,max(sim2[,2:3],na.rm = T)),col='grey',lwd=2,xlab = "",ylab = "")
lines(sim2[,1],sim2[,3],col='black',lwd=2)
text(x=300, y=1250, labels = "Harvest Species 2")
mtext("Year", side = 1, outer = T, line = 2.5)
mtext("Population Size", side = 2, outer = T, line = 2.5)


#### 3D PLOT ####

qEs=seq(0,8,length.out = 30)
sto=seq(0,2000, length.out = 30)
df=expand.grid(X=qEs,Y=sto)
df$A1=numeric(nrow(df))
df$A2=numeric(nrow(df))
df$J1=numeric(nrow(df))
df$J2=numeric(nrow(df))
minDiff=100


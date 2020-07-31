# trying to make some heat maps to explore how stocking can maintain sp1 dominance
library(deSolve)
library(ggplot2)
library(ggpubr)
rm(list=ls())
source('q2Func.R')

tstep=1:300


####  BASIC MODEL BEHAVIOR ####



#approx funs for changing parms through time
qE1Fun=approxfun(x=tstep,y=c(seq(0,by=0.1111,length.out = 113),rep(0,187)))
# qE2Fun=approxfun(x=tstep,y=seq(0,5,length.out = length(tstep)))
# qE1Fun=approxfun(x=tstep,y=rep(1.8, length(tstep)))
qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)))

# h1Fun=approxfun(x=tstep,y=seq(0,10,length.out = length(tstep)))
# h2Fun=approxfun(x=tstep,y=seq(0,10,length.out = length(tstep)))
h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))

st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))

##### OUTCOMES W/HYSTERESIS ####
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
  p=c(s1=0.1,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=2,
      s2=0.1,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=3)
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
  p=c(s1=0.1,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
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
  p=c(s1=0.1,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,f1=4,
      s2=0.1,cJ2A2=0.002,cJ2A1=0.003,cJ2J1=0.003,v2=1,f2=6)
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
  p=c(s1=0.1,cJ1A1=0.004,cJ1A2=0.001,cJ1J2=0.001,v1=1,f1=6,
      s2=0.1,cJ2A2=0.004,cJ2A1=0.001,cJ2J1=0.001,v2=1,f2=6)
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
  p=c(s1=0.1,cJ1A1=0.004,cJ1A2=0.001,cJ1J2=0.001,v1=1,f1=6,
      s2=0.1,cJ2A2=0.004,cJ2A1=0.001,cJ2J1=0.001,v2=1,f2=6)
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
  labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")+
  geom_contour(aes(x=dfwo2$X, y=dfwo2$Y, z=dfwo2$diff), breaks = c(minDiff))
# d=vzwo+geom_tile(aes(x=dfwo2$X, y=dfwo2$Y), fill=dfwo2$outcome)+
#   scale_fill_gradient(low ="blue", high = "yellow", name="")+
#   labs(x="Harvest Rate", y="Stocked fish", title = "No Hysteresis Present")
# vzwo+geom_tile(aes(x=dfwo2$X, y=dfwo2$Y, fill=dfwo2$A2))+
#   scale_fill_gradient(low ="blue", high = "yellow", name="")+
#   labs(x="Harvest Rate", y="Stocked fish", title = "Hysteresis Present")

#comparing with and without hysteresis plots together using cowplot package
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

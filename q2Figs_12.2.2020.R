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
vz=ggplot(data = df)+theme_classic()
# a=vz+geom_tile(aes(x=df$X, y=df$Y, fill=df$diff))+
#   scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(df$diff), 1, max(df$dif)), labels=c("sp2", "even", "sp1"))+
#   labs(x="Harvest Rate", y="Stocked fish", title = "")+
#   geom_contour(aes(x=df$X, y=df$Y, z=df$diff), breaks = c(minDiff))
a=vz+geom_point(aes(x=df$X,df$A1,color=df$Y)) #stocking boosts curve up across the board

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
  scale_fill_gradient2(low="darkred", high="darkgreen", mid="white", midpoint=1, name="Sp1 Dominance", breaks=c(min(dfwo2$diff), 1, max(dfwo2$dif)), labels=c("sp2", "even", "sp1"))+
  labs(x="Harvest Rate", y="Stocked fish", title = "")+
  geom_contour(aes(x=dfwo2$X, y=dfwo2$Y, z=dfwo2$diff), breaks = c(minDiff))

#trying some new plots

#plot isoclines for all scenarios on one plot with x=harvest, y= stocking, line where 1>2
df$mod=rep("Maintain 1, ignore 2",nrow(df))
df2$mod=rep("Flip to 1, ignore 2",nrow(df2))
dfwo$mod=rep("Maintain 1, harv 2",nrow(dfwo))
dfwo2$mod=rep("Flip to 1, harv 2",nrow(dfwo2))
allSen=rbind(df,df2,dfwo,dfwo2)
vzI=ggplot(data=allSen, aes(x=allSen$X,y=allSen$Y,color=allSen$mod))+theme_classic()+
  geom_contour(aes(z=allSen$diff),breaks = c(minDiff))+
  labs(x="Species 1 Harvest Rate", y="Species 1 Stocking",color="Scenario")+
  theme(legend.position = 'bottom')+guides(color=guide_legend(nrow=2,byrow = T))
vzI
#comparing with and without hysteresis plots together 
#dev.new()
ggarrange(a,b,c,d,labels = c("A","B","C","D"), nrow = 2,ncol = 2, common.legend = T, legend = 'bottom')

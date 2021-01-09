## code to think about new ways to visualize some of our results

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



#demonstrate alternative stable states
store=data.frame(qEs=seq(0,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(5000,500,0,0)
tstep=1:300
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
store2=data.frame(qEs=seq(0,8,length.out=30),A1=0,A2=0,J1=0,J2=0)
y0=c(500,5000,0,0)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(1.8,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}

## ggplot version of the normal 'demonstrate alt.stable states' plots
panA=as.data.frame(rbind(cbind(store$qEs,store$A1,rep("A1",nrow(store))),cbind(store$qEs,store$A2,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="",labels=c("Species 1, initially dominant", "Species 2"))+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = c(.5,.75))
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="",labels=c("Species 1", "Species 2, initially dominant"))+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10),
        legend.position = c(.5,.75))
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90, size=14),
                bottom = text_grob("Harvest Rate (qE)", size=14))

##### HEATMAPS #####

tstep=1:300

h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)
h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)),rule = 2)

##### MAINTAIN W/O sp2 harv ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,8,length.out = 30)
sto=seq(0,20000, length.out = 30)
df=expand.grid(X=qEs,Y=sto)
df$A1=numeric(nrow(df))
df$A2=numeric(nrow(df))
df$J1=numeric(nrow(df))
df$J2=numeric(nrow(df))
minDiff=100

for(i in 1:nrow(df)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(df$X[i], length(tstep)),rule = 2)
  qE2Fun=approxfun(x=tstep,y=rep(0, length(tstep)),rule = 2)
  st1Fun=approxfun(x=tstep,y=rep(df$Y[i],length(tstep)),rule = 2)
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  y0=c(5000,500,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  df$A1[i]=sim[nrow(sim)-1,2]
  df$A2[i]=sim[nrow(sim)-1,3]
  df$J1[i]=sim[nrow(sim)-1,4]
  df$J2[i]=sim[nrow(sim)-1,5]
}

df$diff=df$A1-df$A2
df$outcome=ifelse(df$A1-df$A2 < minDiff,"darkgreen", "darkred")



##### MAINTAIN W/ sp2 harv ####
#matrix to hold output, starting with different harvest levels on each species

qEs=seq(0,8,length.out = 30)
sto=seq(0,20000, length.out = 30)
dfwo=expand.grid(X=qEs,Y=sto)
dfwo$A1=numeric(nrow(dfwo))
dfwo$A2=numeric(nrow(dfwo))
dfwo$J1=numeric(nrow(dfwo))
dfwo$J2=numeric(nrow(dfwo))

for(i in 1:nrow(dfwo)){
  tstep=1:100
  qE1Fun=approxfun(x=tstep,y=rep(dfwo$X[i], length(tstep)),rule = 2)
  qE2Fun=approxfun(x=tstep,y=rep(2, length(tstep)),rule = 2)
  st1Fun=approxfun(x=tstep,y=rep(dfwo$Y[i],length(tstep)),rule = 2)
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  y0=c(5000,500,0,0)
  sim=ode(y=y0,times=tstep,func=simBiggsQ2,parms=p)
  dfwo$A1[i]=sim[nrow(sim)-1,2]
  dfwo$A2[i]=sim[nrow(sim)-1,3]
  dfwo$J1[i]=sim[nrow(sim)-1,4]
  dfwo$J2[i]=sim[nrow(sim)-1,5]
}

dfwo$diff=dfwo$A1-dfwo$A2
dfwo$outcome=ifelse(dfwo$A1-dfwo$A2 < minDiff,"darkgreen", "darkred")
df$mod=rep("Maintain sp 1, ignore sp 2",nrow(df))
dfwo$mod=rep("Maintain sp 1, harv sp 2",nrow(dfwo))
allSen=rbind(df,dfwo)

vzI=ggplot(data=allSen, aes(x=allSen$X,y=allSen$Y,linetype=allSen$mod))+theme_classic()+
  geom_contour(aes(z=allSen$diff),breaks = c(minDiff), color='black', size=1.7)+
  labs(x="Species 1 Harvest Rate", y="Species 1 Stocking",linetype="Scenario")+
  theme(legend.position = 'right') + theme(legend.position = 'none')
vzI


### TRADEOFFS FIGURE ####
#new figure 3, isoclines at different harvests

#plot to look at the cost/benefits of stocking or predator reduction for managing a focal species. 

minDiff=100
qEs=rep(seq(0,8,length.out = 30),3)
sto=rep(seq(0,20000, length.out = 30),3)
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
  p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
  y0=c(5000,500,0,0)
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
  theme(legend.position = c(.75,.75))+
  scale_linetype_manual(values = c("solid","dashed","twodash"))
  #xlim(0,8)+
  #ylim(0,2000)
vzT


#### DELAY A TRANSITION ####
#same function as used in other figs but with option to change fecundity now
#comes from 'fecundityDecline.R' script
simBiggsR.a<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+((5000*dFun(t))*A1/(500+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

#running to sp1 dominating before habitat decline
tstep=1:100
h1Fun=approxfun(x=tstep,y=c(rep(8,length(tstep))),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,length(tstep))),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))),rule = 2)
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)),rule = 2)
dFun=approxfun(x=tstep,y=rep(1,length(tstep)),rule = 2)

p=c(s1=0.1,m1=0.1,cJ1A1=0.001,cJ1A2=0.05,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.1,cJ2A2=0.001,cJ2A1=0.03,cJ2J1=0.003,v2=1)
y0=c(5000,500,0,0)
simPre=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#Fecundity decline -> what will happen if no action taken
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
dFun=approxfun(x=tstep,y=c(seq(1,.01,length.out = 100),rep(.01,400)),rule = 2)

p=c(s1=0.1,m1=0.1,cJ1A1=0.002,cJ1A2=0.05,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.1,cJ2A2=0.002,cJ2A1=0.03,cJ2J1=0.003,v2=1)
#y0=simPre[nrow(simPre)-1,2:5]
sim=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(500,500)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 100),rep(.01,400)),rule = 2)

simS=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#harvesting delays the transition to sp2 - no stocking
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0.5,500)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 100),rep(.01,400)),rule = 2)

simH=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#harvesting and stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0.5,500)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(500,500)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 100),rep(.01,400)),rule = 2)

simB=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#plotting

# ggplot version

#do nothing
panA=as.data.frame(rbind(cbind(sim[,c(1,2)],rep("A1",nrow(sim))),cbind(sim[,c(1,3)],rep("A2",nrow(sim)))))
colnames(panA)=c("Time","Abund","sp")
panA$Abund=as.numeric(panA$Abund);panA$Time=as.numeric(panA$Time)

#harvest sp2
panB=as.data.frame(rbind(cbind(simH[,c(1,2)],rep("A1",nrow(simH))),cbind(simH[,c(1,3)],rep("A2",nrow(simH)))))
colnames(panB)=c("Time","Abund","sp")
panB$Time=as.numeric(panB$Time);panB$Abund=as.numeric(panB$Abund)

#stock sp1
panC=as.data.frame(rbind(cbind(simS[,c(1,2)],rep("A1",nrow(simS))),cbind(simS[,c(1,3)],rep("A2",nrow(simS)))))
colnames(panC)=c("Time","Abund","sp")
panC$Abund=as.numeric(panC$Abund);panC$Time=as.numeric(panC$Time)

#stock and harvest
panD=as.data.frame(rbind(cbind(simB[,c(1,2)],rep("A1",nrow(simB))),cbind(simB[,c(1,3)],rep("A2",nrow(simB)))))
colnames(panD)=c("Time","Abund","sp")
panD$Abund=as.numeric(panD$Abund);panD$Time=as.numeric(panD$Time)

a4=ggplot(data = panA,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())+
  ylim(-1,7000)
b4=ggplot(data = panB,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  ylim(-1,7000)
c4=ggplot(data = panC,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())+
  ylim(-1,7000)
d4=ggplot(data = panD,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())+
  ylim(-1,7000)

fig4=ggarrange(a4,b4,c4,d4,labels = c("A","B","C","D"),common.legend = T,legend = "top",label.x = .8,label.y = .9)
annotate_figure(fig4,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Time"))

## calculating the years of delay in each scenario
tp=function(x=matin){
  dat=as.data.frame(x[,1:3]);colnames(dat)=c("time","A1","A2")
  established = dat$A1[1] > dat$A2[1]
  comp = dat$A1>dat$A2
  if(is.finite(min(which(comp != established)))==T){
    return(min(which(comp != established)))
  }else{return(NA)}
}

simTP=tp(sim)
simHTP=tp(simH)
simSTP=tp(simS)
simbTP=tp(simB)

## what would doubling of the stocking or harvest do for delaying the flip?
#stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(1000,500)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 100),rep(.01,400)),rule = 2)

simS=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#harvesting delays the transition to sp2 - no stocking
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(1,500)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 100),rep(.01,400)),rule = 2)

simH=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#harvesting and stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
h2Fun=approxfun(x=tstep,y=c(rep(8,500)),rule = 2)
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2)
qE2Fun=approxfun(x=tstep,y=c(rep(1,500)),rule = 2)
st1Fun=approxfun(x=tstep,y=c(rep(1000,500)),rule = 2)
st2Fun=approxfun(x=tstep,y=c(rep(0,500)),rule = 2) 
dFun=approxfun(x=tstep,y=c(seq(1,0.01,length.out = 100),rep(.01,400)),rule = 2)

simB=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

DsimTP=tp(sim)
DsimHTP=tp(simH)
DsimSTP=tp(simS)
DsimbTP=tp(simB)

#logged version of the fig
# a4=ggplot(data = panA,aes(x=Time,y=log(Abund),color=sp))+theme_classic()+
#   geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.y = element_blank())
# b4=ggplot(data = panB,aes(x=Time,y=log(Abund),color=sp))+theme_classic()+
#   geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
#   theme(axis.title.x = element_blank(),
#         axis.text.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank())
# c4=ggplot(data = panC,aes(x=Time,y=log(Abund),color=sp))+theme_classic()+
#   geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank())
# d4=ggplot(data = panD,aes(x=Time,y=log(Abund),color=sp))+theme_classic()+
#   geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
#   theme(axis.title.x = element_blank(),
#         axis.title.y = element_blank(),
#         axis.text.y = element_blank())

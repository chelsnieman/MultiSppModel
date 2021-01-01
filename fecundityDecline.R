## Looking to try a scenarios figure with recruitment declines instead of refuge declines
## CJD 1.1.2021

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
  dFun=approxfun(x=tstep,y=rep(1,length(tstep)))
  
  p=c(s1=0.1,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  sim=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)
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
  dFun=approxfun(x=tstep,y=rep(1,length(tstep)))
  
  p=c(s1=0.1,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,
      s2=0.1,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1)
  sim=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}

## ggplot version of the normal 'demonstrate alt.stable states plots
panA=as.data.frame(rbind(cbind(store$qEs,store$A1,rep("A1",nrow(store))),cbind(store$qEs,store$A2,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10))
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 10))
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), common.legend = T,legend="top",label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90, size=14),
                bottom = text_grob("Harvest Rate (qE)", size=14))


#running to sp1 dominating before habitat decline
tstep=1:200
h1Fun=approxfun(x=tstep,y=c(rep(8,200)))
h2Fun=approxfun(x=tstep,y=c(rep(8,200)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
qE2Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st1Fun=approxfun(x=tstep,y=c(rep(0,length(tstep))))
st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
dFun=approxfun(x=tstep,y=rep(1,length(tstep)))

p=c(s1=0.1,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1)
y0=c(500,100,0,0)
simPre=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#Fecundity decline -> what will happen if no action taken
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)))
h2Fun=approxfun(x=tstep,y=c(rep(8,500)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)))
st1Fun=approxfun(x=tstep,y=c(rep(0,500)))
st2Fun=approxfun(x=tstep,y=c(rep(0,500)))
dFun=approxfun(x=tstep,y=c(seq(1,.01,length.out = 300),rep(.01,200)))

p=c(s1=0.1,m1=0.5,cJ1A1=0.002,cJ1A2=0.5,cJ1J2=0.003,v1=1,
    s2=0.1,m2=0.5,cJ2A2=0.002,cJ2A1=0.3,cJ2J1=0.003,v2=1)
y0=simPre[199,2:5]
sim=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)))
h2Fun=approxfun(x=tstep,y=c(rep(8,500)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,500)))
st1Fun=approxfun(x=tstep,y=c(rep(0,100),rep(15,400)))
st2Fun=approxfun(x=tstep,y=c(rep(0,500))) 
dFun=approxfun(x=tstep,y=c(seq(1,.01,length.out = 300),rep(.01,200)))

simS=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#harvesting delays the transition to sp2 - no stocking
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)))
h2Fun=approxfun(x=tstep,y=c(rep(8,500)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,1,length.out = 200), rep(1,200)))
st1Fun=approxfun(x=tstep,y=c(rep(0,500)))
st2Fun=approxfun(x=tstep,y=c(rep(0,500))) 
dFun=approxfun(x=tstep,y=c(seq(1,.01,length.out = 300),rep(.01,200)))

simH=ode(y=y0,times=tstep,func=simBiggsR.a,parms=p)

#harvesting and stocking delays the transition to sp2
tstep=1:500
h1Fun=approxfun(x=tstep,y=c(rep(8,500)))
h2Fun=approxfun(x=tstep,y=c(rep(8,500)))
qE1Fun=approxfun(x=tstep,y=c(rep(0,500)))
qE2Fun=approxfun(x=tstep,y=c(rep(0,100), seq(0,1,length.out = 200), rep(1,200)))
st1Fun=approxfun(x=tstep,y=c(rep(0,100),rep(15,400)))
st2Fun=approxfun(x=tstep,y=c(rep(0,500))) 
dFun=approxfun(x=tstep,y=c(seq(1,.01,length.out = 300),rep(.01,200)))

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
        axis.title.y = element_blank())
a4
b4=ggplot(data = panB,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
b4
c4=ggplot(data = panC,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
c4
d4=ggplot(data = panD,aes(x=Time,y=Abund,color=sp))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.y = element_blank())
d4

fig4=ggarrange(a4,b4,c4,d4,labels = c("A","B","C","D"),common.legend = T,legend = "top",label.x = .8,label.y = .9)
annotate_figure(fig4,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Time"))

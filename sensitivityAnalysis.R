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

#### JUVENILE SURVIVAL LOOP ####
#looking to see how variation in the survival rate for species 1 changes whether or not stable states occur
times=1:300

ss=seq(0.05,.9,length.out = 10)
store=data.frame(qEs=rep(seq(0,8,length.out=30),10),A1=0,A2=0,J1=0,J2=0,ss=rep(ss,10))
store2=data.frame(qEs=rep(seq(0,8,length.out=30),10),A1=0,A2=0,J1=0,J2=0,ss=rep(ss,10))

y0=c(10,100,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=store$ss[i],m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
  #store$ss[i]=ss[f]
}
y0=c(100,10,40,40)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(1.8,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=store2$ss[i],m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
  #store2$ss[i]=ss[f]
}


panA=as.data.frame(rbind(cbind(store$qEs,store$A1,store$ss,rep("A1",nrow(store))),cbind(store$qEs,store$A2,store$ss,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","ss","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)#;panA$ss=as.numeric(panA$ss)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,store2$ss,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,store2$ss,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","ss","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)#;panB$ss=as.numeric(panB$ss)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp,linetype=ss))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp,linetype=ss))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), common.legend = T,legend="top",label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Harvest Rate (qE)"))


#### CANNIBALISM LOOP ####
#looking to see how variation in the cannibalism rate for species 1 changes whether or not stable states occur
times=1:300

can=seq(0.0001,.05,length.out = 10)
store=data.frame(qEs=rep(seq(0,8,length.out=30),10),A1=0,A2=0,J1=0,J2=0,can=rep(can,300))
store2=data.frame(qEs=rep(seq(0,8,length.out=30),10),A1=0,A2=0,J1=0,J2=0,can=rep(can,300))


y0=c(10,100,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=store$can[i],cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y0=c(100,10,40,40)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(1.8,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=store2$can[i],cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
  
}

panA=as.data.frame(rbind(cbind(store$qEs,store$A1,store$can,rep("A1",nrow(store))),cbind(store$qEs,store$A2,store$can,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","can","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,store2$can,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,store2$can,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","can","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp,linetype=can))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp,linetype=can))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), common.legend = T,legend="top",label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Harvest Rate (qE)"))

#### ADULT NATURAL MORTALITY LOOP ####
#looking to see how variation in the adult natural mortality rate for species 1 changes whether or not stable states occur
times=1:300

ms=seq(0.05,.9,length.out = 10)
store=data.frame(qEs=rep(seq(0,8,length.out=30),10),A1=0,A2=0,J1=0,J2=0,ms=rep(ms,10))
store2=data.frame(qEs=rep(seq(0,8,length.out=30),10),A1=0,A2=0,J1=0,J2=0,ms=rep(ms,10))

y0=c(10,100,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.5,m1=store$ms[i],cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  store$A1[i]=sim[nrow(sim)-1,2]
  store$A2[i]=sim[nrow(sim)-1,3]
  store$J1[i]=sim[nrow(sim)-1,4]
  store$J2[i]=sim[nrow(sim)-1,5]
}
y0=c(100,10,40,40)
for(i in 1:nrow(store2)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store2$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(1.8,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.5,m1=store2$ms[i],cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


panA=as.data.frame(rbind(cbind(store$qEs,store$A1,store$ms,rep("A1",nrow(store))),cbind(store$qEs,store$A2,store$ms,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","ms","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,store2$ms,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,store2$ms,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","ms","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp,linetype=ms))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp,linetype=ms))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), common.legend = T,legend="top",label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Harvest Rate (qE)"))
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
tstep=times=1:300

ss=seq(0.05,.9,length.out = 5)
qEs=seq(0,8,length.out=30)
combos=expand.grid(qEs,ss);colnames(combos)=c("qEs","ss")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

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
  scale_linetype_discrete(name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp,linetype=ss))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  scale_linetype_discrete(name="")+
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

can=seq(0.0001,.05,length.out = 5)
qEs=seq(0,8,length.out=30)
combos=expand.grid(qEs,can);colnames(combos)=c("qEs","can")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)


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

ms=seq(0.05,.9,length.out = 5)
qEs=seq(0,8,length.out=30)
combos=expand.grid(qEs,ms);colnames(combos)=c("qEs","ms")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

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

#### ADULT PREDATION ON J2 LOOP ####
#looking to see how variation in the adult predation j2 for species 1 changes whether or not stable states occur
times=1:300

cj2a1s=seq(0,1,length.out = 5)
qEs=seq(0,8,length.out=30)
combos=expand.grid(qEs,cj2a1s);colnames(combos)=c("qEs","cj2a1s")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y0=c(10,100,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=store$cj2a1s[i],cJ2J1=0.003,v2=1))
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
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=store2$cj2a1s[i],cJ2J1=0.003,v2=1))
  sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


panA=as.data.frame(rbind(cbind(store$qEs,store$A1,store$cj2a1s,rep("A1",nrow(store))),cbind(store$qEs,store$A2,store$cj2a1s,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","cj2a1s","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,store2$cj2a1s,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,store2$cj2a1s,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","cj2a1s","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp,linetype=cj2a1s))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp,linetype=cj2a1s))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), common.legend = T,legend="top",label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Harvest Rate (qE)"))

#### J1 EFFECT ON J2 LOOP ####
#looking to see how variation in the j1 effect on j2 changes whether or not stable states occur
times=1:300

cj2j1s=seq(0,0.5,length.out = 5)
qEs=seq(0,8,length.out=30)
combos=expand.grid(qEs,cj2a1s);colnames(combos)=c("qEs","cj2j1s")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y0=c(10,100,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=store$cj2j1s[i],v2=1))
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
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=store2$cj2j1s[i],v2=1))
  sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


panA=as.data.frame(rbind(cbind(store$qEs,store$A1,store$cj2j1s,rep("A1",nrow(store))),cbind(store$qEs,store$A2,store$cj2j1s,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","cj2j1s","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,store2$cj2j1s,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,store2$cj2j1s,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","cj2j1s","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp,linetype=cj2j1s))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp,linetype=cj2j1s))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), common.legend = T,legend="top",label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Harvest Rate (qE)"))

#### A2 EFFECT ON J1 LOOP ####
#looking to see how variation in the adult predation j1 for species 2 changes whether or not stable states occur
times=1:300

cj1a2s=seq(0,1,length.out = 5)
qEs=seq(0,8,length.out=30)
combos=expand.grid(qEs,cj2a1s);colnames(combos)=c("qEs","cj1a2s")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y0=c(10,100,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=store$cj1a2s[i],cJ1J2=0.003,v1=1),
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
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=store2$cj1a2s[i],cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


panA=as.data.frame(rbind(cbind(store$qEs,store$A1,store$cj1a2s,rep("A1",nrow(store))),cbind(store$qEs,store$A2,store$cj1a2s,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","cj1a2s","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,store2$cj1a2s,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,store2$cj1a2s,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","cj1a2s","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp,linetype=cj1a2s))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp,linetype=cj1a2s))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), common.legend = T,legend="top",label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Harvest Rate (qE)"))


#### RICKER SENSITIVITY ####
#varying the ricker parms a & b for species 1 to see the effect of differing recruitment on stable states
# I need to create a few functions that have different parms since the function has them coded into it.
# based on Hilborn and Walters book ricker parm a is maybe more important to vary
#species 2 stock recruitment will stay the same through all variants of speces 1 recruitment

#original a=5000; b=500
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

#var1 a=2500; b=500
var1<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(2500*A1/(500+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var2 a=1250; b=500
var2<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(1250*A1/(500+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 3 a=7500; b=500
var3<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(7500*A1/(500+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 4 a=4500; b=500
var4<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(4500*A1/(500+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 5 a=5500; b=500
var5<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(5500*A1/(500+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

times=1:300
funlist=list(simBiggsQ2,var1,var2,var3,var4,var5)
as=c(2500,1250,7500,4500,5500,5000)
qEs=seq(0,8,length.out=30)
combos=expand.grid(qEs,as);colnames(combos)=c("qEs","as")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y0=c(10,100,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  if(store$as[i]==5000){
    sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  }else{if(store$as[i]==5500){
    sim=ode(y=y0,times=times,func=var5,parms=p)
  }else{if(store$as[i]==4500){
    sim=ode(y=y0,times=times,func=var4,parms=p)
  }else{if(store$as[i]==7500){
    sim=ode(y=y0,times=times,func=var3,parms=p)
  }else{if(store$as[i]==1250){
    sim=ode(y=y0,times=times,func=var2,parms=p)
  }else{if(store$as[i]==2500){
    sim=ode(y=y0,times=times,func=var1,parms=p)
  }
    }}}}}
  
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
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  if(store2$as[i]==5000){
    sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  }else{if(store2$as[i]==5500){
    sim=ode(y=y0,times=times,func=var5,parms=p)
  }else{if(store2$as[i]==4500){
    sim=ode(y=y0,times=times,func=var4,parms=p)
  }else{if(store2$as[i]==7500){
    sim=ode(y=y0,times=times,func=var3,parms=p)
  }else{if(store2$as[i]==1250){
    sim=ode(y=y0,times=times,func=var2,parms=p)
  }else{if(store2$as[i]==2500){
    sim=ode(y=y0,times=times,func=var1,parms=p)
  }
  }}}}}
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


panA=as.data.frame(rbind(cbind(store$qEs,store$A1,store$as,rep("A1",nrow(store))),cbind(store$qEs,store$A2,store$as,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","as","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,store2$as,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,store2$as,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","as","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp,linetype=as))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp,linetype=as))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), common.legend = T,legend="top",label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Harvest Rate (qE)"))


#### RICKER PARM B SENSITIVITY ####
#original a=5000; b=500
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

#var1 a=5000; b=250
var1b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(5000*A1/(250+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var2 a=5000; b=125
var2b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(1250*A1/(125+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 3 a=5000; b=750
var3b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(5000*A1/(750+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 4 a=5000; b=450
var4b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(5000*A1/(450+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

# var 5 a=5000; b=550
var5b<-function(t,y,params){
  A1<-y[1]
  A2<-y[2]
  J1<-y[3]
  J2<-y[4]
  with(as.list(params),{
    dA1dt=-qE1Fun(t)*A1-m1*A1+s1*J1
    dA2dt=-qE2Fun(t)*A2-m2*A2+s2*J2
    dJ1dt=-cJ1J2*J2*J1-((cJ1A2*v1*A2*J1)/(h1Fun(t)+v1+cJ1A2*A2))-((cJ1A1*v1*A1*J1)/(h1Fun(t)+v1+cJ1A1*A1))-s1*J1+(5000*A1/(550+A1))+st1Fun(t)
    dJ2dt=-cJ2J1*J1*J2-((cJ2A1*v2*A1*J2)/(h2Fun(t)+v2+cJ2A1*A1))-((cJ2A2*v2*A2*J2)/(h2Fun(t)+v2+cJ2A2*A2))-s2*J2+(5000*A2/(500+A2))+st2Fun(t)
    return(list(c(dA1dt,dA2dt,dJ1dt,dJ2dt)))
  })
}

times=1:300
bs=c(250,125,750,450,550,500)
qEs=seq(0,8,length.out=30)
combos=expand.grid(qEs,bs);colnames(combos)=c("qEs","bs")
store=data.frame(A1=0,A2=0,J1=0,J2=0);store=cbind(combos,store)
store2=data.frame(A1=0,A2=0,J1=0,J2=0);store2=cbind(combos,store2)

y0=c(10,100,0,0)
for(i in 1:nrow(store)){
  qE1Fun=approxfun(x=tstep,y=c(rep(store$qEs[i],length(tstep))))
  qE2Fun=approxfun(x=tstep,y=rep(2,length(tstep)))
  h1Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  h2Fun=approxfun(x=tstep,y=rep(8,length(tstep)))
  st1Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  st2Fun=approxfun(x=tstep,y=rep(0,length(tstep)))
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  if(store$bs[i]==500){
    sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  }else{if(store$bs[i]==550){
    sim=ode(y=y0,times=times,func=var5b,parms=p)
  }else{if(store$bs[i]==450){
    sim=ode(y=y0,times=times,func=var4b,parms=p)
  }else{if(store$bs[i]==750){
    sim=ode(y=y0,times=times,func=var3b,parms=p)
  }else{if(store$bs[i]==125){
    sim=ode(y=y0,times=times,func=var2b,parms=p)
  }else{if(store$bs[i]==250){
    sim=ode(y=y0,times=times,func=var1b,parms=p)
  }
  }}}}}
  
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
  
  p=c(c(s1=0.5,m1=0.5,cJ1A1=0.001,cJ1A2=0.5,cJ1J2=0.003,v1=1),
      c(s2=0.5,m2=0.5,cJ2A2=0.001,cJ2A1=0.3,cJ2J1=0.003,v2=1))
  if(store2$bs[i]==500){
    sim=ode(y=y0,times=times,func=simBiggsQ2,parms=p)
  }else{if(store2$bs[i]==550){
    sim=ode(y=y0,times=times,func=var5b,parms=p)
  }else{if(store2$bs[i]==450){
    sim=ode(y=y0,times=times,func=var4b,parms=p)
  }else{if(store2$bs[i]==750){
    sim=ode(y=y0,times=times,func=var3b,parms=p)
  }else{if(store2$bs[i]==125){
    sim=ode(y=y0,times=times,func=var2b,parms=p)
  }else{if(store2$bs[i]==250){
    sim=ode(y=y0,times=times,func=var1b,parms=p)
  }
  }}}}}
  store2$A1[i]=sim[nrow(sim)-1,2]
  store2$A2[i]=sim[nrow(sim)-1,3]
  store2$J1[i]=sim[nrow(sim)-1,4]
  store2$J2[i]=sim[nrow(sim)-1,5]
}


panA=as.data.frame(rbind(cbind(store$qEs,store$A1,store$bs,rep("A1",nrow(store))),cbind(store$qEs,store$A2,store$bs,rep("A2",nrow(store)))))
colnames(panA)=c("qEs","Abund","bs","sp")
panA$qEs=as.numeric(panA$qEs);panA$Abund=as.numeric(panA$Abund)
panB=as.data.frame(rbind(cbind(store2$qEs,store2$A1,store2$bs,rep("A1",nrow(store2))),cbind(store2$qEs,store2$A2,store2$bs,rep("A2",nrow(store2)))))
colnames(panB)=c("qEs","Abund","bs","sp")
panB$qEs=as.numeric(panB$qEs);panB$Abund=as.numeric(panB$Abund)

a=ggplot(data = panA,aes(x=qEs,y=Abund,color=sp,linetype=bs))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.y = element_blank())
a
b=ggplot(data = panB,aes(x=qEs,y=Abund,color=sp,linetype=bs))+theme_classic()+
  geom_line(size=1)+scale_color_manual(values = c("black","grey"),name="")+
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank())
b
fig1=ggarrange(a,b,ncol=1,labels = c("A","B"), common.legend = T,legend="top",label.x = 0.9)
annotate_figure(fig1,
                left = text_grob("Adult Abundance", rot = 90),
                bottom = text_grob("Harvest Rate (qE)"))



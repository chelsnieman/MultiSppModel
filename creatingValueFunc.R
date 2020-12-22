## CD 12.22.2020
## working to put some values on fishery outcomes

## still need the bse q2 function

## consider making another function to add values to the fishery each year of the matrix after running the model
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

#value function

value=function(sim,harv1=qE1Fun(tstep),harv2=qE2Fun(tstep)){
  x=as.data.frame(sim);colnames(x)=c("time","A1","A2","J1","J2")
  x$A1_val=(x$A1*harv1*.9)+(x$A1*.1)
  x$A2_val=(x$A2*harv2*0.25)+(x$A2*0.75)
  
  x$totVal=rowSums(x[,6:7],na.rm = T)
  x=x[is.na(x$A1)==F,] #just to get rid of the last row of NAs that all sims have
  return(x)
}

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
val=value(sim)
plot(sim[,1],sim[,2],type = 'l',lwd=2,ylim = c(0,max(sim[,2:3],na.rm = T)),xlab = "time",ylab="adult abund.")
legend('topleft',legend = c('walleye','bass'),col = c('black','grey'),lwd = 2)
lines(sim[,1],sim[,3],col='gray',lwd=2)
plot(val$time,val$A1_val,type = 'l',lwd=2,ylim=c(0,max(val[,c(6:8)])),xlab = 'time',ylab = 'value')
lines(val$time,val$A2_val,col='grey',lwd=2)
lines(val$time,val$totVal,col='blue',lwd=2)
legend('top',legend = c('sp1','sp2','total'),col = c('black','grey','blue'),lwd=2)
